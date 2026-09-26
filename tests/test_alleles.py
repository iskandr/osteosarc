"""The allele classifier on hand-built alignments over a known reference."""

import pytest

from osteosarc import CoordinateError
from osteosarc.alleles import allele_window, read_allele, template_allele

pysam = pytest.importorskip("pysam")

#                 0         1         2         3
#                 0123456789012345678901234567890123456789
REFERENCE = "GATTACAGCTAAAAAGCGTCAGTCAGTTGACCATGGCATG"
HEADER = pysam.AlignmentHeader.from_dict(dict(SQ=[dict(SN="chr1", LN=1000)]))


def read(start, cigar, sequence, *, flag=0):
    record = pysam.AlignedSegment(HEADER)
    record.query_name, record.flag, record.reference_id = "r", flag, 0
    record.reference_start, record.cigarstring, record.query_sequence = start, cigar, sequence
    return record


def window(position, ref, alt):
    return allele_window("chr1", position, ref, alt, REFERENCE, 0)


def test_a_snv_is_read_base_by_base():
    snv = window(8, "G", "T")  # zero-based 7
    assert (snv.start, snv.end, snv.ref, snv.alt) == (7, 8, "G", "T")
    assert read_allele(read(2, "10M", REFERENCE[2:12]), snv) == "ref"
    assert read_allele(read(2, "10M", REFERENCE[2:7] + "T" + REFERENCE[8:12]), snv) == "alt"
    assert read_allele(read(2, "10M", REFERENCE[2:7] + "A" + REFERENCE[8:12]), snv) == "other"
    assert read_allele(read(2, "10M", REFERENCE[2:7] + "N" + REFERENCE[8:12]), snv) == "uncallable"
    # The bases either side must be aligned too.
    assert read_allele(read(7, "5M", REFERENCE[7:12]), snv) == "uncallable"
    assert read_allele(read(2, "6M", REFERENCE[2:8]), snv) == "uncallable"
    assert read_allele(read(2, "10M", REFERENCE[2:12], flag=256), snv) == "uncallable"


def test_a_deletion_counts_wherever_the_aligner_put_it_in_the_repeat():
    # Deleting one A from the run AAAAA at zero-based 10..14 (VCF: T>... anchored at 9).
    deletion = window(10, "TA", "T")
    assert (deletion.start, deletion.end) == (10, 15)
    assert (deletion.ref, deletion.alt) == ("AAAAA", "AAAA")
    shortened = REFERENCE[5:10] + "AAAA" + REFERENCE[15:20]
    for cigar in ("5M1D9M", "9M1D5M", "7M1D7M"):  # the same haplotype, deleted at different places
        assert read_allele(read(5, cigar, shortened), deletion) == "alt"
    assert read_allele(read(5, "15M", REFERENCE[5:20]), deletion) == "ref"
    assert read_allele(read(5, "5M2D8M", REFERENCE[5:10] + "AAA" + REFERENCE[15:20]), deletion) == "other"
    # Reads that stop inside the repeat can't tell.
    assert read_allele(read(5, "8M", REFERENCE[5:13]), deletion) == "uncallable"


def test_insertions_splices_and_complex_changes():
    # Inserting CG after zero-based 15 (G), inside the GCG run at 15..17.
    insertion = window(17, "C", "CGC")
    assert (insertion.start, insertion.end, insertion.ref, insertion.alt) == (15, 18, "GCG", "GCGCG")
    inserted = REFERENCE[10:17] + "GC" + REFERENCE[17:25]
    assert read_allele(read(10, "7M2I8M", inserted), insertion) == "alt"
    assert read_allele(read(10, "8M2I7M", inserted), insertion) == "alt"  # the same haplotype, shifted
    assert read_allele(read(10, "15M", REFERENCE[10:25]), insertion) == "ref"
    # A splice across the window hides it.
    spliced = read(2, "4M20N4M", REFERENCE[2:6] + REFERENCE[26:30])
    assert read_allele(spliced, window(8, "G", "T")) == "uncallable"
    # A complex change is compared as a whole, however the aligner wrote it.
    complex_change = window(20, "CAGT", "GG")
    assert (complex_change.start, complex_change.end, complex_change.ref, complex_change.alt) == (19, 23, "CAGT", "GG")
    changed = REFERENCE[15:19] + "GG" + REFERENCE[23:28]
    assert read_allele(read(15, "4M4D2I5M", changed), complex_change) == "alt"
    assert read_allele(read(15, "4M2I4D5M", changed), complex_change) == "alt"
    assert read_allele(read(15, "13M", REFERENCE[15:28]), complex_change) == "ref"
    assert read_allele(read(15, "4M4D1I5M", REFERENCE[15:19] + "G" + REFERENCE[23:28]), complex_change) == "other"


def test_windows_check_the_reference_and_the_context():
    with pytest.raises(CoordinateError, match="doesn't match"):
        window(8, "A", "T")
    with pytest.raises(CoordinateError, match="small variant"):
        window(8, "G", "G")
    with pytest.raises(CoordinateError):
        allele_window("chr1", 11, "TA", "T", REFERENCE[:13], 0)  # the repeat runs past the context


def test_a_template_shows_what_its_alignments_agree_on():
    assert template_allele(["alt", "uncallable"]) == "alt"
    assert template_allele(["ref", "ref"]) == "ref"
    assert template_allele(["alt", "ref"]) == "other"
    assert template_allele(["uncallable"]) == template_allele([]) == "uncallable"
