"""Data identities, source claims, and explicit coordinate conventions."""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from functools import cached_property

from .curation import check_filter
from .errors import CoordinateError


@dataclass(frozen=True)
class Region:
    """Zero-based, half-open interval on a named assembly, or with assembly None,
    on whatever genome build the alignment it's read from uses (its header says).

    reference_length optionally pins the contig length, and is mandatory for
    GRCh37 mitochondrial extraction because hg19 and hs37d5 differ there.
    """

    contig: str
    start: int
    end: int
    assembly: str
    reference_length: int | None = None

    def __post_init__(self):
        if (not self.contig or self.contig in ("*", ".") or re.search(r"[\s:]", self.contig)
                or not isinstance(self.start, int) or not isinstance(self.end, int)
                or self.start < 0 or self.end <= self.start or self.assembly == ""):
            raise CoordinateError(f"Invalid zero-based, half-open region: {self}")
        if self.reference_length is not None and self.end > self.reference_length:
            raise CoordinateError("Region extends beyond its reference length")

    @classmethod
    def from_samtools(cls, text, *, assembly, reference_length=None):
        """Parse contig:start-end with one-based, inclusive coordinates."""
        match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", text)
        if match is None:
            raise CoordinateError("Expected contig:start-end (one-based inclusive)")
        return cls(match[1], int(match[2]) - 1, int(match[3]), assembly, reference_length)


@dataclass(frozen=True)
class SampleClaim:
    """Published sample/library metadata, not an assertion of biological independence."""

    source: str
    label: str = ""
    timepoint: str | None = None
    date: str | None = None
    assay: str | None = None
    platform: str | None = None
    tissue: str | None = None
    provider: str | None = None
    library: str | None = None
    basis: str = "published"


@dataclass(frozen=True)
class File:
    """One file in the bucket, or one of the site's tables. Each processing
    product is a separate file with its own ID."""

    id: str
    key: str
    url: str
    kind: str
    format: str
    size: int | None = None
    modified: str | int | None = None
    index_urls: tuple[str, ...] = ()
    claims: tuple[SampleClaim, ...] = ()
    metadata: dict = field(default_factory=dict, compare=False, repr=False)

    @property
    def samples(self):
        """IDs of the samples this file belongs to (see Dataset.samples)."""
        return self.metadata.get("samples", ())

    def values(self, name, *, include_inferred=False):
        """Distinct nonmissing source assertions for a sample field.

        A less precise date (2025-04) consistent with a more precise one
        (2025-04-09) is not a separate value; claims retain both.
        """
        if name not in SampleClaim.__dataclass_fields__:
            raise ValueError(f"Unknown sample field: {name}")
        values = {getattr(c, name) for c in self.claims
                  if (include_inferred or c.basis == "published") and getattr(c, name)}
        if name == "date":
            values = {v for v in values if not any(o.startswith(v + "-") for o in values)}
        return tuple(sorted(values))

    def resolved(self, name):
        """Return a unique published assertion, otherwise None."""
        values = self.values(name)
        return values[0] if len(values) == 1 else None

    @property
    def conflicts(self):
        return {name: self.values(name) for name in
                ("timepoint", "date", "assay", "platform", "tissue", "provider", "library")
                if len(self.values(name)) > 1}


@dataclass(frozen=True)
class Variant:
    """Website entry plus original allele candidates and source annotations."""

    id: str
    gene: str
    assembly: str
    alleles: tuple[tuple[str, int, str, str], ...]
    status: str
    on_site: bool = True
    vaccine_count: int | None = None
    vaccines: tuple[str, ...] = ()
    pipelines: tuple[str, ...] = ()
    annotations: dict = field(default_factory=dict, compare=False, repr=False)

    @property
    def allele(self):
        if self.status != "ready":
            raise CoordinateError(f"{self.id}: {self.status}")
        return self.alleles[0]

    def region(self, *, padding=0):
        """Span the original VCF reference allele, retaining its anchor."""
        if not isinstance(padding, int) or padding < 0:
            raise ValueError("padding must be a nonnegative integer")
        chrom, pos, ref, _ = self.allele
        return Region(chrom, max(0, pos - 1 - padding), pos - 1 + len(ref) + padding,
                      self.assembly)


class Collection(Sequence):
    """Immutable sequence with explicit, composable equality filters."""

    def __init__(self, items=(), *, source=None):
        self._items = tuple(items)
        self.source = dict(source or {})

    def __len__(self):
        return len(self._items)

    def __iter__(self):
        # Sequence's own __iter__ calls __getitem__ once per item, in Python.
        return iter(self._items)

    def __getitem__(self, key):
        if isinstance(key, slice):
            return type(self)(self._items[key], source=self.source)
        if isinstance(key, str):
            matches = [item for item in self if item.id == key]
            if len(matches) != 1:
                raise KeyError(key)
            return matches[0]
        return self._items[key]

    def where(self, predicate):
        return type(self)((item for item in self if predicate(item)), source=self.source)

    def to_records(self):
        return [asdict(item) for item in self]

    def __getstate__(self):  # indexes are built again, never saved
        return {"_items": self._items, "source": self.source}

    def __repr__(self):
        return f"<{type(self).__name__}: {len(self)} items>"


class Files(Collection):
    """Files in the bucket and the site's tables. Index by key, URL or ID."""

    def __repr__(self):
        from .display import preview
        from .views import size_text
        return preview(f"{len(self):,} files. Narrow them with .select(kind=, sample=, prefix=, ...).",
                       self, ("kind", "size", "key"),
                       lambda f: dict(kind=f.kind, size=size_text(f.size), key=f.key), fixed=("key",))

    def __getitem__(self, key):
        if isinstance(key, str):
            matches = self._names.get(key, resource=False)
            if len(matches) != 1:
                raise KeyError(key)
            return matches[0]
        return super().__getitem__(key)

    @cached_property
    def _names(self):
        return _NameIndex(self)

    def select(self, *, kind=None, format=None, prefix=None, contains=None, sample=None, timepoint=None,
               assay=None, platform=None, tissue=None, provider=None, library=None,
               include_conflicts=False, include_inferred=False):
        """Filter files; unresolved metadata does not match by default.

        sample keeps a sample's BAMs and the files in its FASTQ folders (see
        Dataset.samples). include_conflicts matches any of several published
        assertions; include_inferred also admits explicitly marked path
        inferences. String matching is case-sensitive and uses normalized
        assay names.
        """
        filters = dict(timepoint=timepoint, assay=assay, platform=platform,
                       tissue=tissue, provider=provider, library=library)
        for name in ("assay", "platform", "tissue"):
            check_filter(name, filters[name], lambda name=name: {
                v for f in self for v in f.values(name, include_inferred=True)})

        def match(file):
            if any(value is not None and getattr(file, name) != value
                   for name, value in (("kind", kind), ("format", format))):
                return False
            if prefix is not None and not file.key.startswith(prefix):
                return False
            if contains is not None and contains not in file.key:
                return False
            for name, value in filters.items():
                if value is None:
                    continue
                values = file.values(name, include_inferred=include_inferred)
                if value not in values or (len(values) > 1 and not include_conflicts):
                    return False
            return True
        pool = self if sample is None else self._by_sample.get(sample, ()) if isinstance(sample, str) else ()
        return type(self)((file for file in pool if match(file)), source=self.source)

    @cached_property
    def _by_sample(self):
        by_sample = defaultdict(list)
        for file in self:
            for sample in dict.fromkeys(file.samples):
                by_sample[sample].append(file)
        return by_sample


class _NameIndex:
    """Files by ID, key, URL or resource name. Each kind of name gets its own
    dictionary when first looked up: building one from the files takes about a
    tenth of a second, less than loading a saved one."""

    def __init__(self, files):
        self._files, self._by = files, {}

    def _names(self, kind):
        """{name: file}, with a list for a name several files share."""
        if kind not in self._by:
            files = list(self._files)
            names = [f.metadata.get("resource") if kind == "resource" else getattr(f, kind) for f in files]
            by = dict(zip(names, files))  # built in C: a tenth of the time of a loop
            by.pop(None, None), by.pop("", None)
            if len(by) < sum(1 for name in names if name):
                for name, count in Counter(name for name in names if name).items():
                    if count > 1:
                        by[name] = [f for n, f in zip(names, files) if n == name]
            self._by[kind] = by
        return self._by[kind]

    def find(self, kind, name):
        """The files with this name of this kind (id, key, url or resource)."""
        found = self._names(kind).get(name)
        return [] if found is None else found if isinstance(found, list) else [found]

    def get(self, name, *, resource=True):
        """Every file with this ID, key, URL or (with resource) resource name."""
        if not isinstance(name, str):
            return []
        kinds = ["id", "key", *(["resource"] if resource else []), *(["url"] if "://" in name else [])]
        return list({id(file): file for kind in kinds for file in self.find(kind, name)}.values())


@dataclass(frozen=True)
class Sample:
    """One collected sample: a tumor biopsy or resection, an organoid, or a blood draw.

    sequencing lists (assay, platform) pairs from the site's sample registry
    and from the sample's FASTQ folders. bams are file keys; fastq_folders are
    bucket folders. missing_bams are BAMs the registry names but the bucket
    doesn't have.
    """

    id: str
    timepoint: str | None = None
    date: str | None = None
    tissue: str | None = None
    site: str | None = None
    description: str | None = None
    providers: tuple[str, ...] = ()
    sequencing: tuple[tuple[str, str | None], ...] = ()
    bams: tuple[str, ...] = ()
    fastq_folders: tuple[str, ...] = ()
    missing_bams: tuple[str, ...] = ()
    notes: str = ""
    disagreements: tuple[dict, ...] = field(default=(), compare=False)
    corrections: tuple[str, ...] = ()
    details: dict = field(default_factory=dict, compare=False, repr=False)

    @property
    def assays(self):
        """Distinct assay names, such as ('rna-seq', 'wgs', 'scrna-seq')."""
        return tuple(dict.fromkeys(assay for assay, _ in self.sequencing))

    @property
    def files(self):
        """This sample's BAMs and the files in its FASTQ folders."""
        dataset = getattr(self, "_dataset", None)
        if dataset is None:
            raise ValueError(f"{self.id} is not attached to a Dataset")
        return dataset.files.select(sample=self.id)

    def __repr__(self):
        from .display import Text
        from .views import sample_view
        return str(Text(sample_view(self, getattr(self, "_dataset", None), python=True)))


class Samples(Collection):
    """Samples in collection order. Index by ID: data.samples["T1_tumor"]."""

    def __getitem__(self, key):
        if isinstance(key, str):
            matches = [s for s in self if s.id == key]
            if not matches:
                raise KeyError(f"No sample {key!r}; samples: {', '.join(s.id for s in self)}")
            return matches[0]
        return super().__getitem__(key)

    def select(self, *, timepoint=None, tissue=None, assay=None, platform=None):
        """Samples with this timepoint, tissue, or sequencing (assay and platform
        use the file filter names: rna-seq, scrna-seq, ont, ...)."""
        for name, value in (("tissue", tissue), ("assay", assay), ("platform", platform)):
            check_filter(name, value)
        return self.where(lambda s: (timepoint is None or s.timepoint == timepoint)
                          and (tissue is None or s.tissue == tissue)
                          and (assay is None and platform is None or any(
                              (assay is None or a == assay) and (platform is None or p == platform)
                              for a, p in s.sequencing)))

    def to_records(self):
        return [asdict(sample) for sample in self]

    def __repr__(self):
        from .display import Text
        from .views import samples_view
        return str(Text(samples_view(self)))


class Variants(Collection):
    def __repr__(self):
        from .display import preview
        from .views import allele_text
        ready = sum(v.status == "ready" for v in self)
        return preview(f"{len(self)} variants, {ready} with a ready allele. Narrow them with "
                       ".select(gene=, status=, vaccine=, pipeline=); index by ID.", self,
                       ("id", "gene", "status", "allele"),
                       lambda v: dict(id=v.id, gene=v.gene, status=v.status, allele=allele_text(v)))

    def select(self, *, gene=None, ids=None, vaccine=None, vaccinated=None, pipeline=None,
               status=None, on_site=None, vaccine_source="overlap"):
        if vaccine_source not in ("overlap", "source_variants"):
            raise ValueError("vaccine_source must be 'overlap' or 'source_variants'")
        ids = None if ids is None else {ids} if isinstance(ids, str) else set(ids)
        return self.where(lambda v:
                          (gene is None or v.gene == gene)
                          and (ids is None or v.id in ids)
                          and (vaccine is None or vaccine in (v.vaccines if vaccine_source == "overlap"
                               else v.annotations.get("source_vaccines", ())))
                          and (vaccinated is None or (v.vaccine_count is not None
                               and (v.vaccine_count > 0) == vaccinated))
                          and (pipeline is None or pipeline in v.pipelines)
                          and (status is None or v.status == status)
                          and (on_site is None or v.on_site == on_site))

    def regions(self, *, padding=0):
        """Convert all selected literal alleles; unresolved entries raise."""
        return tuple(v.region(padding=padding) for v in self)

    def to_varcode(self, *, genome, assembly=None, convert_ucsc_contig_names=True,
                   normalize_contig_names=True):
        """Return native variants on the caller's PyEnsembl reference.

        Custom-named subset genomes require an explicit assembly. Their unique
        reference_name is preserved. Every selected entry must have a ready
        allele; metadata retains original entries even if Varcode merges them.
        No annotation data are downloaded by this adapter.

        Varcode converts primary UCSC names (chr1 -> 1, chrM -> MT) by default
        and retains original_contig. Set both contig options to False to use
        literal names in a custom reference. Renaming is not assembly liftover.
        """
        from varcode import Variant as NativeVariant
        from varcode import VariantCollection

        from .reads import ASSEMBLY_LENGTHS, normalize_assembly
        named = normalize_assembly(genome.reference_name)
        declared = normalize_assembly(assembly) if assembly else named
        if declared not in ASSEMBLY_LENGTHS:
            raise ValueError("Supply assembly='GRCh38' or 'GRCh37' for a custom-named genome")
        if named in ASSEMBLY_LENGTHS and named != declared:
            raise ValueError("Declared assembly conflicts with the genome reference name")
        result, metadata = [], {}
        for item in self:
            if normalize_assembly(item.assembly) != declared:
                raise ValueError("Genome assembly must match variant assembly")
            chrom, pos, ref, alt = item.allele
            variant = NativeVariant(chrom, pos, ref, alt, ensembl=genome,
                                    convert_ucsc_contig_names=convert_ucsc_contig_names,
                                    normalize_contig_names=normalize_contig_names)
            result.append(variant)
            metadata.setdefault(variant, dict(source=dict(self.source), entries=[]))["entries"].append(asdict(item))
        source = "osteosarc:" + self.source.get("snapshot_id", "unversioned")
        return VariantCollection(result, source_to_metadata_dict={source: metadata})
