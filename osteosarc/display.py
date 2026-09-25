"""Plain-text previews for exploring the data in a REPL or notebook."""


class Text(str):
    """A string that displays as text, without quotes or escapes, when echoed."""

    def __repr__(self):
        return str(self)


def preview(title, rows, columns, *, total, limit=8, fixed=(), width=None):
    """A title line and the first few rows as a table, noting how many more there are."""
    from .explore import table
    if not total:
        return f"{title}\n(none)"
    return f"{title}\n" + table(list(rows)[:limit], columns, width=width, fixed=fixed) + (
        f"\n... {total - limit} more" if total > limit else "")
