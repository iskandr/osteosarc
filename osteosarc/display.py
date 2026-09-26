"""Plain-text previews for exploring the data in a REPL or notebook."""

PREVIEW_ROWS = 8


class Text(str):
    """A string that displays as text, without quotes or escapes, when echoed."""

    def __repr__(self):
        return str(self)


def preview(title, items, columns, row, *, fixed=()):
    """A title line and the first few items as table rows, noting how many more there are."""
    from .views import table
    if not len(items):
        return Text(f"{title}\n(none)")
    rows = [row(item) for item in items[:PREVIEW_ROWS]]
    return Text(f"{title}\n" + table(rows, columns, limit=PREVIEW_ROWS, total=len(items), fixed=fixed))
