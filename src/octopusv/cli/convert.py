from pathlib import Path

import typer


def convert(
    input: Path = typer.Argument(..., exists=True, dir_okay=False, resolve_path=True),
    output: Path = typer.Argument(..., dir_okay=False, resolve_path=True),
    pos_tolerance: int = typer.Option(
        3,
        "--pos-tolerance",
        "-pt",
        help="Position tolerance.",
    ),
):
    """Correct SV events."""
