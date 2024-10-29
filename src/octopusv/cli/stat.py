from pathlib import Path

import typer
from octopusv.stater.sv_stater import SVStater


def stat(
    input_file: Path = typer.Argument(..., help="Input SVCF file to analyze."),
    output_file: Path = typer.Option(..., "--output-file", "-o", help="Output file for statistics."),
    min_size: int = typer.Option(50, "--min-size", help="Minimum SV size to consider."),
    max_size: int = typer.Option(None, "--max-size", help="Maximum SV size to consider."),
):
    """Analyze a single SVCF file and generate comprehensive statistics."""
    sv_stater = SVStater(str(input_file), min_size=min_size, max_size=max_size)
    sv_stater.analyze()
    sv_stater.write_results(output_file)
    typer.echo(f"Analysis results written to {output_file}")


if __name__ == "__main__":
    typer.run(stat)
