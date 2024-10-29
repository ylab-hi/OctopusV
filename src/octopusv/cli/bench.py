from pathlib import Path

import typer
from octopusv.bencher.sv_bencher import SVBencher


def bench(
    truth_file: Path = typer.Argument(..., help="Path to the truth (ground truth) SVCF file."),
    call_file: Path = typer.Argument(..., help="Path to the call (test) SVCF file."),
    output_dir: Path = typer.Option(..., "--output-dir", "-o", help="Output directory for benchmark results."),
    max_distance: int = typer.Option(
        50, "--max-distance", help="Maximum allowed distance between start or end positions for matching events."
    ),
    max_length_ratio: float = typer.Option(
        1.3, "--max-length-ratio", help="Maximum allowed ratio between event lengths for matching events."
    ),
    min_jaccard: float = typer.Option(
        0.7, "--min-jaccard", help="Minimum required Jaccard index for overlap to match events."
    ),
    tra_delta: int = typer.Option(
        50, "--tra-delta", help="Position uncertainty threshold for TRA events (in base pairs)."
    ),
    tra_min_overlap_ratio: float = typer.Option(0.5, "--tra-min-overlap", help="Minimum overlap ratio for TRA events."),
    tra_strand_consistency: bool = typer.Option(
        True, "--tra-strand-consistency", help="Whether to require strand consistency for TRA events."
    ),
):
    """Benchmark structural variation calls against a truth set."""
    bencher = SVBencher(
        truth_file,
        call_file,
        output_dir,
        max_distance=max_distance,
        max_length_ratio=max_length_ratio,
        min_jaccard=min_jaccard,
        tra_delta=tra_delta,
        tra_min_overlap_ratio=tra_min_overlap_ratio,
        tra_strand_consistency=tra_strand_consistency,
    )

    bencher.run_benchmark()
    typer.echo(f"Benchmark results written to {output_dir}")


if __name__ == "__main__":
    typer.run(bench)
