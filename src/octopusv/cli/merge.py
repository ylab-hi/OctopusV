from pathlib import Path

from typing import List

import typer

from octopusv.merger.base import SVMerger

from octopusv.merger.merge_strategies import UnionStrategy, IntersectionStrategy, SpecificMergeStrategy

def merge(
    input_files: List[Path] = typer.Argument(..., help="List of input SVCF files to merge."),
    output_file: Path = typer.Option(..., help="Output file for merged SV data."),
    intersect: bool = typer.Option(False, "--intersect", help="Apply intersection strategy for merging."),
    union: bool = typer.Option(False, "--union", help="Apply union strategy for merging."),
    specific: List[Path] = typer.Option(None, "--specific", help="Extract SVs that are specifically supported by provided files."),
    overlap: int = typer.Option(None, "--overlap", help="Minimum number of files that must support an SV to be included in the output."),
):
    merger = SVMerger()
    merger.load_data([str(file) for file in input_files])  # Convert Path objects to strings

    if intersect:
        strategy = IntersectionStrategy()
    elif union:
        strategy = UnionStrategy()
    elif specific:
        strategy = SpecificMergeStrategy(specific)
    else:
        strategy = SomeOverlapStrategy(overlap)  # This is hypothetical; you need to implement this

    merged_tree = merger.merge_by_type("any", strategy)  # Replace "any" with actual SV type
    merger.write_result_by_type(output_file, "any")  # Replace "any" with actual SV type or handle dynamically
