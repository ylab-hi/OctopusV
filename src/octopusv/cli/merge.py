from pathlib import Path
import typer
from typing import List, Optional
from octopusv.merger.sv_merger import SVMerger
from octopusv.utils.SV_classifier_by_chromosome import SVClassifiedByChromosome
from octopusv.utils.SV_classifier_by_type import SVClassifierByType
from octopusv.utils.svcf_parser import SVCFFileEventCreator


def merge(
        input_files: Optional[List[Path]] = typer.Argument(None, help="List of input SVCF files to merge."),
        input_option: Optional[List[Path]] = typer.Option(None, "--input-file", "-i", help="Input SVCF files to merge."),
        output_file: Path = typer.Option(..., "--output-file", "-o", help="Output file for merged SV data."),
        intersect: bool = typer.Option(False, "--intersect", help="Apply intersection strategy for merging."),
        union: bool = typer.Option(False, "--union", help="Apply union strategy for merging."),
        specific: List[Path] = typer.Option(
            None, "--specific", help="Extract SVs that are specifically supported by provided files."
        ),
        overlap: int = typer.Option(
            None, "--overlap", help="Minimum number of files that must support an SV to be included in the output."
        ),
        tra_distance_threshold: int = typer.Option(
            100, "--tra-threshold", help="Distance threshold for merging TRA events (in base pairs)."
        ),
):
    """
    Merge multiple SVCF files based on specified strategy.
    """
    # Combine input files from arguments and options
    all_input_files = (input_files or []) + (input_option or [])

    if not all_input_files:
        typer.echo("Error: No input files provided.", err=True)
        raise typer.Exit(code=1)
    if specific and not specific[0]:
        typer.echo("Error: --specific option requires at least one file.", err=True)
        raise typer.Exit(code=1)
    if overlap is not None and overlap < 1:
        typer.echo("Error: --overlap must be a positive integer.", err=True)
        raise typer.Exit(code=1)

    sv_event_creator = SVCFFileEventCreator([str(file) for file in all_input_files])
    sv_event_creator.parse()
    classifier = SVClassifierByType(sv_event_creator.events)
    classifier.classify()
    chromosome_classifier = SVClassifiedByChromosome(classifier.get_classified_events())
    chromosome_classifier.classify()
    sv_merger = SVMerger(chromosome_classifier.get_classified_events(), tra_distance_threshold)
    sv_merger.merge()

    if intersect:
        results = sv_merger.get_events_by_source([str(file) for file in all_input_files], operation='intersection')
    elif union:
        results = sv_merger.get_events_by_source([str(file) for file in all_input_files], operation='union')
    elif specific:
        results = sv_merger.get_events_by_source([str(file) for file in specific], operation='specific')
    elif overlap is not None:
        results = sv_merger.get_events_by_overlap(overlap)
    else:
        raise ValueError("No merge strategy specified. Please use --intersect, --union, --specific, or --overlap.")

    sv_merger.write_results(output_file, results)
    typer.echo(f"Merged results written to {output_file}")