from pathlib import Path

import typer
from octopusv.merger.sv_merger import SVMerger
from octopusv.utils.SV_classifier_by_chromosome import SVClassifiedByChromosome
from octopusv.utils.SV_classifier_by_type import SVClassifierByType
from octopusv.utils.svcf_parser import SVCFFileEventCreator

def merge(
    input_files: list[Path] = typer.Argument(..., help="List of input SVCF files to merge."),
    output_file: Path = typer.Option(..., help="Output file for merged SV data."),
    intersect: bool = typer.Option(False, "--intersect", help="Apply intersection strategy for merging."),
    union: bool = typer.Option(False, "--union", help="Apply union strategy for merging."),
    specific: list[Path] = typer.Option(
        None, "--specific", help="Extract SVs that are specifically supported by provided files."
    ),
    overlap: int = typer.Option(
        None, "--overlap", help="Minimum number of files that must support an SV to be included in the output."
    ),
):
    if not input_files:
        typer.echo("Error: No input files provided.", err=True)
        raise typer.Exit(code=1)

    if specific and not specific[0]:
        typer.echo("Error: --specific option requires at least one file.", err=True)
        raise typer.Exit(code=1)

    if overlap is not None and overlap < 1:
        typer.echo("Error: --overlap must be a positive integer.", err=True)
        raise typer.Exit(code=1)

    sv_event_creator = SVCFFileEventCreator([str(file) for file in input_files])
    sv_event_creator.parse()

    classifier = SVClassifierByType(sv_event_creator.events)
    classifier.classify()

    chromosome_classifier = SVClassifiedByChromosome(classifier.get_classified_events())
    chromosome_classifier.classify()

    sv_merger = SVMerger(chromosome_classifier.get_classified_events())
    sv_merger.merge()

    if intersect:
        results = sv_merger.get_events_by_source([str(file) for file in input_files], operation='intersection')
    elif union:
        results = sv_merger.get_events_by_source([str(file) for file in input_files], operation='union')
    elif specific:
        results = sv_merger.get_events_by_source([str(file) for file in specific], operation='specific')
    elif overlap is not None:
        results = sv_merger.get_events_by_overlap(overlap)
    else:
        raise ValueError("No merge strategy specified. Please use --intersect, --union, --specific, or --overlap.")

    sv_merger.write_results(output_file, results)

    typer.echo(f"Merged results written to {output_file}")

if __name__ == "__main__":
    typer.run(merge)