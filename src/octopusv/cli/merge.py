from pathlib import Path
from typing import List
import typer
from octopusv.merger.sv_interval_tree_merger import SVIntervalTreeMerger
from octopusv.merger.merge_strategies import UnionStrategy, IntersectionStrategy, SpecificMergeStrategy
from octopusv.utils.svcf_parser import SVCFFileEventCreator
from octopusv.utils.SV_classifier_by_type import SVClassifierByType

def merge(
    input_files: List[Path] = typer.Argument(..., help="List of input SVCF files to merge."),
    output_file: Path = typer.Option(..., help="Output file for merged SV data."),
    intersect: bool = typer.Option(False, "--intersect", help="Apply intersection strategy for merging."),
    union: bool = typer.Option(False, "--union", help="Apply union strategy for merging."),
    specific: List[Path] = typer.Option(None, "--specific", help="Extract SVs that are specifically supported by provided files."),
    overlap: int = typer.Option(None, "--overlap", help="Minimum number of files that must support an SV to be included in the output."),
):
    sv_event_creator = SVCFFileEventCreator([str(file) for file in input_files]) # Create a lot of SVCF events
    sv_event_creator.parse()  # This function will add all the events to the event attribute.

    # Classify events by type
    classifier = SVClassifierByType(sv_event_creator.events)
    classifier.classify()

    sv_interval_tree_merger = SVIntervalTreeMerger()
    sv_interval_tree_merger.load_events_into_tree(sv_event_creator.events)  # Feed the events into merger and load into interval tree

    if intersect:
        strategy = IntersectionStrategy()
    elif union:
        strategy = UnionStrategy()
    elif specific:
        strategy = SpecificMergeStrategy(specific)
    else:
        strategy = SomeOverlapStrategy(overlap)  # This is hypothetical; you need to implement this

    merged_tree = sv_interval_tree_merger.merge_by_type("any", strategy)  # Replace "any" with actual SV type
    sv_interval_tree_merger.write_result_by_type(output_file, "any")  # Replace "any" with actual SV type or handle dynamically
