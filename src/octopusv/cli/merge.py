from pathlib import Path

import typer
from octopusv.merger.merge_strategies import IntersectionStrategy, SpecificMergeStrategy, UnionStrategy
from octopusv.merger.sv_interval_tree_merger import SVIntervalTreeMerger
from octopusv.utils.SV_classified_by_chromosome import SVClassifiedByChromosome
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
    sv_event_creator = SVCFFileEventCreator([str(file) for file in input_files])  # Create a lot of SVCF events
    sv_event_creator.parse()  # This function will add all the events to the event attribute.

    # Classify events by type
    classifier = SVClassifierByType(sv_event_creator.events)
    classifier.classify()

    # Further classify events by chromosome
    chromosome_classifier = SVClassifiedByChromosome(classifier.get_classified_events())
    chromosome_classifier.classify()

    sv_interval_tree_merger = SVIntervalTreeMerger()

    # Load events into tree for each SV type and chromosome
    for sv_type, chrom_events in chromosome_classifier.get_classified_events().items():
        for chrom, events in chrom_events.items():  # like { "chr1": [event1, event2], "chr5": [event3]}
            sv_interval_tree_merger.load_events_into_tree(events, sv_type, chrom)  # SVs with the same type and same chr

    # After all events are loaded, perform the merging of overlaps
    sv_interval_tree_merger.merge_by_type_and_chromosome()

    if intersect:
        strategy = IntersectionStrategy()
    elif union:
        strategy = UnionStrategy()
    elif specific:
        strategy = SpecificMergeStrategy(specific)
    else:
        strategy = SomeOverlapStrategy(overlap)  # This is hypothetical; you need to implement this

    merged_tree = sv_interval_tree_merger.merge_by_type("any", strategy)  # Replace "any" with actual SV type
    sv_interval_tree_merger.write_result_by_type(
        output_file, "any"
    )  # Replace "any" with actual SV type or handle dynamically


"""
# Before merging overlaps:
{
    "INS": {
        "chr1": IntervalTree([
            Interval(100, 200, {'fileA'}),
            Interval(150, 250, {'fileB'}),
            Interval(180, 300, {'fileC'})
        ]),
        "chr5": IntervalTree([
            Interval(300, 400, {'fileC'})
        ])
    },
    "DEL": {
        "chr5": IntervalTree([
            Interval(500, 600, {'fileD'}),
            Interval(450, 550, {'fileE'})
        ]),
        "chr8": IntervalTree([
            Interval(700, 800, {'fileF'})
        ])
    }
}

# After merging overlaps:
{
    "INS": {
        "chr1": IntervalTree([
            Interval(100, 300, {'fileA', 'fileB', 'fileC'})  # Merged all overlapping intervals
        ]),
        "chr5": IntervalTree([
            Interval(300, 400, {'fileC'})  # No overlaps to merge
        ])
    },
    "DEL": {
        "chr5": IntervalTree([
            Interval(450, 600, {'fileD', 'fileE'})  # Merged overlapping intervals
        ]),
        "chr8": IntervalTree([
            Interval(700, 800, {'fileF'})  # No overlaps
        ])
    }
}

"""
