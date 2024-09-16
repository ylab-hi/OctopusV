from .sv_interval_tree import SVIntervalTree

# SVMerger: A class for merging and managing structural variant (SV) events
class SVMerger:
    def __init__(self, classified_events):
        # Store the classified events and initialize the interval tree
        self.classified_events = classified_events
        self.interval_tree = SVIntervalTree()

    def merge(self):
        # Merge all classified events into the interval tree
        for sv_type, chromosomes in self.classified_events.items():
            for chromosome, events in chromosomes.items():
                for event in events:
                    # Add each event to the interval tree
                    self.interval_tree.add_event(
                        sv_type, chromosome, event.start_pos, event.end_pos, event.source_file
                    )
        # Merge overlapping intervals in the tree
        self.interval_tree.merge_overlaps()

    def get_events(self, sv_type, chromosome, start, end):
        # Query the interval tree for events within a specific range
        return self.interval_tree.query(sv_type, chromosome, start, end)

    def get_events_by_source(self, sources, operation='union'):
        # Retrieve events based on their source files and specified operation
        return self.interval_tree.get_events_by_source(sources, operation)

    def get_events_by_overlap(self, min_overlap):
        # Get events that are supported by at least a minimum number of sources
        return self.interval_tree.get_events_by_overlap(min_overlap)

    def write_results(self, output_file, events):
        # Write merged events to an output file
        with open(output_file, 'w') as f:
            for sv_type, chromosome, start, end, sources in events:
                # Format each event as a tab-separated line
                f.write(f"{sv_type}\t{chromosome}\t{start}\t{end}\t{','.join(sources)}\n")