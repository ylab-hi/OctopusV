from .sv_interval_tree import SVIntervalTree
from .tra_merger import TRAMerger

# SVMerger: A class for merging and managing structural variant (SV) events
class SVMerger:
    def __init__(self, classified_events):
        # Store the classified events and initialize the interval tree
        self.classified_events = classified_events
        self.interval_tree = SVIntervalTree()
        self.tra_merger = TRAMerger()

    def merge(self):
        for sv_type, chromosomes in self.classified_events.items():
            if sv_type == 'TRA':
                for (chr1, chr2), events in chromosomes.items():
                    for event in events:
                        self.tra_merger.add_event(chr1, event.start_pos, chr2, event.end_pos, event.source_file)
            else:
                for chromosome, events in chromosomes.items():
                    # Merge all classified events (except TRA) into the interval tree
                    for event in events:
                        # Add each event to the interval tree
                        self.interval_tree.add_event(
                            sv_type, chromosome, event.start_pos, event.end_pos, event.source_file
                        )
        # Merge overlapping intervals in the tree
        self.interval_tree.merge_overlaps()

    def get_events(self, sv_type, chromosome, start, end):
        if sv_type == 'TRA':
            return self.tra_merger.get_merged_events()
        # Query the interval tree for events within a specific range
        return self.interval_tree.query(sv_type, chromosome, start, end)

    def get_events_by_source(self, sources, operation='union'):
        # Retrieve events based on their source files and specified operation
        tra_events = self.tra_merger.get_merged_events()
        other_events = self.interval_tree.get_events_by_source(sources, operation)

        if operation == 'union':
            tra_filtered = [event for event in tra_events if any(source in event[5] for source in sources)]
        elif operation == 'intersection':
            tra_filtered = [event for event in tra_events if all(source in event[5] for source in sources)]
        elif operation == 'specific':
            tra_filtered = [event for event in tra_events if set(event[5]) == set(sources)]
        else:
            raise ValueError(f"Unsupported operation: {operation}")

        return other_events + tra_filtered

    def get_events_by_overlap(self, min_overlap):
        # Get events that are supported by at least a minimum number of sources
        tra_events = self.tra_merger.get_merged_events()
        other_events = self.interval_tree.get_events_by_overlap(min_overlap)
        tra_filtered = [event for event in tra_events if len(event[5]) >= min_overlap]
        return other_events + tra_filtered

    def write_results(self, output_file, events):
        # Write merged events to an output file
        with open(output_file, 'w') as f:
            for event in events:
                if event[0] == 'TRA':
                    f.write(f"{event[0]}\t{event[1]}\t{event[2]}\t{event[3]}\t{event[4]}\t{','.join(event[5])}\n")
                else:
                    f.write(f"{event[0]}\t{event[1]}\t{event[2]}\t{event[3]}\t{','.join(event[4])}\n")