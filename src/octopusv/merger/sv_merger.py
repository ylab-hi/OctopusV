from .tra_merger import TRAMerger
from .sv_merge_logic import should_merge
from typing import List, Dict

class SVMerger:
    def __init__(self, classified_events, tra_distance_threshold=100, max_distance=50, max_length_ratio=1.3, min_jaccard=0.7):
        self.classified_events = classified_events
        self.merged_events: Dict[str, Dict[str, List]] = {}
        self.tra_merger = TRAMerger(tra_distance_threshold)
        self.max_distance = max_distance
        self.max_length_ratio = max_length_ratio
        self.min_jaccard = min_jaccard

    def merge(self):
        for sv_type, chromosomes in self.classified_events.items():
            if sv_type == 'TRA':
                for (chr1, chr2), events in chromosomes.items():
                    for event in events:
                        self.tra_merger.add_event(chr1, event.start_pos, chr2, event.end_pos, event.source_file, event.bnd_pattern)
            else:
                if sv_type not in self.merged_events:
                    self.merged_events[sv_type] = {}
                for chromosome, events in chromosomes.items():
                    if chromosome not in self.merged_events[sv_type]:
                        self.merged_events[sv_type][chromosome] = []
                    for event in events:
                        self.add_and_merge_event(sv_type, chromosome, event)

    def add_and_merge_event(self, sv_type, chromosome, new_event):
        events = self.merged_events[sv_type][chromosome]
        for existing_event in events:
            if should_merge(existing_event, new_event, self.max_distance, self.max_length_ratio, self.min_jaccard):
                self.merge_events(existing_event, new_event)
                return
        events.append(new_event)

    def merge_events(self, event1, event2):
        event1.start_pos = min(event1.start_pos, event2.start_pos)
        event1.end_pos = max(event1.end_pos, event2.end_pos)
        event1.source_file = f"{event1.source_file},{event2.source_file}"

    def get_events(self, sv_type, chromosome, start, end):
        if sv_type == 'TRA':
            return self.tra_merger.get_merged_events()
        if sv_type in self.merged_events and chromosome in self.merged_events[sv_type]:
            return [event for event in self.merged_events[sv_type][chromosome] if event.start_pos <= end and event.end_pos >= start]
        return []

    def get_events_by_source(self, sources, operation='union'):
        tra_events = self.tra_merger.get_merged_events()
        other_events = []
        for sv_type, chromosomes in self.merged_events.items():
            for chromosome, events in chromosomes.items():
                other_events.extend(events)

        if operation == 'union':
            tra_filtered = [event for event in tra_events if any(source in event[5] for source in sources)]
            other_filtered = [event for event in other_events if any(source in event.source_file.split(',') for source in sources)]
        elif operation == 'intersection':
            tra_filtered = [event for event in tra_events if all(source in event[5] for source in sources)]
            other_filtered = [event for event in other_events if all(source in event.source_file.split(',') for source in sources)]
        elif operation == 'specific':
            tra_filtered = [event for event in tra_events if set(event[5]) == set(sources)]
            other_filtered = [event for event in other_events if set(event.source_file.split(',')) == set(sources)]
        else:
            raise ValueError(f"Unsupported operation: {operation}")

        return other_filtered + tra_filtered

    def get_events_by_overlap(self, min_overlap):
        tra_events = self.tra_merger.get_merged_events()
        other_events = []
        for sv_type, chromosomes in self.merged_events.items():
            for chromosome, events in chromosomes.items():
                other_events.extend(events)

        tra_filtered = [event for event in tra_events if len(event[5]) >= min_overlap]
        other_filtered = [event for event in other_events if len(event.source_file.split(',')) >= min_overlap]
        return other_filtered + tra_filtered

    def write_results(self, output_file, events):
        with open(output_file, 'w') as f:
            for event in events:
                if isinstance(event, tuple):  # TRA event
                    f.write(f"{event[0]}\t{event[1]}\t{event[2]}\t{event[3]}\t{event[4]}\t{','.join(event[5])}\t{event[6]}\n")
                else:  # Other SV events
                    f.write(f"{event.sv_type}\t{event.chrom}\t{event.start_pos}\t{event.end_pos}\t{event.source_file}\n")