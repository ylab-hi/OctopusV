from .tra_merger import TRAMerger
from .sv_merge_logic import should_merge
from .sv_selector import select_representative_sv
from typing import List, Dict

class SVMerger:
    def __init__(self, classified_events, tra_delta=50, tra_min_overlap_ratio=0.5, tra_strand_consistency=True,
                max_distance=50, max_length_ratio=1.3, min_jaccard=0.7):
        self.classified_events = classified_events
        self.merged_events: Dict[str, Dict[str, List]] = {}
        self.event_groups: Dict[str, Dict[str, List[List]]] = {}
        self.tra_merger = TRAMerger(tra_delta, tra_min_overlap_ratio, tra_strand_consistency)
        self.max_distance = max_distance
        self.max_length_ratio = max_length_ratio
        self.min_jaccard = min_jaccard

    def merge(self):
        for sv_type, chromosomes in self.classified_events.items():
            if sv_type == 'TRA':
                for (chr1, chr2), events in chromosomes.items():
                    for event in events:
                        self.tra_merger.add_event(event)
            else:
                if sv_type not in self.merged_events:
                    self.merged_events[sv_type] = {}
                    self.event_groups[sv_type] = {}
                for chromosome, events in chromosomes.items():
                    if chromosome not in self.merged_events[sv_type]:
                        self.merged_events[sv_type][chromosome] = []
                        self.event_groups[sv_type][chromosome] = []
                    for event in events:
                        self.add_and_merge_event(sv_type, chromosome, event)

    def add_and_merge_event(self, sv_type, chromosome, new_event):
        events = self.merged_events[sv_type][chromosome]
        event_groups = self.event_groups[sv_type][chromosome]
        for idx, existing_event in enumerate(events):
            if should_merge(existing_event, new_event, self.max_distance, self.max_length_ratio, self.min_jaccard):
                # Combind
                event_groups[idx].append(new_event)
                return
        events.append(new_event)
        event_groups.append([new_event])

    def get_events(self, sv_type, chromosome, start, end):
        if sv_type == 'TRA':
            return self.tra_merger.get_merged_events()
        if sv_type in self.event_groups and chromosome in self.event_groups[sv_type]:
            events = []
            for sv_group in self.event_groups[sv_type][chromosome]:
                representative_sv = select_representative_sv(sv_group)
                if representative_sv.start_pos <= end and representative_sv.end_pos >= start:
                    events.append(representative_sv)
            return events
        return []

    def get_all_merged_events(self):
        merged_events = []
        for sv_type in self.event_groups:
            for chromosome in self.event_groups[sv_type]:
                for sv_group in self.event_groups[sv_type][chromosome]:
                    representative_sv = select_representative_sv(sv_group)
                    merged_events.append(representative_sv)
        return merged_events

    def get_events_by_source(self, sources, operation='union'):
        tra_events = self.tra_merger.get_merged_events()
        other_events = self.get_all_merged_events()

        if operation == 'union':
            tra_filtered = [event for event in tra_events if any(source in event.source_file for source in sources)]
            other_filtered = [event for event in other_events if any(source in event.source_file.split(',') for source in sources)]
        elif operation == 'intersection':
            tra_filtered = [event for event in tra_events if all(source in event.source_file for source in sources)]
            other_filtered = [event for event in other_events if all(source in event.source_file.split(',') for source in sources)]
        elif operation == 'specific':
            tra_filtered = [event for event in tra_events if set(event.source_file.split(',')) == set(sources)]
            other_filtered = [event for event in other_events if set(event.source_file.split(',')) == set(sources)]
        else:
            raise ValueError(f"Unsupported operation: {operation}")

        return other_filtered + tra_filtered

    def get_events_by_overlap(self, min_overlap):
        tra_events = self.tra_merger.get_merged_events()
        other_events = self.get_all_merged_events()

        tra_filtered = [event for event in tra_events if len(set(event.source_file.split(','))) >= min_overlap]
        other_filtered = [event for event in other_events if len(set(event.source_file.split(','))) >= min_overlap]
        return other_filtered + tra_filtered

    def write_results(self, output_file, events):
        with open(output_file, 'w') as f:
            # Write VCF header if needed
            # f.write('##fileformat=VCFv4.2\n')
            # f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
            for event in events:
                if isinstance(event, tuple):  # TRA event
                    # Add source_file info
                    f.write(
                        f"{event.chrom}\t{event.pos}\t{event.sv_id}\t{event.ref}\t{event.alt}\t{event.quality}\t{event.filter}\t"
                        f"{';'.join([f'{k}={v}' for k, v in event.info.items()])};SOURCES={event.source_file}\t"
                        f"{event.format}\t{':'.join(event.sample.values())}\n"
                    )
                else:  # Other SV events
                    # Add SOURCES in INFO
                    info_field = ';'.join([f'{k}={v}' for k, v in event.info.items()])
                    info_field += f";SOURCES={event.source_file}"
                    sample_field = ':'.join(event.sample.values())
                    f.write(
                        f"{event.chrom}\t{event.pos}\t{event.sv_id}\t{event.ref}\t{event.alt}\t{event.quality}\t{event.filter}\t"
                        f"{info_field}\t{event.format}\t{sample_field}\n"
                    )
