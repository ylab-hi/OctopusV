from .TRA_merge_logic import should_merge_tra

class TRAMerger:
    def __init__(self, delta=50, min_overlap_ratio=0.5, strand_consistency=True):
        self.tra_events = {}
        self.delta = delta
        self.min_overlap_ratio = min_overlap_ratio
        self.strand_consistency = strand_consistency

    def add_event(self, event):
        key = tuple(sorted([event.start_chrom, event.end_chrom]))
        if key not in self.tra_events:
            self.tra_events[key] = []
        self.tra_events[key].append(event)

    def merge_events(self):
        all_chromosome_pair_events = {}
        for chromosome_pair, unmerged_events in self.tra_events.items():
            merged_events_for_certain_chrom_pair = []
            for current_event in unmerged_events:
                event_was_merged = False
                for existing_merged_event in merged_events_for_certain_chrom_pair:
                    if should_merge_tra(existing_merged_event, current_event, self.delta, self.min_overlap_ratio, self.strand_consistency):
                        # 合并事件
                        existing_merged_event.info['SUPPORT'] = str(int(existing_merged_event.info.get('SUPPORT', 1)) + int(current_event.info.get('SUPPORT', 1)))
                        # 正确处理文件名
                        existing_merged_event.source_file = ','.join(sorted(set(existing_merged_event.source_file.split(',') + [current_event.source_file])))
                        event_was_merged = True
                        break
                if not event_was_merged:
                    merged_events_for_certain_chrom_pair.append(current_event)
            all_chromosome_pair_events[chromosome_pair] = merged_events_for_certain_chrom_pair
        return all_chromosome_pair_events

    def get_merged_events(self):
        all_chromosome_pair_events = self.merge_events()
        merged_events = []
        for chromosome_pair, events in all_chromosome_pair_events.items():
            for event in events:
                merged_events.append(('TRA', event.start_chrom, event.end_chrom, event.start_pos, event.end_pos, event.source_file, event.bnd_pattern))
        return merged_events


"""
self.tra_events = {
    ('chr1', 'chr2'): [
        (1000, 2000, 'FileA'),
        (1050, 2030, 'FileB'),
        (1500, 2500, 'FileC'),
        (1520, 2510, 'FileD')
    ],
    ('chr3', 'chr4'): [
        (3000, 4000, 'FileA'),
        (3100, 4100, 'FileB')
    ]
}

self.distance_threshold = 100

all_chromosome_pair_events = {
    ('chr1', 'chr2'): [
        [1000, 2000, {'FileA', 'FileB'}],
        [1500, 2500, {'FileC', 'FileD'}]
    ],
    ('chr3', 'chr4'): [
        [3000, 4000, {'FileA', 'FileB'}]
    ]
}

merged_events_for_certain_chrom_pair = [
    [1000, 2000, {'FileA', 'FileB'}],
    [1500, 2500, {'FileC', 'FileD'}]
]

"""