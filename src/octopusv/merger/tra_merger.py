class TRAMerger:
    def __init__(self):
        self.tra_events = {}

    def add_event(self, chr1, pos1, chr2, pos2, source_file):
        key = (chr1, chr2)
        if key not in self.tra_events:
            self.tra_events[key] = []
        self.tra_events[key].append((pos1, pos2, source_file))

    def merge_events(self, distance_threshold=100):
        merged_events = {}
        for key, events in self.tra_events.items():
            merged = []
            for event in events:
                found_match = False
                for merged_event in merged:
                    if (abs(event[0] - merged_event[0]) <= distance_threshold and
                            abs(event[1] - merged_event[1]) <= distance_threshold):
                        merged_event[2].add(event[2])
                        found_match = True
                        break
                if not found_match:
                    merged.append([event[0], event[1], {event[2]}])
            merged_events[key] = merged
        return merged_events

    def get_merged_events(self, distance_threshold=100):
        merged = self.merge_events(distance_threshold)
        return [
            ('TRA', chr1, chr2, pos1, pos2, sources)
            for (chr1, chr2), events in merged.items()
            for pos1, pos2, sources in events
        ]