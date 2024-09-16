from .sv_interval_tree import SVIntervalTree

class SVMerger:
    def __init__(self, classified_events):
        self.classified_events = classified_events
        self.interval_tree = SVIntervalTree()

    def merge(self):
        for sv_type, chromosomes in self.classified_events.items():
            for chromosome, events in chromosomes.items():
                for event in events:
                    self.interval_tree.add_event(
                        sv_type, chromosome, event.start_pos, event.end_pos, event.source_file
                    )
        self.interval_tree.merge_overlaps()

    def get_events(self, sv_type, chromosome, start, end):
        return self.interval_tree.query(sv_type, chromosome, start, end)

    def get_events_by_source(self, sources, operation='union'):
        return self.interval_tree.get_events_by_source(sources, operation)

    def get_events_by_overlap(self, min_overlap):
        return self.interval_tree.get_events_by_overlap(min_overlap)

    def write_results(self, output_file, events):
        with open(output_file, 'w') as f:
            for sv_type, chromosome, start, end, sources in events:
                f.write(f"{sv_type}\t{chromosome}\t{start}\t{end}\t{','.join(sources)}\n")