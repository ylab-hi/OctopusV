class SVCFtoBEDConverter:
    def __init__(self, events):
        self.events = events

    def convert(self):
        bed_content = ""
        for event in self.events:
            bed_content += self._convert_event_to_bed(event)
        return bed_content

    def _convert_event_to_bed(self, event):
        chrom = event.chrom
        start = event.start_pos - 1  # BED format is 0-based
        end = event.end_pos
        name = f"{event.id}_{event.sv_type}"
        score = '.'
        strand = event.strand if hasattr(event, 'strand') else '.'

        if event.sv_type == 'INS':
            end = start + 1  # Insertions are represented as a single base pair in BED

        return f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n"