class SVCFtoBEDPEConverter:
    def __init__(self, events):
        self.events = events

    def convert(self):
        bedpe_content = "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\n"
        for event in self.events:
            bedpe_content += self._convert_event_to_bedpe(event)
        return bedpe_content

    def _convert_event_to_bedpe(self, event):
        chrom1 = event.chrom
        start1 = event.start_pos - 1  # BEDPE format is 0-based
        end1 = event.start_pos
        chrom2 = event.chr2 if event.sv_type == 'TRA' else chrom1
        start2 = event.end_pos - 1
        end2 = event.end_pos
        name = f"{event.id}_{event.sv_type}"
        score = '.'
        strand1 = '+'
        strand2 = '-' if event.sv_type in ['INV', 'TRA'] else '+'

        if event.sv_type == 'INS':
            start2 = start1
            end2 = end1

        return f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{name}\t{score}\t{strand1}\t{strand2}\n"