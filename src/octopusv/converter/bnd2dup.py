from .base import Converter, get_alt_chrom_pos, get_bnd_pattern


class BND_to_DUP_Converter(Converter):
    """This class inherits from the `Converter` base class and implements the conversion logic for BND to DUP conversion."""

    def convert(self, event):
        try:
            if event.is_BND():
                pattern = get_bnd_pattern(event.alt)
                if pattern in ["]p]t", "t[p["]:
                    chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
                    if chrom_alt is None:
                        print("Failed to get ALT chrom and pos")
                    elif event.chrom == chrom_alt:
                        if pattern == "]p]t" and event.pos < pos_alt:
                            end = pos_alt
                            svlen = abs(end - event.pos)
                            self.convert_to_DUP(event, end, svlen)
                        elif pattern == "t[p[" and event.pos > pos_alt:
                            end = event.pos
                            svlen = abs(event.pos - pos_alt)
                            event.pos = pos_alt
                            self.convert_to_DUP(event, end, svlen)
        except Exception as e:
            print("Failed to convert BND to Duplication: ", e)

    def convert_to_DUP(self, event, end, svlen):
        event.alt = "<DUP>"
        event.info["SVTYPE"] = "DUP"
        event.info["END"] = end
        event.info["SVLEN"] = svlen
