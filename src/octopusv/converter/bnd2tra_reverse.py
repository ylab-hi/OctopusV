from .base import Converter, get_alt_chrom_pos, get_bnd_pattern


class BND_to_TRA_Reverse_Converter(Converter):
    """This class inherits from the `Converter` base class and implements the conversion logic for BND to TRA reverse conversion."""

    def convert(self, event):
        try:
            if event.is_BND():
                pattern = get_bnd_pattern(event.alt)
                if pattern in ["[p[t", "t]p]"]:
                    chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
                    if chrom_alt is None:
                        print("Failed to get ALT chrom and pos")
                    elif event.chrom == chrom_alt:
                        if pattern == "t]p]" and event.pos > pos_alt:
                            end = pos_alt
                            self.convert_to_TRA_reverse(event, end)
                        elif pattern == "[p[t" and event.pos < pos_alt:
                            end = pos_alt
                            self.convert_to_TRA_reverse(event, end)
        except Exception as e:
            print("Failed to convert BND to Translocation: ", e)

    def convert_to_TRA_reverse(self, event, end):
        event.info["SVTYPE"] = "TRA"
        event.info["END"] = end
        event.info["SVLEN"] = 0
        event.info["CHR2"] = event.chrom
