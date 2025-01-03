import logging

from .base import Converter, get_alt_chrom_pos, get_bnd_pattern


class BND_to_TRA_Forward_Converter(Converter):
    """This class inherits from the `Converter` base class and implements.

    the conversion logic for BND to TRA Forward translocation cut-paste.
    """

    def convert(self, event):
        try:
            if event.is_BND():
                pattern = get_bnd_pattern(event.alt)
                if pattern in ["t[p[", "]p]t"]:
                    chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
                    if chrom_alt is None:
                        logging.info("Failed to get ALT chrom and pos")
                    elif event.chrom == chrom_alt:
                        if (pattern == "t[p[" and event.pos < pos_alt) or (pattern == "]p]t" and event.pos > pos_alt):
                            end = pos_alt
                            self.convert_to_TRA_forward(event, end)
        except Exception as e:
            logging.error("Failed to convert BND to Translocation: ", e)

    def convert_to_TRA_forward(self, event, end):
        event.info["SVTYPE"] = "TRA"
        event.info["END"] = end
        event.info["SVLEN"] = "."
        event.info["CHR2"] = event.chrom
        event.info["SVMETHOD"] = "octopusV"
