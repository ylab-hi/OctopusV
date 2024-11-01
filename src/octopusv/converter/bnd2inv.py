import logging

from .base import Converter, get_alt_chrom_pos, get_bnd_pattern


class BND_to_INV_Converter(Converter):
    """This class inherits from the `Converter` base class and implements the conversion logic for BND to INV conversion."""

    def convert(self, event):  # Each converter has a convert function
        try:
            if not event.is_BND():
                return

            pattern = get_bnd_pattern(event.alt)
            if pattern in ("t]p]", "[p[t"):
                chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
                if chrom_alt is None:
                    logging.info("Failed to get ALT chrom and pos")
                elif event.chrom == chrom_alt:  # Do this only when same chromosome
                    if pattern == "t]p]" and event.pos < pos_alt:
                        end = pos_alt
                        svlen = abs(event.pos - pos_alt)
                        self.convert_to_inv(
                            event, end, svlen, chrom_alt
                        )  # self will be passed to convert_to_INV method.
                    elif pattern == "[p[t" and event.pos > pos_alt:
                        end = event.pos
                        event.pos = pos_alt
                        svlen = abs(event.pos - end)
                        self.convert_to_inv(event, end, svlen, chrom_alt)

        except Exception as e:
            logging.error("Failed to convert BND to Inversion: ", e)

    def convert_to_inv(self, event, end, svlen, chrom_alt):
        # Detail logic changing BND event to INV event
        event.alt = "<INV>"
        event.info["SVTYPE"] = "INV"
        event.info["END"] = end
        event.info["SVLEN"] = svlen
        event.info["CHR2"] = chrom_alt
        event.info["SVMETHOD"] = "octopusV"
