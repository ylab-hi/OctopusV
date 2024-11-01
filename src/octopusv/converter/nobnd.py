from .base import Converter


class NonBNDConverter(Converter):
    """A converter for non-BND events that enriches the event with additional info."""

    def convert(self, event):
        """Simplifies SVTYPE to one of the specified types and sets CHR2 and SVMETHOD for non-BND events."""
        # Check and simplify SVTYPE
        svtype = event.info.get("SVTYPE", "")
        if ":" in svtype:
            # Only take the part before the colon as SVTYPE
            svtype = svtype.split(":")[0]
        # Ensure the simplified SVTYPE is one of the allowed five types, otherwise keep the original SVTYPE unchanged
        if svtype in ["INV", "INS", "DEL", "TRA", "DUP"]:
            event.info["SVTYPE"] = svtype
        else:
            # If SVTYPE is not among the specified types, print a warning or keep it unchanged
            logging.warning(f"Simplified SVTYPE '{svtype}' is not one of the allowed types. Keeping original SVTYPE.")

        # Only set CHR2 if it does not already exist
        if "CHR2" not in event.info:
            event.info["CHR2"] = event.chrom

        # Set SVMETHOD to "octopusV"
        event.info["SVMETHOD"] = "octopusV"

        # Calculate SVLEN if not present and event is not TRA
        if "SVLEN" not in event.info and svtype in ["DEL", "DUP", "INV"]:
            try:
                pos = int(event.pos)
                end = int(event.info["END"])
                event.info["SVLEN"] = abs(end - pos)
            except (ValueError, KeyError):
                logging.error(f"Unable to calculate SVLEN for event {event.id}")

        if svtype == "TRA":
            event.info["SVLEN"] = "."

        # Deal with the "END" problem in INS
        if svtype == "INS":
            pos = int(event.pos)
            svlen = int(event.info.get("SVLEN", "0"))
            event.info["END"] = str(pos + svlen)
