from enum import Enum


class SVType(Enum):
    """A class to represent the SV type."""

    DEL = "DEL"
    DUP = "DUP"
    INV = "INV"
    INS = "INS"
    TRA = "TRA"
    BND = "BND"


class SVEvent:
    """represent each SV event."""

    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, format, sample):
        self.chrom = chrom
        self.pos = int(pos)
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = self._parse_info(info)  # self.info will be a dict
        self.format = format
        self.sample = sample

    def _parse_info(self, info):  # Parse the info field according to your requirement
        info_dict = {}
        for item in info.split(";"):
            parts = item.split("=")
            if len(parts) == 2:  # Just in case that the INFO doesn't contain = symbol
                key, value = parts
                info_dict[key] = value
            else:
                info_dict[parts[0]] = None
        return info_dict

    def __getitem__(self, key):
        # enable you to access object as dict
        return self.info[key]

    def __setitem__(self, key, value):
        self.info[key] = value

    def is_duplication(self):
        # Whether is duplication
        return self.info["SVTYPE"] == "DUP"

    def is_inversion(self):
        # Whether is inversion
        return self.info["SVTYPE"] == "INV"

    def is_insertion(self):
        #  Whether is insertion
        return self.info["SVTYPE"] == "INS"

    def is_TRA(self):
        #  Whether is TRA
        return self.info["SVTYPE"] == "TRA"

    def is_BND(self):
        #  Whether is BND
        return self.info["SVTYPE"] == "BND"

    def __str__(self):
        # Fixed order for INFO fields, using '.' as a placeholder for missing values
        info_order = ["SVTYPE", "END", "SVLEN", "CHR2", "SUPPORT", "SVMETHOD", "RTID", "AF", "STRAND"]
        info_str = ";".join(f"{key}={self.info.get(key, '.')}" for key in info_order)

        # Fixed order for FORMAT fields
        format_order = ["GT", "AD", "LN", "ST", "QV", "TY", "ID", "SC", "REF", "ALT", "CO"]
        format_str = ":".join(format_order)

        # Correcting the logic for processing self.sample to match the FORMAT field order
        # Split the input format and sample strings into lists
        input_format_parts = self.format.split(":")
        input_sample_parts = self.sample.split(":")

        # Map the input sample to a dictionary for easy key access
        sample_dict = dict(zip(input_format_parts, input_sample_parts))

        # Generate the corresponding Sample string, filling in missing parts with '.'
        # output_sample_parts is a list of new sample column
        output_sample_parts = [sample_dict.get(part, '.') for part in format_order]
        sample_str = ":".join(output_sample_parts)

        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.chrom,
            self.pos,
            self.id,
            self.ref,
            self.alt,
            self.qual,
            self.filter,
            info_str,
            format_str,
            sample_str,
        )
