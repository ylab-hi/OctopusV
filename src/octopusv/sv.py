from enum import Enum

from octopusv.utils.construct_sample_string import construct_sample_string


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
        self.orig_alt = alt
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
                if key.lower() == "strand":
                    key = "STRAND"
                info_dict[key] = value
            else:
                info_dict[parts[0]] = None

        for support_key in ["SUPP", "SUPPREAD", "WEIGHT", "RE"]:
            if support_key in info_dict:
                info_dict["SUPPORT"] = info_dict[support_key]
                break

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
        # Process 'PR' for 'SUPPORT' if not already set
        if "SUPPORT" not in self.info:
            format_parts = self.format.split(":")
            sample_parts = self.sample.split(":")
            format_sample_dict = dict(zip(format_parts, sample_parts))
            # Handling both 'PR' only and 'PR:SR' cases.
            pr_key = "PR" if "PR" in format_sample_dict else None
            if pr_key:
                pr_values = format_sample_dict[pr_key].split(",")
                if len(pr_values) > 1:
                    self.info["SUPPORT"] = pr_values[1]  # Use the alt allele's paired-read count

        # Fixed order for INFO fields, using '.' as a placeholder for missing values
        info_order = ["SVTYPE", "END", "SVLEN", "CHR2", "SUPPORT", "SVMETHOD", "RTID", "AF", "STRAND"]
        info_str_parts = []
        for key in info_order:
            value = self.info.get(key, ".")
            if key == "SVLEN" and value != ".":  # Check if SVLEN is present and not placeholder
                try:
                    # Ensure the value is an integer and take its absolute value
                    value = str(abs(int(value)))
                except ValueError:
                    # In case value is not an integer, keep the original value (should not happen if your data is correct)
                    pass
            info_str_parts.append(f"{key}={value}")
            info_str = ";".join(info_str_parts)

        # Fixed order for FORMAT fields
        format_order = ["GT", "AD", "LN", "ST", "QV", "TY", "ID", "SC", "REF", "ALT", "CO"]
        format_str = ":".join(format_order)

        # Use construct_sample_string to build sample string
        sample_str = construct_sample_string(self)

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
