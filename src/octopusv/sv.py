from enum import Enum


class SVType(Enum):
    DEL = "DEL"
    DUP = "DUP"
    INV = "INV"
    INS = "INS"
    TRA = "TRA"
    BND = "BND"


class SVEvent:  # A class to represent each SV event
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
        # enable you to output object as vcf string format
        # Define your desired order
        ordered_keys = ["SVTYPE", "END", "SVLEN", "SUPPORT", "COVERAGE", "STRAND"]
        info_items = []

        # First add the ordered keys
        for key in ordered_keys:
            if key in self.info:
                info_items.append(f"{key}={self.info[key]}")

        # Then add the remaining keys
        for key in self.info:
            if key not in ordered_keys:
                info_items.append(f"{key}={self.info[key]}")

        info_str = ";".join(info_items)  # SVTYPE=INV;END=69650961;SVLEN=1835141

        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.chrom,
            self.pos,
            self.id,
            self.ref,
            self.alt,
            self.qual,
            self.filter,
            info_str,
            self.format,
            self.sample,
        )
