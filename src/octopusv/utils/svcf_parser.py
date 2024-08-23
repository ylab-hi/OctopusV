class SVCFEvent:
    """
    Represents a structural variant (SV) event parsed from an SVCF file.

    Attributes:
        chrom (str): Chromosome on which the SV event occurs.
        pos (int): Starting position of the SV on the chromosome.
        sv_id (str): Unique identifier for the SV event.
        ref (str): Reference bases at the SV position.
        alt (str): Alternate bases indicating the SV.
        quality (str): Quality score of the SV event.
        filter (str): Filter status of the SV event.
        info (dict): Dictionary containing additional information about the SV event.
        format (str): Format of the sample data related to the SV event.
        sample (str): Sample-specific data for the SV event.
        source_file (str): The source file from which this SV event was parsed.
    """

    def __init__(self, chrom, pos, sv_id, ref, alt, quality, filter, info, format, sample, source_file):
        self.chrom = chrom
        self.pos = int(pos)
        self.sv_id = sv_id
        self.ref = ref
        self.alt = alt
        self.quality = quality
        self.filter = filter
        self.info = self._parse_info(info)
        self.format = format
        self.sample = self._parse_sample(sample)
        self.source_file = source_file  # Store the source file name for tracking

        self.sv_type = self.info.get('SVTYPE', '')
        self.start_chrom, self.start_pos, self.end_chrom, self.end_pos = self._parse_coordinates()

    def _parse_info(self, info_str):
        """
        Parses the info field from a VCF record into a dictionary.

        Args:
            info_str (str): The raw info string from a VCF record.

        Returns:
            dict: A dictionary where each key is an info field name and the value is the field value.
        """
        info = {}
        for item in info_str.split(';'):
            if '=' in item:
                key, value = item.split('=')
                info[key] = value
            else:
                info[item] = True  # Flags without a value are stored as True
        return info

    def _parse_sample(self, sample):
        """
        Parses the sample string based on the format to extract all relevant data.
        """
        parts = sample.split(':')
        keys = self.format.split(':')
        return dict(zip(keys, parts))

    def _parse_coordinates(self):
        """
        Extracts and parses the coordinates of the SV from the sample info.
        """
        co = self.sample.get('CO', f"{self.chrom}_{self.pos}-{self.chrom}_{self.pos}")
        start, end = co.split('-')
        start_chrom, start_pos = start.split('_')
        end_chrom, end_pos = end.split('_')
        return start_chrom, int(start_pos), end_chrom, int(end_pos)


class SVCFFileEventCreator:
    """
    Parses all SVCF files to create SVCFEvent instances.

    Attributes:
        filenames (list): List of filenames (strings) to be parsed.
        events (list): List of SVCFEvent objects parsed from the files.
    """

    def __init__(self, filenames):
        self.filenames = filenames
        self.events = []

    def parse(self):
        """
        Parses each file in the filenames list, creating SVCFEvent objects for each valid SV record.
        """
        for filename in self.filenames:
            with open(filename, 'r') as file:
                for line in file:
                    if line.startswith('#'):
                        continue  # Skip comment lines that start with '#'
                    parts = line.strip().split()
                    if len(parts) < 10:
                        continue  # Ensure that there are enough parts to form a complete SV event
                    sv_event = SVCFEvent(*parts, source_file=filename)  # Create an SV event object
                    self.events.append(sv_event)  # Add the event to the list of events
