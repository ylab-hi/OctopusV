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
        self.bnd_pattern = self._extract_bnd_pattern()

        try:
            self.start_chrom, self.start_pos, self.end_chrom, self.end_pos = self._parse_coordinates()
        except Exception as e:
            print(f"Warning: Error parsing coordinates for {sv_id} in {source_file}: {e}")
            # # Use default values as fallback options
            self.start_chrom, self.start_pos = chrom, int(pos)
            self.end_chrom = self.info.get('CHR2', chrom)
            self.end_pos = int(self.info.get('END', pos))

    def _extract_bnd_pattern(self):
        if self.sv_type == 'BND' or self.sv_type == 'TRA':
            return self.alt
        return None

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
        Handles fields that may contain colons (like ALT and CO).
        """
        format_keys = self.format.split(':')
        sample_parts = sample.split(':')
        result = {}

        # Special handling of fields that may contain colons
        special_fields = ['ALT', 'CO', 'REF']  # Fields that might need special handling.
        special_field_index = None

        for i, key in enumerate(format_keys):
            if key in special_fields:
                special_field_index = i
                break

        if special_field_index is not None:
            # Regular fields
            for i in range(special_field_index):
                result[format_keys[i]] = sample_parts[i]

            # Process the special field and all fields after it
            result[format_keys[special_field_index]] = ':'.join(sample_parts[special_field_index:])
        else:
            # If there are no special fields, handle them in the usual way
            result = dict(zip(format_keys, sample_parts))

        return result


    def _parse_coordinates(self):
        """
        Extracts and parses the coordinates of the SV from the sample info or INFO field.
        """
        if self.sv_type == 'BND':
            # For BND types, resolve the second position from the ALT field
            alt_parts = self.alt.split(':')
            if len(alt_parts) > 1:
                end_chrom = alt_parts[0].split(']')[-1].split('[')[-1]
                end_pos = alt_parts[1].rstrip('[]')
            else:
                # If the ALT field format is not as expected, try to get the information from the INFO field
                end_chrom = self.info.get('CHR2', self.chrom)
                end_pos = self.info.get('END', self.pos)
            return self.chrom, int(self.pos), end_chrom, int(end_pos)
        else:
            # For other types try to use CO field, if not then use INFO field
            co = self.sample.get('CO')
            if co and '-' in co:
                start, end = co.split('-')
                start_chrom, start_pos = start.split('_')
                end_chrom, end_pos = end.split('_')
            else:
                start_chrom = self.chrom
                start_pos = self.pos
                end_chrom = self.info.get('CHR2', self.chrom)
                end_pos = self.info.get('END', self.pos)
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