import datetime
from .sv_merge_logic import should_merge
from .sv_selector import select_representative_sv
from .tra_merger import TRAMerger


class SVMerger:
    def __init__(
            self,
            classified_events,
            tra_delta=50,
            tra_min_overlap_ratio=0.5,
            tra_strand_consistency=True,
            max_distance=50,
            max_length_ratio=1.3,
            min_jaccard=0.7,
    ):
        self.classified_events = classified_events
        self.merged_events: dict[str, dict[str, list]] = {}
        self.event_groups: dict[str, dict[str, list[list]]] = {}
        self.tra_merger = TRAMerger(tra_delta, tra_min_overlap_ratio, tra_strand_consistency)
        self.max_distance = max_distance
        self.max_length_ratio = max_length_ratio
        self.min_jaccard = min_jaccard

    def merge(self):
        for sv_type, chromosomes in self.classified_events.items():
            if sv_type == "TRA":
                for (_chr1, _chr2), events in chromosomes.items():
                    for event in events:
                        self.tra_merger.add_event(event)
            else:
                if sv_type not in self.merged_events:
                    self.merged_events[sv_type] = {}
                    self.event_groups[sv_type] = {}
                for chromosome, events in chromosomes.items():
                    if chromosome not in self.merged_events[sv_type]:
                        self.merged_events[sv_type][chromosome] = []
                        self.event_groups[sv_type][chromosome] = []
                    for event in events:
                        self.add_and_merge_event(sv_type, chromosome, event)

    def add_and_merge_event(self, sv_type, chromosome, new_event):
        events = self.merged_events[sv_type][chromosome]
        event_groups = self.event_groups[sv_type][chromosome]
        for idx, existing_event in enumerate(events):
            if should_merge(existing_event, new_event, self.max_distance, self.max_length_ratio, self.min_jaccard):
                event_groups[idx].append(new_event)
                return
        events.append(new_event)
        event_groups.append([new_event])

    def get_events(self, sv_type, chromosome, start, end):
        if sv_type == "TRA":
            return self.tra_merger.get_merged_events()
        if sv_type in self.event_groups and chromosome in self.event_groups[sv_type]:
            events = []
            for sv_group in self.event_groups[sv_type][chromosome]:
                representative_sv = select_representative_sv(sv_group)
                if representative_sv.start_pos <= end and representative_sv.end_pos >= start:
                    events.append(representative_sv)
            return events
        return []

    def get_all_merged_events(self):
        merged_events = []
        for sv_type in self.event_groups:
            for chromosome in self.event_groups[sv_type]:
                for sv_group in self.event_groups[sv_type][chromosome]:
                    representative_sv = select_representative_sv(sv_group)
                    merged_events.append(representative_sv)
        return merged_events

    def get_events_by_source(self, sources, operation="union"):
        tra_events = self.tra_merger.get_merged_events()
        other_events = self.get_all_merged_events()

        if operation == "union":
            tra_filtered = [event for event in tra_events if any(source in event.source_file for source in sources)]
            other_filtered = [
                event for event in other_events if any(source in event.source_file.split(",") for source in sources)
            ]
        elif operation == "intersection":
            tra_filtered = [event for event in tra_events if all(source in event.source_file for source in sources)]
            other_filtered = [
                event for event in other_events if all(source in event.source_file.split(",") for source in sources)
            ]
        elif operation == "specific":
            tra_filtered = [event for event in tra_events if set(event.source_file.split(",")) == set(sources)]
            other_filtered = [event for event in other_events if set(event.source_file.split(",")) == set(sources)]
        else:
            msg = f"Unsupported operation: {operation}"
            raise ValueError(msg)

        return other_filtered + tra_filtered

    def get_events_by_exact_support(self, exact_support):
        """Get events supported by exactly N files."""
        tra_events = self.tra_merger.get_merged_events()
        other_events = self.get_all_merged_events()

        tra_filtered = [event for event in tra_events if len(set(event.source_file.split(","))) == exact_support]
        other_filtered = [event for event in other_events if len(set(event.source_file.split(","))) == exact_support]
        return other_filtered + tra_filtered

    def get_events_by_support_range(self, min_support=None, max_support=None):
        """Get events supported by a range of files."""
        tra_events = self.tra_merger.get_merged_events()
        other_events = self.get_all_merged_events()

        def within_range(event):
            support_count = len(set(event.source_file.split(",")))
            if min_support is not None and support_count < min_support:
                return False
            if max_support is not None and support_count > max_support:
                return False
            return True

        tra_filtered = [event for event in tra_events if within_range(event)]
        other_filtered = [event for event in other_events if within_range(event)]
        return other_filtered + tra_filtered

    def evaluate_expression(self, expression, event_sources):
        """Evaluate a logical expression against event sources."""
        # Convert file paths to simple identifiers (A, B, C, etc.)
        source_map = {chr(65 + i): str(source) for i, source in enumerate(set(event_sources))}
        reverse_map = {v: k for k, v in source_map.items()}

        # Replace file paths with identifiers in the expression
        expr = expression
        for file_path, identifier in reverse_map.items():
            expr = expr.replace(str(file_path), identifier)

        # Convert event sources to identifiers
        event_identifiers = {reverse_map[str(source)] for source in event_sources}

        # Create evaluation context
        context = {identifier: identifier in event_identifiers for identifier in source_map.keys()}

        # Convert expression to Python boolean expression
        expr = expr.replace("AND", "and").replace("OR", "or").replace("NOT", "not")

        try:
            return eval(expr, {"__builtins__": {}}, context)
        except Exception as e:
            raise ValueError(f"Invalid expression: {e}")

    def get_events_by_expression(self, expression):
        """Get events that satisfy a logical expression."""
        tra_events = self.tra_merger.get_merged_events()
        other_events = self.get_all_merged_events()

        tra_filtered = [
            event for event in tra_events
            if self.evaluate_expression(expression, event.source_file.split(","))
        ]
        other_filtered = [
            event for event in other_events
            if self.evaluate_expression(expression, event.source_file.split(","))
        ]
        return other_filtered + tra_filtered

    def format_sample_values(self, format_keys, sample_dict):
        """Format sample values according to the FORMAT field."""
        values = []
        for key in format_keys:
            value = sample_dict.get(key, ".")
            if isinstance(value, list):
                value = ",".join(map(str, value))
            elif value is None:
                value = "."
            values.append(str(value))

        # Remove the :.: at the end if it exists
        result = ":".join(values)
        if result.endswith(":.:."):
            result = result[:-4]
        return result

    def write_results(self, output_file, events, contigs):
        with open(output_file, "w") as f:
            # Write VCF header
            f.write("##fileformat=SVCFv1.0\n")
            file_date = datetime.datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|")
            f.write(f"##fileDate={file_date}\n")
            f.write("##source=octopusV\n")

            # Write contig information
            for contig_id, contig_length in contigs.items():
                f.write(f"##contig=<ID={contig_id},length={contig_length}>\n")

            # Write column headers
            if events:
                sample_names = ["SAMPLE"]  # We'll use a single column for all samples
            else:
                sample_names = ["SAMPLE"]

            header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n"
            f.write(header_line)

            for event in events:
                # Prepare INFO field
                info_field = ";".join([f"{k}={v}" for k, v in event.info.items()])
                if "SOURCES" not in info_field:
                    info_field += f";SOURCES={event.source_file}"

                # Get the FORMAT field from the representative event
                format_field = event.format
                format_keys = format_field.split(":")

                # Prepare sample values
                merged_samples = getattr(event, "merged_samples", [])
                if merged_samples:
                    # Find the representative sample and other samples
                    rep_sample_data = None
                    other_samples = []
                    for sample_name, sample_format, sample_data in merged_samples:
                        if isinstance(sample_data, dict) and sample_data.get("ID") == event.sv_id:
                            rep_sample_data = (sample_name, sample_format, sample_data)
                        else:
                            other_samples.append((sample_name, sample_format, sample_data))

                    # Format samples in the correct order (representative first, then others)
                    sample_strings = []

                    # Add representative sample first
                    if rep_sample_data:
                        _, _, sample_data = rep_sample_data
                        if isinstance(sample_data, dict):
                            values = []
                            for key in format_keys:
                                value = sample_data.get(key, ".")
                                values.append(str(value))
                            sample_str = ":".join(values)
                            if sample_str.endswith(":.:."):
                                sample_str = sample_str[:-4]
                            sample_strings.append(sample_str)
                        else:
                            if str(sample_data).endswith(":.:."):
                                sample_data = str(sample_data)[:-4]
                            sample_strings.append(str(sample_data))

                    # Add other samples
                    for _, _, sample_data in other_samples:
                        if isinstance(sample_data, dict):
                            values = []
                            for key in format_keys:
                                value = sample_data.get(key, ".")
                                values.append(str(value))
                            sample_str = ":".join(values)
                            if sample_str.endswith(":.:."):
                                sample_str = sample_str[:-4]
                            sample_strings.append(sample_str)
                        else:
                            if str(sample_data).endswith(":.:."):
                                sample_data = str(sample_data)[:-4]
                            sample_strings.append(str(sample_data))

                    sample_part = "\t".join(sample_strings)
                elif hasattr(event, "sample"):
                    formatted_values = self.format_sample_values(format_keys, event.sample)
                    if formatted_values.endswith(":.:."):
                        formatted_values = formatted_values[:-4]
                    sample_part = formatted_values
                else:
                    sample_part = "./."

                # Write the record
                record_part1 = f"{event.chrom}\t{event.pos}\t{event.sv_id}\t{event.ref}\t{event.alt}\t"
                record_part2 = f"{event.quality}\t{event.filter}\t{info_field}\t{format_field}\t"
                f.write(record_part1 + record_part2 + sample_part + "\n")