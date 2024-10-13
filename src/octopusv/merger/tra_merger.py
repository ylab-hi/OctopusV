class TRAMerger:
    def __init__(self, distance_threshold=100):
        # Dictionary to store TRA events, keyed by chromosome pairs
        self.tra_events = {}
        # Maximum allowed distance between events to be considered for merging
        self.distance_threshold = distance_threshold

    def add_event(self, chr1, pos1, chr2, pos2, source_file, bnd_pattern):
        # Create a key from the chromosome pair
        key = (chr1, chr2)
        # Initialize an empty list for this chromosome pair if it doesn't exist
        if key not in self.tra_events:
            self.tra_events[key] = []
        # Add the new event to the list for this chromosome pair
        self.tra_events[key].append((pos1, pos2, source_file, bnd_pattern))

    def merge_events(self):
        # Dictionary to store all merged events, grouped by chromosome pairs
        all_chromosome_pair_events = {}

        # Iterate through each chromosome pair and its events
        for chromosome_pair, unmerged_events in self.tra_events.items():
            # List to store merged events for this specific chromosome pair
            merged_events_for_certain_chrom_pair = []

            # Process each unmerged event
            for current_event in unmerged_events:
                event_was_merged = False

                # Check if this event can be merged with any existing merged event
                for existing_merged_event in merged_events_for_certain_chrom_pair:
                    # If both positions are within the threshold distance
                    if (abs(current_event[0] - existing_merged_event[0]) <= self.distance_threshold and
                            abs(current_event[1] - existing_merged_event[1]) <= self.distance_threshold and
                            self._is_compatible_bnd_pattern(current_event[3], existing_merged_event[3])):
                        # Add the source file to the merged event
                        # current_event[2] is the new source_file, existing_merged_event[2] is the set of existing source_files
                        existing_merged_event[2].add(current_event[2])
                        event_was_merged = True
                        break

                # If no match found, create a new merged event
                if not event_was_merged:
                    merged_events_for_certain_chrom_pair.append(
                        [current_event[0], current_event[1], {current_event[2]}, current_event[3]])

            # Store the merged events for this chromosome pair
            all_chromosome_pair_events[chromosome_pair] = merged_events_for_certain_chrom_pair

        return all_chromosome_pair_events

    def _is_compatible_bnd_pattern(self, pattern1, pattern2):
        if pattern1 is None or pattern2 is None:
            return False

        def classify_pattern(pattern):
            if pattern.startswith(']') and pattern.endswith('N'):
                return 1  # ]chr:pos]N
            elif pattern.startswith('N[') and pattern.endswith('['):
                return 2  # N[chr:pos[
            elif pattern.startswith('N]') and pattern.endswith(']'):
                return 3  # N]chr:pos]
            elif pattern.startswith('[') and pattern.endswith('N'):
                return 4  # [chr:pos[N
            else:
                return 0  # Unknown pattern

        class1 = classify_pattern(pattern1)
        class2 = classify_pattern(pattern2)

        # If two BND patter are compatible
        return class1 == class2 and class1 != 0

    def get_merged_events(self):
        # Get all merged events for all chromosome pairs
        all_chromosome_pair_events = self.merge_events()

        # Convert merged events into a list of tuples
        # Each tuple contains: ('TRA', chr1, chr2, pos1, pos2, sources, bnd_pattern)
        return [
            ('TRA', chr1, chr2, pos1, pos2, sources, bnd_pattern)
            for (chr1, chr2), merged_events_for_certain_chrom_pair in all_chromosome_pair_events.items()
            for pos1, pos2, sources, bnd_pattern in merged_events_for_certain_chrom_pair
        ]


"""
self.tra_events = {
    ('chr1', 'chr2'): [
        (1000, 2000, 'FileA'),
        (1050, 2030, 'FileB'),
        (1500, 2500, 'FileC'),
        (1520, 2510, 'FileD')
    ],
    ('chr3', 'chr4'): [
        (3000, 4000, 'FileA'),
        (3100, 4100, 'FileB')
    ]
}

self.distance_threshold = 100

all_chromosome_pair_events = {
    ('chr1', 'chr2'): [
        [1000, 2000, {'FileA', 'FileB'}],
        [1500, 2500, {'FileC', 'FileD'}]
    ],
    ('chr3', 'chr4'): [
        [3000, 4000, {'FileA', 'FileB'}]
    ]
}

merged_events_for_certain_chrom_pair = [
    [1000, 2000, {'FileA', 'FileB'}],
    [1500, 2500, {'FileC', 'FileD'}]
]

"""