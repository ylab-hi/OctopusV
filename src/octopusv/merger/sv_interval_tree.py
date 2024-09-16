from intervaltree import Interval, IntervalTree

class SVInterval(Interval):
    """Defines an interval. The interval tree consists of many intervals, and 'data' is the annotation attribute."""
    def __init__(self, begin, end, data):
        super().__init__(begin, end, data)

    def __repr__(self):
        return f"SVInterval({self.begin}, {self.end}, {self.data})"

class SVIntervalTree:
    def __init__(self):
        self.trees = {} # {sv_type: {chromosome: IntervalTree()}}, for the same type and chromosome, there is one tree

    """This add_event method is used to add a structural variation (SV) event to the SVIntervalTree data structure.
        This data structure is a nested dictionary, where the outer dictionary's keys are the structural variation types (sv_type),
        and the values are another dictionary;
        The inner dictionary's keys are chromosome names (chromosome), and the values are IntervalTree objects. 
        Each interval in the interval tree represents a structural variation event,
        and contains data related to that event, such as source file information.
        """

    def add_event(self, sv_type, chromosome, start, end, source_file): # Adding new events and updating overlapping interval data
        if sv_type not in self.trees:
        if sv_type not in self.trees:
            self.trees[sv_type] = {}
        if chromosome not in self.trees[sv_type]:
            self.trees[sv_type][chromosome] = IntervalTree()

        """existing = self.trees[sv_type][chromosome][start:end] Here, [start:end] is not a traditional slicing operation. In IntervalTree,
                this operation is used to query intervals within the given start and end range. If there are overlapping intervals, 
                existing will be an iterable object containing SVInterval objects that meet the criteria,
                each object indeed has three elements: start position, end position, and data. If there are no overlapping intervals, 
                existing will be an empty iterable object, such as an empty list.
                """
        existing = self.trees[sv_type][chromosome][start:end]  # Check if the newly added single SV event (i.e., interval) overlaps with existing ones
        if existing:  # existing contains the overlapping intervals
            # If an overlapping interval exists, update its data
            for interval in existing:
                if isinstance(interval.data, set):  # If data is already a set, directly add the new source
                    interval.data.add(
                        source_file)  # Can directly use the set's add method to quickly add new source file, which is a very efficient operation
                else:  # If not a set (possibly a single value), create a new set containing the original value and the new source
                    interval.data = {interval.data, source_file}

        # Regardless of overlap, add a new interval object to the interval tree, facilitating later merge_overlaps to merge intervals and take the union
        self.trees[sv_type][chromosome].add(
            SVInterval(start, end, {source_file}))  # Created a new SVInterval object and added it to the interval tree

    def merge_overlaps(self):  # Responsible for updating and merging overlapping interval ranges
        for sv_type in self.trees:
            for chromosome in self.trees[sv_type]:  # This is taking each interval tree
                self.trees[sv_type][chromosome].merge_overlaps(
                    data_reducer=lambda x, y: x.union(y) if isinstance(x, set) and isinstance(y, set) else {x, y}
                )  # merge_overlaps is a built-in method of interval tree, it will automatically
                # merge overlapping intervals, data_reducer will define how to handle the "data" part, i.e., the file source part

    def query(self, sv_type, chromosome, start, end):  # All intervals falling between the given start and end points
        if sv_type in self.trees and chromosome in self.trees[sv_type]:
            return self.trees[sv_type][chromosome][start:end]
        return []

    def get_all_events(self):  # Retrieve all intervals from the interval tree and return them as a list
        all_events = []
        for sv_type in self.trees:
            for chromosome in self.trees[sv_type]:
                for interval in self.trees[sv_type][chromosome]:
                    all_events.append((sv_type, chromosome, interval.begin, interval.end, interval.data))
        return all_events  # Each event is represented as a tuple: containing 5 elements,
        # the structure of each event is as follows, (sv_type, chromosome, start, end, sources)

    def get_events_by_source(self, sources, operation='union'):
        all_events = self.get_all_events()
        if operation == 'union':
            return [event for event in all_events if any(source in event[4] for source in sources)]
        elif operation == 'intersection':
            return [event for event in all_events if all(source in event[4] for source in sources)]
        elif operation == 'specific':
            return [event for event in all_events if set(event[4]) == set(sources)]
        else:
            raise ValueError(f"Unsupported operation: {operation}")

    def get_events_by_overlap(self, min_overlap):
        """
        Get events that are supported by at least the specified number of source files.
        Args:
            min_overlap (int): The minimum number of source files that must support an event.
        Returns:
            List: A list of events that meet the overlap criteria.
        """
        all_events = self.get_all_events()
        return [event for event in all_events if len(event[4]) >= min_overlap]


"""
add_event() adds a new interval regardless of whether there's an overlap or not
Using your new `add_event` method and simulated data, let's track the processing step by step and explain the final data structure.
### Input data:
- **sample1.vcf**:
  ```plaintext
  DUP:chr1:100-200
  ```
- **sample2.vcf**:
  ```plaintext
  DUP:chr1:150-250
  ```
- **sample3.vcf**:
  ```plaintext
  DUP:chr1:180-220
  ```
1. **Processing sample1.vcf**:
   - Call `add_event('DUP', 'chr1', 100, 200, 'sample1.vcf')`.
   - Since `'DUP'` and `'chr1'` don't exist in `self.trees` yet, new entries and a new `IntervalTree` will be created.
   - Add new interval `SVInterval(100, 200, {'sample1.vcf'})` to the tree.
2. **Processing sample2.vcf**:
   - Call `add_event('DUP', 'chr1', 150, 250, 'sample2.vcf')`.
   - At this point, there's an existing interval `SVInterval(100, 200, {'sample1.vcf'})` that overlaps with the new interval `150-250`.
   - For the overlapping interval (i.e., `SVInterval(100, 200, {'sample1.vcf'})`), update the data set by adding `'sample2.vcf'`.
   - Add new interval `SVInterval(150, 250, {'sample2.vcf'})` to the tree.
3. **Processing sample3.vcf**:
   - Call `add_event('DUP', 'chr1', 180, 220, 'sample3.vcf')`.
   - At this point, two intervals overlap: `SVInterval(100, 200, {'sample1.vcf', 'sample2.vcf'})` and `SVInterval(150, 250, {'sample2.vcf'})`.
   - Update the data sets of these two overlapping intervals by adding `'sample3.vcf'`.
   - Add new interval `SVInterval(180, 220, {'sample3.vcf'})` to the tree.
### Final data structure (before calling `merge_overlaps`):
In `self.trees['DUP']['chr1']`, you will have the following intervals:
- `SVInterval(100, 200, {'sample1.vcf', 'sample2.vcf', 'sample3.vcf'})`
- `SVInterval(150, 250, {'sample2.vcf', 'sample3.vcf'})`
- `SVInterval(180, 220, {'sample3.vcf'})`
### Result after calling `merge_overlaps`:
- `merge_overlaps` will simplify these intervals by merging overlapping ones.
- It will merge all overlapping intervals and combine their data sets into one.
- The result will be a single interval:
  - `SVInterval(100, 250, {'sample1.vcf', 'sample2.vcf', 'sample3.vcf'})`
"""
"""
# Assume we have the following intervals:
# SVInterval(100, 200, {'sample1.vcf', 'sample2.vcf'})
# SVInterval(150, 250, {'sample2.vcf'})
# merge_overlaps will perform the following operations:
# 1. Detect that these two intervals overlap
# 2. Merge the interval range: (100, 250)
# 3. Use data_reducer to merge the data:
#    {'sample1.vcf', 'sample2.vcf'}.union({'sample2.vcf'})
#    Result: {'sample1.vcf', 'sample2.vcf'}
# Final result:
# SVInterval(100, 250, {'sample1.vcf', 'sample2.vcf'})
"""

"""
After merge_overlaps:
self.trees = {
    'DUP': {
        'chr1': IntervalTree([
            SVInterval(100, 250, {'sample1.vcf', 'sample2.vcf', 'sample3.vcf'})
        ]),
        'chr4': IntervalTree([
            SVInterval(700, 850, {'sample2.vcf', 'sample3.vcf'})
        ])
    },
    'DEL': {
        'chr2': IntervalTree([
            SVInterval(300, 400, {'sample1.vcf', 'sample2.vcf'})
        ])
    },
    'INV': {
        'chr3': IntervalTree([
            SVInterval(500, 600, {'sample1.vcf', 'sample3.vcf'})
        ])
    }
}
"""

