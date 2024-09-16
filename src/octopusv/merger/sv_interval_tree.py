from intervaltree import Interval, IntervalTree

class SVInterval(Interval): # 这里定义了一个区间，区间树是由很多区间组成的，data是注释属性
    def __init__(self, begin, end, data):
        super().__init__(begin, end, data)

    def __repr__(self):
        return f"SVInterval({self.begin}, {self.end}, {self.data})"

class SVIntervalTree:
    def __init__(self):
        self.trees = {}  # {sv_type: {chromosome: IntervalTree()}},同一个类型，一个染色体，有一棵树

    """这个方法add_event是用来向SVIntervalTree数据结构中添加一个结构变异（Structural Variation，SV）事件的。
    这个数据结构是一个嵌套的字典，其中外层字典的键是结构变异类型（sv_type），值是另一个字典；
    内层字典的键是染色体名称（chromosome），值是一个区间树（IntervalTree）对象。每个区间树中的区间代表一个结构变异事件，
    并且包含了与该事件相关的数据，比如源文件信息。
    """

    def add_event(self, sv_type, chromosome, start, end, source_file): # 函数负责添加新事件和更新重叠区间的数据
        if sv_type not in self.trees:
            self.trees[sv_type] = {}
        if chromosome not in self.trees[sv_type]:
            self.trees[sv_type][chromosome] = IntervalTree()

        """existing = self.trees[sv_type][chromosome][start:end]这里的[start:end]不是传统意义上的切片操作。在IntervalTree中，
        这个操作是用于查询在给定的start和end范围内的区间。如果有重叠的区间，existing将是一个可迭代对象，其中包含了符合条件的SVInterval对象，
        每个对象确实有三要素：起始位置、结束位置和数据。如果没有重叠的区间，existing将是一个空的可迭代对象，比如空列表。
"""

        existing = self.trees[sv_type][chromosome][start:end] # 看新加入的单个SV事件（也就是区间）和现有的重合与否
        if existing: # existing是现有的重合的区间
            # If an overlapping interval exists, update its data
            for interval in existing:
                if isinstance(interval.data, set): # 如果data已经是集合，直接添加新来源。
                    interval.data.add(source_file) # 可以直接使用集合的add方法快速添加新的源文件，这是一个非常高效的操作
                else: # 如果不是集合（可能是单一值），创建一个新集合，包含原有值和新来源。
                    interval.data = {interval.data, source_file}

        # 不管重不重合，都向区间树里面添加一个新区间对象，方便后面merge_overlaps合并区间，取并集
        self.trees[sv_type][chromosome].add(SVInterval(start, end, {source_file})) # 创建了一个新的SVInterval对象并将其添加到区间树中。

    def merge_overlaps(self): # 负责更新并合并重叠的区间范围
        for sv_type in self.trees:
            for chromosome in self.trees[sv_type]: # 这里是取每一个区间树
                self.trees[sv_type][chromosome].merge_overlaps(
                    data_reducer=lambda x, y: x.union(y) if isinstance(x, set) and isinstance(y, set) else {x, y}
                ) # merge_overlaps是interval tree自带方法，会自动把重叠区间合并，data_reducer会定义如何处理"数据"部分，也就是文件来源部分

    def query(self, sv_type, chromosome, start, end): # 落在给定起始点（start）和终止点（end）之间的所有区间
        if sv_type in self.trees and chromosome in self.trees[sv_type]:
            return self.trees[sv_type][chromosome][start:end]
        return []

    def get_all_events(self): # 从区间树中检索出所有的区间，并将它们作为一个列表返回
        all_events = []
        for sv_type in self.trees:
            for chromosome in self.trees[sv_type]:
                for interval in self.trees[sv_type][chromosome]:
                    all_events.append((sv_type, chromosome, interval.begin, interval.end, interval.data))
        return all_events # 每个事件被表示为一个元组：包含5个元素,每个事件的结构如下,(sv_type, chromosome, start, end, sources)

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
        获取至少被指定数量的源文件支持的事件。

        Args:
            min_overlap (int): 支持事件的最小源文件数量。

        Returns:
            List: 符合重叠条件的事件列表。
        """
        all_events = self.get_all_events()
        return [event for event in all_events if len(event[4]) >= min_overlap]


"""
add_event()无论有没有重合都添加一个新区间

使用你的新的 `add_event` 方法和模拟数据，让我们逐步跟踪处理过程并解释最终的数据结构。

### 输入数据：

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

1. **处理 sample1.vcf**:
   - 调用 `add_event('DUP', 'chr1', 100, 200, 'sample1.vcf')`。
   - 由于 `'DUP'` 和 `'chr1'` 在 `self.trees` 中尚不存在，会创建新的条目和新的 `IntervalTree`。
   - 添加新区间 `SVInterval(100, 200, {'sample1.vcf'})` 到树中。

2. **处理 sample2.vcf**:
   - 调用 `add_event('DUP', 'chr1', 150, 250, 'sample2.vcf')`。
   - 此时，已有一个区间 `SVInterval(100, 200, {'sample1.vcf'})` 存在，且与新的区间 `150-250` 重叠。
   - 对于重叠的区间（即 `SVInterval(100, 200, {'sample1.vcf'})`），更新数据集合，添加 `'sample2.vcf'`。
   - 添加新区间 `SVInterval(150, 250, {'sample2.vcf'})` 到树中。

3. **处理 sample3.vcf**:
   - 调用 `add_event('DUP', 'chr1', 180, 220, 'sample3.vcf')`。
   - 此时，两个区间存在重叠：`SVInterval(100, 200, {'sample1.vcf', 'sample2.vcf'})` 和 `SVInterval(150, 250, {'sample2.vcf'})`。
   - 更新这两个重叠区间的数据集合，添加 `'sample3.vcf'`。
   - 添加新区间 `SVInterval(180, 220, {'sample3.vcf'})` 到树中。

### 最终数据结构（在调用 `merge_overlaps` 之前）：

在 `self.trees['DUP']['chr1']` 中，你会有以下的区间：

- `SVInterval(100, 200, {'sample1.vcf', 'sample2.vcf', 'sample3.vcf'})`
- `SVInterval(150, 250, {'sample2.vcf', 'sample3.vcf'})`
- `SVInterval(180, 220, {'sample3.vcf'})`

### 调用 `merge_overlaps` 后的结果：

- `merge_overlaps` 将通过合并重叠区间来简化这些区间。
- 它将合并所有重叠的区间，并将其数据集合合并为一个。
- 结果将是单个区间：
  - `SVInterval(100, 250, {'sample1.vcf', 'sample2.vcf', 'sample3.vcf'})`

"""

"""
# 假设我们有以下区间：
# SVInterval(100, 200, {'sample1.vcf', 'sample2.vcf'})
# SVInterval(150, 250, {'sample2.vcf'})

# merge_overlaps 会执行以下操作：
# 1. 检测到这两个区间重叠
# 2. 合并区间范围：(100, 250)
# 3. 使用 data_reducer 合并数据：
#    {'sample1.vcf', 'sample2.vcf'}.union({'sample2.vcf'})
#    结果：{'sample1.vcf', 'sample2.vcf'}

# 最终结果：
# SVInterval(100, 250, {'sample1.vcf', 'sample2.vcf'})
"""

"""
在 merge_overlaps 之后:
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

