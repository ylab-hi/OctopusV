# octopusV: Correcting and Standardize benchmarking of Structural Variant (SV) Results from SV Callers

## What is octopusV?

octopusV is a toolkit for efficient manipulation and comparative analysis of structural variant (SV) callsets from multiple approaches. It streamlines SV analysis workflows by providing a standardized input/output interface and a comprehensive set of commands for filtration, merging, intersection, benchmarking and visualization of VCF from various SV callers. Octopus aims to offer a holistic solution for researchers venturing into the complex landscape of structural variations.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

- Python 3

### Installation

```bash
pip install octopusv
```

## Usage

### Standardize input with `convert`

Convert VCFs from multiple callers to standardized SVCF format:

```
octopusv convert -i caller1.vcf -o caller1.svcf
octopusv convert -i caller2.vcf -o caller2.svcf
```

### Merge callsets with `merge`

Extract intersection or union of multiple callsets:

```
octopusv merge -i caller1.svcf caller2.svcf -o intersect.svcf --intersect
octopusv merge -i caller1.svcf caller2.svcf -o union.svcf --union
```

Extract SVs specifically supported by vcf1:

```
octopusv merge -i vcf1.svcf vcf2.svcf -o uniq_vcf1.svcf --specific vcf1.svcf
```

Extract SVs specifically supported by vcf1 and vcf2:

```
octopusv merge -i vcf1.svcf vcf2.svcf vcf3.svcf -o common_vcf1_vcf2.svcf --specific vcf1.svcf vcf2.svcf
```

Extract SVs supported by 2+ callers:

```
octopusv merge -i vcf1.svcf vcf2.svcf vcf3.svcf -o overlap.svcf --overlap 2
```

### Evaluate calls against truthsets using `bench`

```
octopusv bench -i candidate.svcf -g truthset.svcf -o bench_stat.txt
```

### Generate statistical reports for the called SVs using `stat`

```
octopusv stat -i merged.svcf -o sv_statistics.txt
```

### Plotting the stat results using `plot`

```
octopusv plot -i sv_statistics.svcf -o plot
```

```
octopusv plot -i bench_statistics.txt -o plot
```

## Contribution

We welcome contributions from the community. Feel free to submit issues, feature requests, and pull requests.

- Install the library

```bash
poetry install
```

- Made and commit your changes

```bash
git add .
git commit -m "feat: xxx"
```

- Execute pre-commit

```bash
pre-commit run -a
```

## License

octopusv is licensed under MIT License.

## Acknowledgments

-

## Anecdote

Agent Octopus Code V: Dive deep into the structural variations' ocean! A state-of-the-art tool equipped with tentacles of precision to merge, correct, and visualize structural results from a multitude of callers

In a world where structural variations are as vast and varied as the ocean, one sentinel rises from the depths to make sense of it all: Agent Octopus, codenamed octopusV.
Imagine a vast digital sea where every wave represents a different structural variation caller, each with its own unique patterns and anomalies.
Navigating through this maze, ensuring seamless merges, and pinpointing corrections used to be daunting tasks. That was until octopusV emerged.
Endowed with eight intelligent tentacles, octopusV has an unparalleled ability to reach out to multiple structural results, embrace them, and bring them together.
Each tentacle is designed to understand, analyze, and adapt to the nuances of different callers, ensuring that the merged structure is flawless.
But the prowess of Agent Octopus doesn’t end there.

Its core, Code V, is a visual marvel.
It's not just correct; it also brings the results to life.
Researchers can now see their structural variation results in vibrant visualizations, making comprehension and further analysis a breeze.

In the expansive sea of structural variations, challenges may be as vast as the ocean, but with octopusV, clarity is but a tentacle touch away.
With Agent Octopus Code V on their side, researchers now had a trusted ally, ensuring excellence in merging, correction, and visualization.
A beacon in the vast ocean of structural variations
