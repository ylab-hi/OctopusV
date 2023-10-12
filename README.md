# Octopus: A Tool for Comparing and Merging Structural Variant (SV) Results from Various SV Callers

## What is Octopus?
Octopus is a toolkit for efficient manipulation and comparative analysis of structural variant (SV) callsets from multiple approaches. It streamlines SV analysis workflows by providing a standardized input/output interface and a comprehensive set of commands for filtration, merging, intersection, clustering and visualization of SV calls. Octopus aims to offer a holistic solution for researchers venturing into the complex landscape of structural variations.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites
- Python 3
- Required software/libraries: []

### Installation

```bash
pip install octopus
```

## Usage

### Standardize input with `convert`

Convert VCFs from multiple callers to standardized OCTOPUS format:

```
octopus convert -i caller1.vcf -o caller1.std.vcf 
octopus convert -i caller2.vcf -o caller2.std.vcf
```

### Merge callsets with `merge` 

Extract intersection or union of multiple callsets:

```
octopus merge -i caller1.std.vcf caller2.std.vcf -o intersect.vcf --intersect
octopus merge -i caller1.std.vcf caller2.std.vcf -o union.vcf --union
```

Extract SVs specifically supported by vcf1:

```
octopus merge -i vcf1.std.vcf vcf2.std.vcf -o uniq_vcf1.vcf --specific vcf1.std.vcf
```

Extract SVs specifically supported by vcf1 and vcf2:

```
octopus merge -i vcf1.std.vcf vcf2.std.vcf vcf3.std.vcf -o common_vcf1_vcf2.vcf --specific vcf1.std.vcf vcf2.std.vcf
```

Extract SVs supported by 2+ callers:

```
octopus merge -i vcf1.std.vcf vcf2.std.vcf vcf3.std.vcf -o overlap.vcf --overlap 2
```

### Evaluate calls against truthsets using `bench`

```
octopus bench -i candidate.std.vcf -g truthset.std.vcf -o bench_stat.txt
```

### Generate statistical reports for the called SVs using `stat`

```
octopus stat -i merged.std.vcf -o sv_statistics.txt
```

### Plotting the stat results using `plot`

```
octopus plot -i sv_statistics.vcf -o plot
```

```
octopus plot -i bench_statistics.txt -o plot
```

## Contribution

We welcome contributions from the community. Feel free to submit issues, feature requests, and pull requests.

## License

Octopus is licensed under MIT License.

## Acknowledgments

- 
