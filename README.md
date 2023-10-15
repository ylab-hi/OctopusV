# SVoctopus: A Tool for Comparing and Merging Structural Variant (SV) Results from Various SV Callers

## What is SVoctopus?
SVoctopus is a toolkit for efficient manipulation and comparative analysis of structural variant (SV) callsets from multiple approaches. It streamlines SV analysis workflows by providing a standardized input/output interface and a comprehensive set of commands for filtration, merging, intersection, benchmarking and visualization of VCF from various SV callers. Octopus aims to offer a holistic solution for researchers venturing into the complex landscape of structural variations.

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

Convert VCFs from multiple callers to standardized SVCF format:

```
svoctopus convert -i caller1.vcf -o caller1.svcf 
svoctopus convert -i caller2.vcf -o caller2.svcf
```

### Merge callsets with `merge` 

Extract intersection or union of multiple callsets:

```
svoctopus merge -i caller1.svcf caller2.svcf -o intersect.svcf --intersect
svoctopus merge -i caller1.svcf caller2.svcf -o union.svcf --union
```

Extract SVs specifically supported by vcf1:

```
svoctopus merge -i vcf1.svcf vcf2.svcf -o uniq_vcf1.svcf --specific vcf1.svcf
```

Extract SVs specifically supported by vcf1 and vcf2:

```
svoctopus merge -i vcf1.svcf vcf2.svcf vcf3.svcf -o common_vcf1_vcf2.svcf --specific vcf1.svcf vcf2.svcf
```

Extract SVs supported by 2+ callers:

```
svoctopus merge -i vcf1.svcf vcf2.svcf vcf3.svcf -o overlap.svcf --overlap 2
```

### Evaluate calls against truthsets using `bench`

```
svoctopus bench -i candidate.svcf -g truthset.svcf -o bench_stat.txt
```

### Generate statistical reports for the called SVs using `stat`

```
svoctopus stat -i merged.svcf -o sv_statistics.txt
```

### Plotting the stat results using `plot`

```
svoctopus plot -i sv_statistics.svcf -o plot
```

```
svoctopus plot -i bench_statistics.txt -o plot
```

## Contribution

We welcome contributions from the community. Feel free to submit issues, feature requests, and pull requests.

## License

svoctopus is licensed under MIT License.

## Acknowledgments

- 
