# OctopusV: Advanced Structural Variant Analysis Toolkit 🐙

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/octopusV_logo.png" width="50%" height="50%">
</p>

[![PyPI version](https://badge.fury.io/py/octopusv.svg)](https://badge.fury.io/py/octopusv)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

OctopusV is a comprehensive toolkit for analyzing, comparing, and standardizing structural variant (SV) calls from multiple SV callers. It provides a unified framework for SV analysis, offering tools for standardization, merging, benchmarking, and visualization of structural variants.

## Key Features

- **Standardization**: Convert various SV caller outputs to a unified [SVCF format](https://github.com/ylab-hi/octopusV/blob/main/docs/SVCF_specifications.md)
- **Merging**: Combine SV calls from multiple callers with flexible strategies
- **Benchmarking**: Compare SV calls against truth sets using industry-standard metrics
- **Visualization**: Generate informative plots for SV analysis
- **Format Conversion**: Convert between different SV formats (VCF, BED, BEDPE)

## Installation

```bash
pip install octopusv
```

## Quick Start

### 1. Convert SV Calls to Standardized Format

Convert VCF files from different callers to the standardized SVCF format:

```bash
# Basic conversion
octopusv convert -i input.vcf -o output.svcf

# With position tolerance adjustment
octopusv convert -i input.vcf -o output.svcf --pos-tolerance 5
```

### 2. Merge SV Calls

Combine SV calls from multiple callers with various strategies:

```bash
# Intersect calls from multiple callers
octopusv merge -i caller1.svcf caller2.svcf -o merged.svcf --intersect

# Union of calls
octopusv merge -i caller1.svcf caller2.svcf -o merged.svcf --union

# Get SVs supported by at least 2 callers
octopusv merge -i caller1.svcf caller2.svcf caller3.svcf -o merged.svcf --min-support 2

# Extract caller-specific SVs
octopusv merge -i caller1.svcf caller2.svcf -o specific.svcf --specific caller1.svcf

# Complex merging with expression
octopusv merge -i caller1.svcf caller2.svcf caller3.svcf -o merged.svcf \
    --expression "(caller1 AND caller2) AND NOT caller3"

# Generate UpSet plot of intersections
octopusv merge -i caller1.svcf caller2.svcf caller3.svcf -o merged.svcf \
    --intersect --upsetr --upsetr-output intersections.png
```

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/up_upset.png" width="70%" height="70%">
</p>

### 3. Benchmark Against Truth Sets

Evaluate SV calls against a truth set:

```bash
octopusv benchmark \
    truth.svcf \
    test.svcf \
    -o benchmark_results \
    --reference-distance 500 \
    --size-similarity 0.7 \
    --reciprocal-overlap 0.0 \
    --size-min 50 \
    --size-max 50000
```

### 4. Generate Statistics

Analyze SV characteristics with comprehensive statistics and visualizations:

```bash
# Basic statistics output to text file
octopusv stat -i variants.svcf -o stats.txt

# With size filters
octopusv stat -i variants.svcf -o stats.txt --min-size 50 --max-size 10000

# Generate comprehensive HTML report with visualizations
octopusv stat -i variants.svcf -o stats.txt --report

# Generate only plots
octopusv plot stats.txt -o plots
```

The `--report` option generates an interactive HTML report that includes:
- Summary statistics for all SV types
- SV type distribution visualization
- Size distribution analysis
- Chromosome coverage plots
- Quality metrics and genotype information

![HTML Report Example](imgs/html_example.png)


### 5. Create Visualizations

Generate plots from statistics:

```bash
octopusv plot -i stats.txt -o output_prefix
```

This will create:

- `output_prefix_chromosome_distribution.png`: SV distribution across chromosomes

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/chromosome_distribution.png" width="50%" height="50%">
</p>

- `output_prefix_sv_types.png`: Distribution of SV types

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/sv_types.png" width="50%" height="50%">
</p>

- `output_prefix_sv_sizes.png`: SV size distribution

<p align="center">
  <img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/sv_sizes.png" width="50%" height="50%">
</p>

### 6. Format Conversion

Convert between different SV formats:

```bash
# Convert SVCF to BED
octopusv svcf2bed -i variants.svcf -o variants.bed

# Convert SVCF to BEDPE
octopusv svcf2bedpe -i variants.svcf -o variants.bedpe

# Convert SVCF back to VCF
octopusv svcf2vcf -i variants.svcf -o variants.vcf
```

## Advanced Usage

### Merge Configuration Options

- `--max-distance`: Maximum distance between start/end positions (default: 50)
- `--max-length-ratio`: Maximum ratio between event lengths (default: 1.3)
- `--min-jaccard`: Minimum Jaccard index for overlap (default: 0.7)
- `--tra-delta`: Position uncertainty for translocations (default: 50)
- `--tra-min-overlap`: Minimum overlap ratio for translocations (default: 0.5)
- `--tra-strand-consistency`: Require strand consistency (default: True)

### Benchmark Configuration Options

- `--reference-distance`: Maximum reference location distance
- `--size-similarity`: Minimum size similarity ratio
- `--reciprocal-overlap`: Minimum reciprocal overlap
- `--type-ignore`: Allow type mismatches
- `--enable-sequence-comparison`: Enable sequence similarity comparison
- `--sequence-similarity`: Minimum sequence similarity threshold

## Contributing

We welcome contributions! Here's how you can help:

1. Clone the repository:

   ```bash
   git clone https://github.com/ylab-hi/octopusV.git
   ```

2. Create a feature branch:

   ```bash
   git checkout -b feature/your-feature-name
   ```

3. Install development dependencies:

   ```bash
   poetry install
   ```

4. Make your changes and ensure all checks pass:

   ```bash
   pre-commit run -a
   ```

5. Submit a pull request with a clear description of your changes

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions and feedback:

- GitHub Issues: [https://github.com/ylab-hi/octopusV/issues](https://github.com/ylab-hi/octopusV/issues)
- Email: [qingxiang.guo@northwestern.edu]
- Email: [yangyang.li@northwestern.edu]
