# octopusV: Correcting and Standardize benchmarking of Structural Variant (SV) Results from SV Callers

<img src="https://github.com/ylab-hi/octopusV/blob/main/imgs/octopusV_logo.png" width=50% height=50%>

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

## Contribution Guide

We appreciate and welcome contributions from the community. Whether it's submitting issues, feature requests, or pull requests, we're always looking forward to your collaboration.

### Steps to Contribute:

1. **Clone the Repository**

   Begin by cloning the repository to your local machine:

   ```bash
   git clone https://github.com/ylab-hi/octopusV.git
   ```

2. **Create a New Branch**

   Keep your work organized by creating a new branch for each contribution:

   ```bash
   git checkout -b feature/descriptive-branch-name
   ```

3. **Install Dependencies (if required)**

   If you're adding new features or making changes that require additional dependencies, make sure to install them:

   ```bash
   poetry install
   ```

4. **Make Your Changes**

5. **Run pre-commit Hooks**

   ```bash
   pre-commit run -a
   ```

   If all checks pass, proceed to the next step.

6. **Commit and Push Your Changes**

   ```bash
   git add .
   git commit -m "A descriptive message about your changes"
   git push origin feature/descriptive-branch-name
   ```

7. **Create a Pull Request**

   Head over to the `ylab-hi/octopusV` repository on GitHub and click "New Pull Request". Choose your branch from the dropdown and submit your PR with a clear description of the changes you've made.

## License

octopusV is licensed under MIT License.

## Acknowledgments

## Anecdote

**Agent Octopus Code V: The Sentinel of the Structural Variations Sea!**
Plunge deep into the ocean of structural variations with the most advanced tool yet: Agent Octopus Code V. Equipped with precision tentacles, it's designed to merge, correct, and beautifully visualize results from diverse callers.

In the vast realm of structural variations, where each wave symbolizes a unique caller with its distinct characteristics, navigating and deciphering these waves can be overwhelming. Enter OctopusV, the guardian of this intricate digital ocean. Boasting eight versatile tentacles, OctopusV has the unmatched capability to synchronize with a multitude of structural results, encapsulate them, and craft a coherent narrative. Every tentacle is meticulously engineered to grasp, dissect, and harmonize with the intricacies of each caller, ensuring that the final structure stands unblemished.

But the genius of Agent Octopus doesn't halt there. Embedded within its core, Code V, lies a visual spectacle. It's not merely about precision; it's about illuminating your results. Dive into vivid visual representations of your structural variations, simplifying comprehension and paving the way for insightful analysis. In the vast expanse of structural variations, while the challenges are as deep and endless as the sea, with OctopusV, clarity is merely a tentacle's reach away. With the formidable Agent Octopus Code V by their side, researchers have found an unwavering ally, promising impeccable merging, meticulous correction, and dynamic visualization.

OctopusV - Your beacon in the boundless ocean of structural variations.
