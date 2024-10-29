import subprocess
import filecmp
import pytest
import os

# Define directory structure
TEST_DATA_DIR = os.path.join("tests", "data")
INPUT_DIR = os.path.join(TEST_DATA_DIR, "input")
STANDARD_DIR = os.path.join(TEST_DATA_DIR, "standard")
OUTPUT_DIR = os.path.join("tests", "output")

# Create directories if they don't exist
os.makedirs(INPUT_DIR, exist_ok=True)
os.makedirs(STANDARD_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)


def run_octopusv(command, *args):
    """
    Run octopusv with specified command and arguments.
    """
    cmd = ["octopusv", command] + list(args)
    subprocess.run(cmd, check=True)


class TestOctopusV:
    def test_correct_svim(self):
        """Test correct command with SVIM input"""
        input_vcf = os.path.join(INPUT_DIR, "svim_example.vcf")
        output_svcf = os.path.join(OUTPUT_DIR, "out_svim_example.svcf")
        standard_svcf = os.path.join(STANDARD_DIR, "standard_svim_example.svcf")

        run_octopusv("correct", "-i", input_vcf, "-o", output_svcf)

        assert filecmp.cmp(output_svcf, standard_svcf), \
            "SVIM correction output doesn't match standard"

    def test_correct_pbsv(self):
        """Test correct command with PBSV input"""
        input_vcf = os.path.join(INPUT_DIR, "pbsv_example.vcf")
        output_svcf = os.path.join(OUTPUT_DIR, "out_pbsv_example.svcf")
        standard_svcf = os.path.join(STANDARD_DIR, "standard_pbsv_example.svcf")

        run_octopusv("correct", "-i", input_vcf, "-o", output_svcf)

        assert filecmp.cmp(output_svcf, standard_svcf), \
            "PBSV correction output doesn't match standard"

    def test_merge_intersect(self):
        """Test merge command with intersect option"""
        input_svcf1 = os.path.join(STANDARD_DIR, "standard_pbsv_example.svcf")
        input_svcf2 = os.path.join(STANDARD_DIR, "standard_svim_example.svcf")
        output_svcf = os.path.join(OUTPUT_DIR, "out_inter.svcf")
        standard_svcf = os.path.join(STANDARD_DIR, "standard_inter.svcf")

        run_octopusv("merge", "-i", input_svcf1, input_svcf2,
                    "--intersect", "-o", output_svcf)

        assert filecmp.cmp(output_svcf, standard_svcf), \
            "Merge intersect output doesn't match standard"

    def setup_method(self):
        """Setup method that runs before each test"""
        # Clear output directory before each test
        for file in os.listdir(OUTPUT_DIR):
            os.remove(os.path.join(OUTPUT_DIR, file))