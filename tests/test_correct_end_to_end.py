import subprocess
import filecmp
import pytest
import os

# The test dir
TEST_DATA_DIR = "tests/data/VCF_for_testing_correct"
OUTPUT_DIR = "tests/output"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def run_octopusv(input_file, output_file):
    """
    Run octopusv.
    """
    cmd = f"octopusv correct {input_file} {output_file}"
    subprocess.run(cmd, shell=True, check=True)


def test_mate_bnd_independent():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_mate_bnd_independent.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_mate_bnd_independent.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_mate_bnd_independent.vcf")

    run_octopusv(input_vcf, output_vcf)

    # Compare output and expected file
    assert filecmp.cmp(output_vcf, expected_output_vcf)


def test_mate_bnd_merge():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_mate_bnd_merge.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_mate_bnd_merge.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_mate_bnd_merge.vcf")

    run_octopusv(input_vcf, output_vcf)

    assert filecmp.cmp(output_vcf, expected_output_vcf)


def test_mate_bnd_reciprocal():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_mate_bnd_reciprocal.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_mate_bnd_reciprocal.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_mate_bnd_reciprocal.vcf")

    run_octopusv(input_vcf, output_vcf)

    # Compare output and expected file
    assert filecmp.cmp(output_vcf, expected_output_vcf)


def test_non_bnd_events():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_non_bnd_events.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_non_bnd_events.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_non_bnd_events.vcf")

    run_octopusv(input_vcf, output_vcf)

    # Compare output and expected file
    assert filecmp.cmp(output_vcf, expected_output_vcf)


def test_same_bnd_to_dup():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_same_bnd_to_dup.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_same_bnd_to_dup.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_same_bnd_to_dup.vcf")

    run_octopusv(input_vcf, output_vcf)

    # Compare output and expected file
    assert filecmp.cmp(output_vcf, expected_output_vcf)


def test_same_bnd_to_FTRA():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_same_bnd_to_FTRA.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_same_bnd_to_FTRA.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_same_bnd_to_FTRA.vcf")

    run_octopusv(input_vcf, output_vcf)

    # Compare output and expected file
    assert filecmp.cmp(output_vcf, expected_output_vcf)


def test_same_bnd_to_INV():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_same_bnd_to_INV.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_same_bnd_to_INV.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_same_bnd_to_INV.vcf")

    run_octopusv(input_vcf, output_vcf)

    # Compare output and expected file
    assert filecmp.cmp(output_vcf, expected_output_vcf)


def test_same_bnd_to_RTRA():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_same_bnd_to_RTRA.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_same_bnd_to_RTRA.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_same_bnd_to_RTRA.vcf")

    run_octopusv(input_vcf, output_vcf)

    # Compare output and expected file
    assert filecmp.cmp(output_vcf, expected_output_vcf)


def test_single_TRA_to_TRA():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_single_TRA_to_TRA.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_single_TRA_to_TRA.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_single_TRA_to_TRA.vcf")

    run_octopusv(input_vcf, output_vcf)

    # Compare output and expected file
    assert filecmp.cmp(output_vcf, expected_output_vcf)


def test_special_no_mate_diff_bnd_pair_to_independent():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_special_no_mate_diff_bnd_pair_to_independent.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_special_no_mate_diff_bnd_pair_to_independent.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_special_no_mate_diff_bnd_pair_to_independent.vcf")

    run_octopusv(input_vcf, output_vcf)

    # Compare output and expected file
    assert filecmp.cmp(output_vcf, expected_output_vcf)


def test_special_no_mate_diff_bnd_pair_to_reciprocal():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_special_no_mate_diff_bnd_pair_to_reciprocal.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_special_no_mate_diff_bnd_pair_to_reciprocal.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_special_no_mate_diff_bnd_pair_to_reciprocal.vcf")

    run_octopusv(input_vcf, output_vcf)

    # Compare output and expected file
    assert filecmp.cmp(output_vcf, expected_output_vcf)
