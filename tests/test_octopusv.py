import subprocess
import pytest
import os
import difflib
import logging

# 配置日志
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# 获取项目根目录的绝对路径
ROOT_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

# 使用绝对路径定义目录结构
TEST_DATA_DIR = os.path.join(ROOT_DIR, "tests", "data")
INPUT_DIR = os.path.join(TEST_DATA_DIR, "input")
STANDARD_DIR = os.path.join(TEST_DATA_DIR, "standard")
OUTPUT_DIR = os.path.join(TEST_DATA_DIR, "output")


def run_octopusv(command, *args, verbose=False):
    """
    运行 octopusv 命令，并捕获输出。
    如果命令失败，则打印详细的错误信息并抛出异常。
    """
    cmd = ["octopusv", command] + list(args)
    if verbose:
        logger.info(f"执行命令: {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=ROOT_DIR
        )
        if verbose and result.stdout:
            logger.info(f"命令 stdout:\n{result.stdout}")
        if result.stderr:
            logger.warning(f"命令 stderr:\n{result.stderr}")
        result.check_returncode()
    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败，退出码 {e.returncode}")
        logger.error(f"错误输出:\n{e.stderr}")
        raise


def compare_files(file1, file2, verbose=False):
    """
    Enhanced file comparison function that handles various edge cases
    and provides detailed difference information.
    """

    def normalize_line(line):
        # 移除常见的路径模式
        import re
        # 移除绝对路径
        line = re.sub(r'/[^,\s;]*/([^/,\s;]+)', r'\1', line)
        # 移除相对路径
        line = re.sub(r'\.\.?/[^,\s;]+/([^/,\s;]+)', r'\1', line)
        return line.rstrip('\r\n')

    def read_file(filename):
        with open(filename, 'r', encoding='utf-8-sig') as f:  # handles BOM
            return [normalize_line(line) for line in f
                    if not line.startswith('##')]

    content1 = read_file(file1)
    content2 = read_file(file2)

    if len(content1) != len(content2):
        if verbose:
            print(f"Files have different number of lines: {len(content1)} vs {len(content2)}")
            print("\nContent of file1:")
            for line in content1:
                print(repr(line))
            print("\nContent of file2:")
            for line in content2:
                print(repr(line))
        return False

    differences = []
    for i, (line1, line2) in enumerate(zip(content1, content2), 1):
        if line1 != line2:
            if verbose:
                print(f"\nDifference at line {i}:")
                print(f"File 1: {repr(line1)}")
                print(f"File 2: {repr(line2)}")

                # Show detailed character comparison
                if len(line1) != len(line2):
                    print(f"Line lengths differ: {len(line1)} vs {len(line2)}")

                # Use difflib to show exact differences
                for i, s in enumerate(difflib.ndiff(line1, line2)):
                    if s[0] == ' ':
                        continue
                    elif s[0] == '-':
                        print(f"Delete '{s[-1]}' from position {i}")
                    elif s[0] == '+':
                        print(f"Add '{s[-1]}' to position {i}")

            differences.append((i, line1, line2))

    return len(differences) == 0


class TestOctopusV:
    def setup_method(self):
        """每个测试运行前的设置方法"""
        # 确保输出目录存在
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    def test_correct_svim(self):
        """测试使用 SVIM 输入的 correct 命令"""
        input_vcf = os.path.join(INPUT_DIR, "svim_example.vcf")
        output_svcf = os.path.join(OUTPUT_DIR, "out_svim_example.svcf")
        standard_svcf = os.path.join(STANDARD_DIR, "standard_svim.svcf")

        # 检查输入和标准文件是否存在
        assert os.path.exists(input_vcf), f"输入文件未找到: {input_vcf}"
        assert os.path.exists(standard_svcf), f"标准文件未找到: {standard_svcf}"

        # 运行 octopusv correct 命令，使用 -i 和 -o 选项
        run_octopusv("correct", "-i", input_vcf, "-o", output_svcf, verbose=False)

        # 检查输出文件是否创建
        assert os.path.exists(output_svcf), f"输出文件未创建: {output_svcf}"

        # 比较输出文件与标准文件的内容
        assert compare_files(output_svcf, standard_svcf, verbose=True), \
            "SVIM correction 输出与标准不匹配"

    def test_correct_pbsv(self):
        """测试使用 PBSV 输入的 correct 命令"""
        input_vcf = os.path.join(INPUT_DIR, "pbsv_example.vcf")
        output_svcf = os.path.join(OUTPUT_DIR, "out_pbsv_example.svcf")
        standard_svcf = os.path.join(STANDARD_DIR, "standard_pbsv.svcf")

        # 检查输入和标准文件是否存在
        assert os.path.exists(input_vcf), f"输入文件未找到: {input_vcf}"
        assert os.path.exists(standard_svcf), f"标准文件未找到: {standard_svcf}"

        # 运行 octopusv correct 命令，使用 -i 和 -o 选项
        run_octopusv("correct", "-i", input_vcf, "-o", output_svcf, verbose=False)

        # 检查输出文件是否创建
        assert os.path.exists(output_svcf), f"输出文件未创建: {output_svcf}"

        # 比较输出文件与标准文件的内容
        assert compare_files(output_svcf, standard_svcf, verbose=True), \
            "PBSV correction 输出与标准不匹配"

    def test_merge_intersect(self):
        """测试 merge 命令的 intersect 选项"""
        input_svcf1 = os.path.join(INPUT_DIR, "pbsv.svcf")
        input_svcf2 = os.path.join(INPUT_DIR, "svim.svcf")
        output_svcf = os.path.join(OUTPUT_DIR, "out_inter.svcf")
        standard_svcf = os.path.join(STANDARD_DIR, "standard_intersection.svcf")

        # 检查输入和标准文件是否存在
        assert os.path.exists(input_svcf1), f"输入文件1未找到: {input_svcf1}"
        assert os.path.exists(input_svcf2), f"输入文件2未找到: {input_svcf2}"
        assert os.path.exists(standard_svcf), f"标准文件未找到: {standard_svcf}"

        # 运行 octopusv merge 命令
        run_octopusv(
            "merge",
            "-i", input_svcf1, input_svcf2,
            "--intersect",
            "-o", output_svcf,
            verbose=True
        )

        # 检查输出文件是否创建
        assert os.path.exists(output_svcf), f"输出文件未创建: {output_svcf}"

        # 使用详细模式比较文件
        is_same = compare_files(output_svcf, standard_svcf, verbose=True)

        # 如果文件不同，打印额外的调试信息
        if not is_same:
            print("\nRaw hex dump of first few lines:")
            for filename in [output_svcf, standard_svcf]:
                print(f"\n{filename}:")
                with open(filename, 'rb') as f:
                    print(f.read().hex())

        assert is_same, "Merge intersect 输出与标准不匹配"