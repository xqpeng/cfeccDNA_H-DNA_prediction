import os
import glob
import argparse

def get_reverse_complement_strictly(seq):
    """
    Computes the reverse complement of a DNA sequence.

    The process follows the standard genomic convention:
    1. Reverses the sequence direction (5' -> 3' becomes 3' -> 5').
    2. Complements each base (A-T, C-G) to represent the opposite strand.

    Args:
        seq (str): Input DNA sequence (supports A, T, C, G, N and lowercase).

    Returns:
        str: The strictly reversed and complemented sequence.
    """
    if not seq: return ""

    reversed_seq = seq[::-1]

    trans_table = str.maketrans("ATCGNatcgn", "TAGCNtagcn")
    final_seq = reversed_seq.translate(trans_table)

    return final_seq


def calculate_metrics(seq):
    """
    Calculates key genomic features for the sequence.

    These metrics (GC content and Purine ratio) are essential for analyzing
    DNA secondary structures, particularly for H-DNA site prediction.

    Returns:
        tuple: (length, gc_content, purine_content)
    """
    length = len(seq)
    if length == 0:
        return 0, 0.0, 0.0

    seq_upper = seq.upper()
    g = seq_upper.count('G')
    c = seq_upper.count('C')
    a = seq_upper.count('A')

    gc_content = ((g + c) / length) * 100

    purine_content = ((a + g) / length) * 100

    return length, gc_content, purine_content


def format_fasta_width(seq, width=80):

    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


def process_line_data(line):
    """
        Parses and merges paired-end read data from a raw data line.

        This function intelligently identifies the positions of R1, R2, and
        the strand orientation (+/-) to reconstruct the full sequence.

        The merging logic: Full_Sequence = R1 + RC(R2)
        """
    parts = line.strip().split()

    if len(parts) < 4:
        return None

    seq_id = parts[0]
    r1_seq = parts[1]

    col2 = parts[2]
    col3 = parts[3]

    r2_seq = ""
    strand = ""

    if col2 in ['+', '-']:
        strand = col2
        r2_seq = col3
    elif col3 in ['+', '-']:
        strand = col3
        r2_seq = col2
    else:
        r2_seq = col2
        strand = col3

    r2_rc = get_reverse_complement_strictly(r2_seq)

    merged_seq = r1_seq + r2_rc

    length, gc, purine = calculate_metrics(merged_seq)

    header = (f">{seq_id}|Length:{length}|Strand:{strand}|"
              f"GC:{gc:.2f}%|Purine:{purine:.2f}%")

    return header, merged_seq


def main():
    parser = argparse.ArgumentParser(description="Convert raw sequence files to FASTA.")
    parser.add_argument("input_dir", help="Path to the directory containing input files.")
    args = parser.parse_args()

    INPUT_FOLDER = args.input_dir
    OUTPUT_FOLDER = os.path.join(INPUT_FOLDER, "output_fasta_final")

    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    all_files = glob.glob(os.path.join(INPUT_FOLDER, "*"))
    input_files = [f for f in all_files if os.path.isfile(f)]

    print(f"目标文件夹: {INPUT_FOLDER}")
    print(f"检测到 {len(input_files)} 个文件，开始处理...\n")

    count = 0
    for file_path in input_files:
        file_name = os.path.basename(file_path)

        if file_name.endswith(".fasta") or file_name.startswith("."):
            continue

        output_path = os.path.join(OUTPUT_FOLDER, f"{file_name}.fasta")

        try:
            with open(file_path, 'r', encoding='utf-8') as f_in, \
                    open(output_path, 'w', encoding='utf-8') as f_out:

                for line in f_in:
                    if not line.strip(): continue

                    try:
                        result = process_line_data(line)
                        if result:
                            header, sequence = result
                            formatted_seq = format_fasta_width(sequence, width=80)

                            f_out.write(f"{header}\n")
                            f_out.write(f"{formatted_seq}\n")

                    except Exception as e:
                        print(f"  [警告] 解析行错误: {e}")

            print(f"已生成: {output_path}")
            count += 1

        except UnicodeDecodeError:
            print(f"  [错误] 文件 {file_name} 编码不是 UTF-8，跳过。")
        except Exception as e:
            print(f"  [错误] 处理文件 {file_name} 失败: {e}")

    print(f"\n处理完成！共处理 {count} 个文件。")
    print(f"结果文件夹: {OUTPUT_FOLDER}")


if __name__ == "__main__":
    main()