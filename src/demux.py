import os

import argparse
from collections import Counter
from gzip import open as gzopen
from annotate_molecule import Molecule


def readST(st_f_name: str) -> dict:
    st = {}
    with open(st_f_name, 'r') as st_f:
        for tag in st_f:
            tag_token = tag.split("\t")
            # {SAMPLE_TAG} -> SAMPLE_TAG_NAME
            st[tag_token[0]] = tag_token[1].strip()

    st_f.close()
    return st


def readCartridgeIndex(ind_f_name: str) -> dict:
    ind = {}
    with open(ind_f_name, 'r') as ind_f:
        for index in ind_f:
            index_token = index.split("\t")
            # {INDEX} -> {INDEX_NAME}
            ind[index_token[0]] = index_token[1].strip()

    ind_f.close()
    return ind


def has_minimum_read_length(read_seq: str, length_limit: int) -> bool:
    """
    Read length: If the length of R1 read is <66 or R2 read is <64, the R1/R2 read pair is dropped.
    :param read_seq:
    :param length_limit:
    :return:
    """
    if len(read_seq) < length_limit:
        return False
    else:
        return True


def check_snf(read_seq: str, snf_limit: float) -> bool:
    """
    Highest Single Nucleotide Frequency (SNF) observed across the bases of the read: If SNF is ≥0.55 for the R1 read
    or SNF ≥0.80 for the R2 read, the read pair is dropped.
    This criterion removes reads with low complexity such as strings of identical bases and tandem repeats.
    :param read_seq:
    :param snf_limit:
    :return:
    """
    ncounts = Counter(read_seq)
    ncounts_freq = [x[1] / len(read_seq) for x in ncounts.most_common()]
    if any(i >= snf_limit for i in ncounts_freq):
        return False
    else:
        return True


def has_minimum_quality_value(read_qual: str, qual_limit: int) -> bool:
    """
    Mean base quality score of the read: If the mean base quality score of either R1 read or R2 read is <20,
    the read pair is dropped
    :param read_qual:
    :param qual_limit:
    :return:
    """
    read_qual_decoded = [ord(c) - 33 for c in read_qual]
    if sum(read_qual_decoded) / len(read_qual_decoded) < qual_limit:
        return False
    else:
        return True


def is_valid_polyt(read_seq_after_umi: str, rate_limit: float) -> bool:
    """
    Each read with a valid cell label is kept for further consideration only if ≥6 out of 8 bases after UMI are
    found to be Ts.
    :param read_seq_after_umi:
    :param rate_limit:
    :return:
    """
    c = Counter(read_seq_after_umi)
    numt = c["T"]
    if numt / len(read_seq_after_umi) >= rate_limit:
        return True
    else:
        return False


def analyze_read(read_seq: str, unaltered_read_length: int) -> tuple:
    cl2 = cl3 = umi = polyt = None

    if len(read_seq) == unaltered_read_length:
        altered_read = False
    else:
        altered_read = True

    cl1 = Molecule.search_cl(read_seq, 1, altered_read)  # 1-9
    if cl1:
        cl2 = Molecule.search_cl(read_seq, 2, altered_read)  # 22–30
        if cl2:
            cl3 = Molecule.search_cl(read_seq, 3, altered_read)  # 44–52
            if cl3:
                umi = Molecule.search_umi(read_seq, altered_read, cl3.last_matched_pos)
                polyt = read_seq[umi.last_matched_pos:]

    return cl1.orig_seq if cl1 else None, \
           cl2.orig_seq if cl2 else None, \
           cl3.orig_seq if cl3 else None, \
           umi.umi_seq if umi else None, \
           polyt


def get_file_name(file_name: str) -> str:
    # get file name without extension
    base_name = os.path.basename(file_name)
    if base_name.endswith(".fastq.gz"):
        base_out_file_name = base_name.replace('.fastq.gz', '')
    else:
        raise TypeError("The fastq file must be gzipped and the file extension must be '.fastq.gz'")

    return base_out_file_name


def readR2(r2_file: str, stags: dict, cart_idx: dict) -> tuple:
    # get file name without extension
    try:
        base_out_file_name = get_file_name(r2_file)
    except TypeError as te:
        raise

    r2_dict = {}
    r2_dropped = []
    with gzopen(r2_file, 'rt') as r2_f:
        for header in r2_f:
            header = header.strip()
            # the read index is located in the last 8 character of the header
            read_idx = header[-8:]
            if read_idx in cart_idx:
                read_seq = r2_f.readline().strip()

                # Check Read length and Highest Single Nucleotide Frequency (SNF)
                if not has_minimum_read_length(read_seq, 64) or not check_snf(read_seq, 0.80):
                    r2_dropped.append(header)
                else:
                    temp_st_name = [st_names for st, st_names in stags.items() if st in read_seq]

                    if len(temp_st_name) == 1:
                        r2_dict[header] = base_out_file_name + "_" + temp_st_name[0] + "_" + cart_idx[read_idx]
                    elif len(temp_st_name) > 1:
                        r2_dropped.append(header)
            else:
                r2_dropped.append(header)
                # skip fasta line
                r2_f.readline()

            # skip "+"
            r2_f.readline()
            # skip quality line
            r2_f.readline()

    return r2_dict, r2_dropped


def readR1(r1_file: str, unaltered_read_length: int) -> tuple:
    r1_passed = {}
    r1_dropped = []
    with gzopen(r1_file, 'rt') as r1_f:
        for header in r1_f:
            header = header.strip()
            read_seq = r1_f.readline().strip()

            # Check Read length and Highest Single Nucleotide Frequency (SNF)
            if not has_minimum_read_length(read_seq, 66) or not check_snf(read_seq, 0.55):
                r1_dropped.append(header)
                # skipp remaining lines concerning this read
                r1_f.readline()
                r1_f.readline()
            else:
                #  skip "+" separator line
                r1_f.readline()

                # Check mean base quality score
                read_qual = r1_f.readline().strip()
                if not has_minimum_quality_value(read_qual, 20):
                    r1_dropped.append(header)
                else:
                    cl1, cl2, cl3, umi, polyt = analyze_read(read_seq, unaltered_read_length)

                    if (cl1 and cl2 and cl3 and umi and polyt) and is_valid_polyt(polyt, 0.75) and "N" not in umi:
                        r1_passed[header] = (cl1, cl2, cl3, umi)
                    else:
                        r1_dropped.append(header)

    return r1_passed, r1_dropped


def main(r1: str, r2: str, read_length: int, stags_file_name: str, index_file_name: str):
    # read sample tags
    stags = readST(stags_file_name)

    # read cartridge index
    cart_idx = readCartridgeIndex(index_file_name)

    # process sample tags and cartridge indices
    # r2_passed, r2_dropped = readR2(r2, stags, cart_idx)
    r1_passed, r1_dropped = readR1(r1, read_length)
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Input and output file arguments')
    parser.add_argument("--r1", help="R1 fastq file")
    parser.add_argument("--r2", help="R2 fastq file")
    parser.add_argument("--tag", help="Sample tags file")
    parser.add_argument("--index", help="Cartridge indices file")
    parser.add_argument("--read_length", help="Sequencing full read length", type=int)
    args = parser.parse_args()

    r1 = args.r1
    r2 = args.r2
    st_f_name = args.tag
    ind_f_name = args.index
    read_length = args.read_length

    from . import log
    log.warning('This is a Warning')
    main(r1, r2, read_length, st_f_name, ind_f_name)
