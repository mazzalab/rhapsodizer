import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import argparse
from rhapsodizer import log
from rhapsodizer.r2 import R2
from rhapsodizer.r1 import R1


def main(r1: str, r2: str, r2_map: str, bed: str, read_length: int, stags_file_name: str, index_file_name: str):
    # read alignment passed/dropped read headers
    r2_map_passed: dict
    r2_map_dropped: set
    r2_map_passed, r2_map_dropped = R2.parse_bam(r2_map, bed)

    log.info('Processing R2')
    r2_passed: dict
    r2_dropped: set
    r2_passed, r2_dropped = R2.parse_fastq(r2, stags_file_name, index_file_name)

    # merge dropped reads in R2 and dispose variables
    r2_all_dropped = r2_map_dropped.union(r2_dropped)
    del r2_map_dropped
    del r2_dropped
    
    #r2_all_dropped.clear()  # needed for testing r1.py

    log.info('Processing R1')
    r1_passed, r1_dropped = R1.parse_fastq(r1, read_length, r2_all_dropped)

    # merge dropped reads
    # r2_passed_good = R2.purge_reads(r2_passed, r1_dropped)
    # r1_passed_good = R1.purge_reads(r1_passed, r2_dropped)
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Input and output file arguments')
    parser.add_argument("--r1", help="R1 fastq file")
    parser.add_argument("--r2", help="R2 fastq file")
    parser.add_argument("--r2_bam", help="R2 bam file")
    parser.add_argument("--bed", help="Target genes coordinates")
    parser.add_argument("--tag", help="Sample tags file")
    parser.add_argument("--index", help="Cartridge indices file")
    parser.add_argument("--read_length", help="Sequencing full read length", type=int)
    args = parser.parse_args()

    r1 = args.r1
    r2 = args.r2
    r2_bam = args.r2_bam
    bed = args.bed
    st_f_name = args.tag
    ind_f_name = args.index
    read_length = args.read_length

    log.info('Start processing BD Rhapsody reads')
    main(r1, r2, r2_bam, bed, read_length, st_f_name, ind_f_name)
