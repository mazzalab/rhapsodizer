import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import argparse
from rhapsodizer import log
from rhapsodizer.r2 import R2
from rhapsodizer.r1 import R1


def main(r1: str, r2: str, read_length: int, stags_file_name: str, index_file_name: str):
    # read sample tags
    log.info('Reading sample tags file')
    stags = R2.readST(stags_file_name)

    # read cartridge index
    log.info('Reading cardridge index list')
    cart_idx = R2.readCartridgeIndex(index_file_name)

    # process sample tags and cartridge indices
    log.info('Processing R2')
    # r2_passed, r2_dropped = R2.readR2(r2, stags, cart_idx)
    log.info('Processing R1')
    r1_passed, r1_dropped = R1.parse(r1, read_length)
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

    log.info('Start processing Raphsody BD reads')
    main(r1, r2, read_length, st_f_name, ind_f_name)
