""""
File description
"""

from tqdm import tqdm
from gzip import open as gzopen
from collections import Counter
from rhapsodizer.read import Read
from rhapsodizer.annotate_molecule import Molecule

__author__ = "Tommaso Mazza"
__copyright__ = "Copyright 2020, The rhapsodizer_project Project"
__version__ = "0.0.1"
__maintainer__ = "Tommaso Mazza"
__email__ = "bioinformatics@css-mendel.it"
__status__ = "Development"
__date__ = "30/10/2020"
__creator__ = "t.mazza"
__license__ = u"""
  Copyright (C) 2016-2020  Tommaso Mazza <t,mazza@css-mendel.it>
  Viale Regina Margherita 261, 00198 Rome, Italy

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301 USA
  """


class R1(Read):
    cell_labels_1 = ["GTCGCTATA",
                    "CTTGTACTA",
                    "CTTCACATA",
                    "ACACGCCGG",
                    "CGGTCCAGG",
                    "AATCGAATG",
                    "CCTAGTATA",
                    "ATTGGCTAA",
                    "AAGACATGC",
                    "AAGGCGATC",
                    "GTGTCCTTA",
                    "GGATTAGGA",
                    "ATGGATCCA",
                    "ACATAAGCG",
                    "AACTGTATT",
                    "ACCTTGCGG",
                    "CAGGTGTAG",
                    "AGGAGATTA",
                    "GCGATTACA",
                    "ACCGGATAG",
                    "CCACTTGGA",
                    "AGAGAAGTT",
                    "TAAGTTCGA",
                    "ACGGATATT",
                    "TGGCTCAGA",
                    "GAATCTGTA",
                    "ACCAAGGAC",
                    "AGTATCTGT",
                    "CACACACTA",
                    "ATTAAGTGC",
                    "AAGTAACCC",
                    "AAATCCTGT",
                    "CACATTGCA",
                    "GCACTGTCA",
                    "ATACTTAGG",
                    "GCAATCCGA",
                    "ACGCAATCA",
                    "GAGTATTAG",
                    "GACGGATTA",
                    "CAGCTGACA",
                    "CAACATATT",
                    "AACTTCTCC",
                    "CTATGAAAT",
                    "ATTATTACC",
                    "TACCGAGCA",
                    "TCTCTTCAA",
                    "TAAGCGTTA",
                    "GCCTTACAA",
                    "AGCACACAG",
                    "ACAGTTCCG",
                    "AGTAAAGCC",
                    "CAGTTTCAC",
                    "CGTTACTAA",
                    "TTGTTCCAA",
                    "AGAAGCACT",
                    "CAGCAAGAT",
                    "CAAACCGCC",
                    "CTAACTCGC",
                    "AATATTGGG",
                    "AGAACTTCC",
                    "CAAAGGCAC",
                    "AAGCTCAAC",
                    "TCCAGTCGA",
                    "AGCCATCAC",
                    "AACGAGAAG",
                    "CTACAGAAC",
                    "AGAGCTATG",
                    "GAGGATGGA",
                    "TGTACCTTA",
                    "ACACACAAA",
                    "TCAGGAGGA",
                    "GAGGTGCTA",
                    "ACCCTGACC",
                    "ACAAGGATC",
                    "ATCCCGGAG",
                    "TATGTGGCA",
                    "GCTGCCAAT",
                    "ATCAGAGCT",
                    "TCGAAGTGA",
                    "ATAGACGAG",
                    "AGCCCAATC",
                    "CAGAATCGT",
                    "ATCTCCACA",
                    "ACGAAAGGT",
                    "TAGCTTGTA",
                    "ACACGAGAT",
                    "AACCGCCTC",
                    "ATTTAGATG",
                    "CAAGCAAGC",
                    "CAAAGTGTG",
                    "GGCAAGCAA",
                    "GAGCCAATA",
                    "ATGTAATGG",
                    "CCTGAGCAA",
                    "GAGTACATT",
                    "TGCGATCTA",
                    "ATCACGTTA"
                    ]
    cell_labels_2 = ["TACAGGATA",
                    "CACCAGGTA",
                    "TGTGAAGAA",
                    "GATTCATCA",
                    "CACCCAAAG",
                    "CACAAAGGC",
                    "GTGTGTCGA",
                    "CTAGGTCCT",
                    "ACAGTGGTA",
                    "TCGTTAGCA",
                    "AGCGACACC",
                    "AAGCTACTT",
                    "TGTTCTCCA",
                    "ACGCGAAGC",
                    "CAGAAATCG",
                    "ACCAAAATG",
                    "AGTGTTGTC",
                    "TAGGGATAC",
                    "AGGGCTGGT",
                    "TCATCCTAA",
                    "AATCCTGAA",
                    "ATCCTAGGA",
                    "ACGACCACC",
                    "TTCCATTGA",
                    "TAGTCTTGA",
                    "ACTGTTAGA",
                    "ATTCATCGT",
                    "ACTTCGAGC",
                    "TTGCGTACA",
                    "CAGTGCCCG",
                    "GACACTTAA",
                    "AGGAGGCGC",
                    "GCCTGTTCA",
                    "GTACATCTA",
                    "AATCAGTTT",
                    "ACGATGAAT",
                    "TGACAGACA",
                    "ATTAGGCAT",
                    "GGAGTCTAA",
                    "TAGAACACA",
                    "AAATAAATA",
                    "CCGACAAGA",
                    "CACCTACCC",
                    "AAGAGTAGA",
                    "TCATTGAGA",
                    "GACCTTAGA",
                    "CAAGACCTA",
                    "GGAATGATA",
                    "AAACGTACC",
                    "ACTATCCTC",
                    "CCGTATCTA",
                    "ACACATGTC",
                    "TTGGTATGA",
                    "GTGCAGTAA",
                    "AGGATTCAA",
                    "AGAATGGAG",
                    "CTCTCTCAA",
                    "GCTAACTCA",
                    "ATCAACCGA",
                    "ATGAGTTAC",
                    "ACTTGATGA",
                    "ACTTTAACT",
                    "TTGGAGGTA",
                    "GCCAATGTA",
                    "ATCCAACCG",
                    "GATGAACTG",
                    "CCATGCACA",
                    "TAGTGACTA",
                    "AAACTGCGC",
                    "ATTACCAAG",
                    "CACTCGAGA",
                    "AACTCATTG",
                    "CTTGCTTCA",
                    "ACCTGAGTC",
                    "AGGTTCGCT",
                    "AAGGACTAT",
                    "CGTTCGGTA",
                    "AGATAGTTC",
                    "CAATTGATC",
                    "GCATGGCTA",
                    "ACCAGGTGT",
                    "AGCTGCCGT",
                    "TATAGCCCT",
                    "AGAGGACCA",
                    "ACAATATGG",
                    "CAGCACTTC",
                    "CACTTATGT",
                    "AGTGAAAGG",
                    "AACCCTCGG",
                    "AGGCAGCTA",
                    "AACCAAAGT",
                    "GAGTGCGAA",
                    "CGCTAAGCA",
                    "AATTATAAC",
                    "TACTAGTCA",
                    "CAACAACGG",
                    "CGATGTTTA"
                    ]
    cell_labels_3 = ["AAGCCTTCT",
                    "ATCATTCTG",
                    "CACAAGTAT",
                    "ACACCTTAG",
                    "GAACGACAA",
                    "AGTCTGTAC",
                    "AAATTACAG",
                    "GGCTACAGA",
                    "AATGTATCG",
                    "CAAGTAGAA",
                    "GATCTCTTA",
                    "AACAACGCG",
                    "GGTGAGTTA",
                    "CAGGGAGGG",
                    "TCCGTCTTA",
                    "TGCATAGTA",
                    "ACTTACGAT",
                    "TGTATGCGA",
                    "GCTCCTTGA",
                    "GGCACAACA",
                    "CTCAAGACA",
                    "ACGCTGTTG",
                    "ATATTGTAA",
                    "AAGTTTACG",
                    "CAGCCTGGC",
                    "CTATTAGCC",
                    "CAAACGTGG",
                    "AAAGTCATT",
                    "GTCTTGGCA",
                    "GATCAGCGA",
                    "ACATTCGGC",
                    "AGTAATTAG",
                    "TGAAGCCAA",
                    "TCTACGACA",
                    "CATAACGTT",
                    "ATGGGACTC",
                    "GATAGAGGA",
                    "CTACATGCG",
                    "CAACGATCT",
                    "GTTAGCCTA",
                    "AGTTGCATC",
                    "AAGGGAACT",
                    "ACTACATAT",
                    "CTAAGCTTC",
                    "ACGAACCAG",
                    "TACTTCGGA",
                    "AACATCCAT",
                    "AGCCTGGTT",
                    "CAAGTTTCC",
                    "CAGGCATTT",
                    "ACGTGGGAG",
                    "TCTCACGGA",
                    "GCAACATTA",
                    "ATGGTCCGT",
                    "CTATCATGA",
                    "CAATACAAG",
                    "AAAGAGGCC",
                    "GTAGAAGCA",
                    "GCTATGGAA",
                    "ACTCCAGGG",
                    "ACAAGTGCA",
                    "GATGGTCCA",
                    "TCCTCAATA",
                    "AATAAACAA",
                    "CTGTACGGA",
                    "CTAGATAGA",
                    "AGCTATGTG",
                    "AAATGGAGG",
                    "AGCCGCAAG",
                    "ACAGTAAAC",
                    "AACGTGTGA",
                    "ACTGAATTC",
                    "AAGGGTCAG",
                    "TGTCTATCA",
                    "TCAGATTCA",
                    "CACGATCCG",
                    "AACAGAAAC",
                    "CATGAATGA",
                    "CGTACTACG",
                    "TTCAGCTCA",
                    "AAGGCCGCA",
                    "GGTTGGACA",
                    "CGTCTAGGT",
                    "AATTCGGCG",
                    "CAACCTCCA",
                    "CAATAGGGT",
                    "ACAGGCTCC",
                    "ACAACTAGT",
                    "AGTTGTTCT",
                    "AATTACCGG",
                    "ACAAACTTT",
                    "TCTCGGTTA",
                    "ACTAGACCG",
                    "ACTCATACG",
                    "ATCGAGTCT",
                    "CATAGGTCA",
                    "TTAGGCATA"
                    ]
    min_read_length = 66
    max_snf = 0.55
    min_qual = 20
    min_t_in_polyt = 6
    cl_length = 9
    umi_length = 8
    cl1_min_coord = 0
    cl1_max_coord = 9
    cl2_min_coord = 21
    cl2_max_coord = 30
    cl3_min_coord = 43
    cl3_max_coord = 52
    umi_max_coord = 60

    @classmethod
    def is_valid_polyt(cls, read_seq_after_umi: str) -> bool:
        """
        Each read with a valid cell label is kept for further consideration only if ≥6 out of 8 bases after UMI are
        found to be Ts.
        :param read_seq_after_umi:
        :return:
        """
        c = Counter(read_seq_after_umi[0:9])
        numt = c["T"]
        return numt >= cls.min_t_in_polyt

    '''
    @classmethod
    def analyze_exact_read(cls, read_seq: str) -> tuple:
        cl1 = cl2 = cl3 = umi = polyt = None
        cl1 = read_seq[cls.cl1_min_coord:cls.cl1_max_coord]
        cl2 = read_seq[cls.cl2_min_coord:cls.cl2_max_coord]
        cl3 = read_seq[cls.cl3_min_coord:cls.cl3_max_coord]
        umi = read_seq[cls.cl3_max_coord:cls.umi_max_coord]
        polyt = read_seq[cls.umi_max_coord:]
        
        return cl1 if cl1 in cls.cell_labels_1 else None, \
               cl2 if cl2 in cls.cell_labels_2 else None, \
               cl3 if cl3 in cls.cell_labels_2 else None, \
               umi if umi else None, \
               polyt
               
    '''
    @classmethod
    def analyze_read(cls, read_seq: str) -> tuple:
        linkers_coords = cl1_matched = cl2_matched = cl3_matched = umi = polyt = None
        linkers_coords = Molecule.get_linkers_coordinates(read_seq)  # ((min1, max1), (min2, max2))

        if linkers_coords:
            cl1 = read_seq[:linkers_coords[0][0]]

            if cls.cl_length-1 <= len(cl1) <= cls.cl_length+1:  # only one indel is allowed; update if search_rule gets modified
                if len(cl1) == cls.cl_length:
                    search_rule = "s<=1:[ATGC]"
                    fuzzy_only = False
                else:
                    search_rule = "i<=1,d<=1,s<=1,2i+2d+1s<=2:[ATGC]"
                    fuzzy_only = True

                cl1_matched = Molecule.cl_match_data(cl1, 1, cls.cell_labels_1, search_rule, fuzzy_only=fuzzy_only)
                if cl1_matched:
                    cl2 = read_seq[linkers_coords[0][1]:linkers_coords[1][0]]

                    if cls.cl_length-1 <= len(cl2) <= cls.cl_length+1:
                        if len(cl2) == cls.cl_length:
                            search_rule = "s<=1:[ATGC]"
                            fuzzy_only = False
                        else:
                            search_rule = "i<=1,d<=1,s<=1,2i+2d+1s<=2:[ATGC]"
                            fuzzy_only = True

                        cl2_matched = Molecule.cl_match_data(cl2, 2, cls.cell_labels_2, search_rule, fuzzy_only=fuzzy_only)
                        if cl2_matched:
                            cl3 = read_seq[linkers_coords[1][1]:linkers_coords[1][1]+cls.cl_length]

                            # cl3 is initialized having len == 9
                            search_rule = "i<=1,d<=1,s<=1,2i+2d+1s<=2:[ATGC]"
                            fuzzy_only = True

                            cl3_matched = Molecule.cl_match_data(cl3, 3, cls.cell_labels_3, search_rule, fuzzy_only=fuzzy_only)
                            if cl3_matched:
                                umi = Molecule.search_umi(read_seq, linkers_coords[1][1], cl3_matched.last_matched_pos)
                                polyt = read_seq[umi.last_matched_pos:]

        return cl1_matched.orig_seq if cl1_matched else None, \
               cl2_matched.orig_seq if cl2_matched else None, \
               cl3_matched.orig_seq if cl3_matched else None, \
               umi.umi_seq if umi else None, \
               polyt

    # @staticmethod
    # def analyze_read(read_seq: str, unaltered_read_length: int) -> tuple:
        # cl2 = cl3 = umi = polyt = None

        # if len(read_seq) == unaltered_read_length:
            # altered_read = False
        # else:
            # altered_read = True

        # cl1 = Molecule.search_cl(read_seq, 1, altered_read)  # 1-9
        # if cl1:
            # cl2 = Molecule.search_cl(read_seq, 2, altered_read)  # 22–30
            # if cl2:
                # cl3 = Molecule.search_cl(read_seq, 3, altered_read)  # 44–52
                # if cl3:
                    # umi = Molecule.search_umi(read_seq, altered_read, cl3.last_matched_pos)
                    # polyt = read_seq[umi.last_matched_pos:]

        # return cl1.orig_seq if cl1 else None, \
               # cl2.orig_seq if cl2 else None, \
               # cl3.orig_seq if cl3 else None, \
               # umi.umi_seq if umi else None, \
               # polyt

    @staticmethod
    def parse_fastq(r1_file: str, unaltered_read_length: int, r2_dropped: set = None) -> tuple:
        r1_passed = {}
        r1_dropped = []

        with gzopen(r1_file, 'rt') as r1_f:
            r1_f.seek(0, 2)
            uncompressed_size = r1_f.tell()
            r1_f.seek(0)

            with tqdm(total=uncompressed_size) as pbar:
                for header in r1_f:
                    header = header.strip()
                    pbar.update(len(header))

                    # continue if this read's mate has been dropped in R2
                    temp = header[:42] + "2" + header[43:]
                    if temp in r2_dropped:
                        # skip following 3 lines
                        temp = len(header) + len(r1_f.readline())  # skip seq
                        r1_f.readline()  # skip '+'
                        temp = temp + len(r1_f.readline() + 1)  # skip qual
                        pbar.update(temp)
                        continue

                    # Check Read length and Highest Single Nucleotide Frequency (SNF)
                    read_seq = r1_f.readline().strip()
                    pbar.update(len(read_seq))
                    if not R1.has_minimum_read_length(read_seq) or not R1.check_snf(read_seq):
                        r1_dropped.append((header, "short", "snf"))
                        # skip remaining info belonging to this read
                        r1_f.readline()

                        temp = r1_f.readline()
                        pbar.update(len(temp) + 1)  # +1 accounts for "+" of the previous line
                    else:
                        # skip "+" separator line
                        r1_f.readline()

                        # Check mean base quality score
                        read_qual = r1_f.readline().strip()
                        pbar.update(len(read_qual) + 1)  # +1 accounts for "+" of the previous line
                        if not R1.has_minimum_quality_value(read_qual):
                            r1_dropped.append((header, "min_qual"))
                        else:
                            cl1, cl2, cl3, umi, polyt = R1.analyze_read(read_seq)

                            if (cl1 and cl2 and cl3 and umi and polyt) and R1.is_valid_polyt(polyt) and "N" not in umi:
                                r1_passed[header] = (cl1, cl2, cl3, umi)
                                # write to file the result
                                #with open('../files/r1_session_1Mtest.txt', 'a') as r1f:
                                    #r1f.write(f'{header} {r1_passed[header]}\n')
                            else:
                                r1_dropped.append((header, cl1, cl2, cl3, umi))

                pbar.update(uncompressed_size)
                pbar.set_description(f"Processed {uncompressed_size} reads")
        return r1_passed, r1_dropped

    # @staticmethod
    # def parse_fastq(r1_file: str, unaltered_read_length: int, r2_dropped: set = None) -> tuple:
        # r1_passed = {}
        # r1_dropped = []

        # with gzopen(r1_file, 'rt') as r1_f:
            # r1_f.seek(0, 2)
            # uncompressed_size = r1_f.tell()
            # r1_f.seek(0)

            # with tqdm(total=uncompressed_size) as pbar:
                # for header in r1_f:
                    # header = header.strip()
                    # pbar.update(len(header))

                    # # continue if this read's mate has been dropped in R2
                    # temp = header[:42] + "2" + header[43:]
                    # if temp in r2_dropped:
                        # # skip following 3 lines
                        # temp = len(header) + len(r1_f.readline())  # skip seq
                        # r1_f.readline()  # skip '+'
                        # temp = temp + len(r1_f.readline() + 1)  # skip qual
                        # pbar.update(temp)
                        # continue

                    # # Check Read length and Highest Single Nucleotide Frequency (SNF)
                    # read_seq = r1_f.readline().strip()
                    # pbar.update(len(read_seq))
                    # if not R1.has_minimum_read_length(read_seq) or not R1.check_snf(read_seq):
                        # r1_dropped.append((header, "short", "snf"))
                        # # skip remaining info belonging to this read
                        # r1_f.readline()

                        # temp = r1_f.readline()
                        # pbar.update(len(temp) + 1)  # +1 accounts for + of the previous line
                    # else:
                        # # skip "+" separator line
                        # r1_f.readline()

                        # # Check mean base quality score
                        # read_qual = r1_f.readline().strip()
                        # pbar.update(len(read_qual) + 1)  # +1 accounts for + of the previous line
                        # if not R1.has_minimum_quality_value(read_qual):
                            # r1_dropped.append((header, "min_qual"))
                        # else:
                            # cl1, cl2, cl3, umi, polyt = R1.analyze_read(read_seq, unaltered_read_length)

                            # if (cl1 and cl2 and cl3 and umi and polyt) and R1.is_valid_polyt(polyt) and "N" not in umi:
                                # r1_passed[header] = (cl1, cl2, cl3, umi)
                            # else:
                                # r1_dropped.append((header, cl1, cl2, cl3, umi))

                # pbar.update(uncompressed_size)
                # pbar.set_description(f"Processed {uncompressed_size} reads")
        # return r1_passed, r1_dropped
