import regex
from typing import NamedTuple


class ClMatchData(NamedTuple):
    """a docstring"""
    orig_seq: str
    matched_seq: str
    last_matched_pos: int


class UmiData(NamedTuple):
    """a docstring"""
    umi_seq: str
    last_matched_pos: int


class Molecule:
    linker1 = ['ACTGGCCTGCGA']
    linker2 = ['GGTAGCGGTGACA']
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
    cl_padding = 2
    cl1_min_coord = 0
    cl1_max_coord = 9
    cl2_min_coord = 21
    cl2_max_coord = 30
    cl3_min_coord = 43
    cl3_max_coord = 52

    def __init__(self, read_length: int):
        self._read_length = read_length

    @classmethod
    def search_umi(cls, read_seq: str, altered_read: bool, last_matched_pos: int) -> UmiData:
        if not altered_read:
            return UmiData(umi_seq=read_seq[52:60], last_matched_pos=60)
        else:
            cl3_left = cls.cl3_min_coord - (cls.cl_padding * 4)
            umi_left = cl3_left+last_matched_pos
            return UmiData(umi_seq=read_seq[umi_left:umi_left+8], last_matched_pos=umi_left+8)

    @classmethod
    def cl_match_data(cls, cl_curr: str, cell_labels: list, matching_rule: str,
                      fuzzy_only: bool = False) -> ClMatchData or None:
        """

        :param fuzzy_only:
        :param cl_curr:
        :param cell_labels:
        :param matching_rule:
        :return:
        """
        if not fuzzy_only and any(cl_curr in s for s in cell_labels):
            return ClMatchData(orig_seq=cl_curr,
                               matched_seq=cl_curr,
                               last_matched_pos=len(cl_curr))
        else:
            # Fuzzy search
            for cl_label in cell_labels:
                res = regex.search(rf"(?:{cl_label}){{{matching_rule}}}", cl_curr)
                if res:
                    return ClMatchData(orig_seq=cl_label,
                                       matched_seq=cl_curr,
                                       last_matched_pos=res.regs[0][1])
                else:
                    return res

    @classmethod
    def search_cl(cls, r1_seq: str, cl_num: int, altered_read: bool) -> ClMatchData or None:
        if altered_read:
            search_rule = "i<=1,d<=1,s<=3,2i+2d+1s<=4"
            fuzzy_only = True
        else:
            search_rule = "s<=3"
            fuzzy_only = False

        if cl_num == 1:
            cl1_left = cls.cl1_min_coord
            if altered_read:
                cl1_right = cls.cl1_max_coord + cls.cl_padding
            else:
                cl1_right = cls.cl1_max_coord

            curr_cl = r1_seq[cl1_left:cl1_right]
            matched_cl = cls.cl_match_data(curr_cl, cls.cell_labels_1, search_rule, fuzzy_only=fuzzy_only)
        elif cl_num == 2:
            if altered_read:
                cl2_left = cls.cl2_min_coord - (cls.cl_padding * 2)
                cl2_right = cls.cl2_max_coord + (cls.cl_padding * 2)
            else:
                cl2_left = cls.cl2_min_coord
                cl2_right = cls.cl2_max_coord

            curr_cl = r1_seq[cl2_left:cl2_right]
            matched_cl = cls.cl_match_data(curr_cl, cls.cell_labels_2, search_rule, fuzzy_only=fuzzy_only)
        elif cl_num == 3:
            if altered_read:
                cl3_left = cls.cl3_min_coord - (cls.cl_padding * 4)
                cl3_right = cls.cl3_max_coord + (cls.cl_padding * 4)
            else:
                cl3_left = cls.cl3_min_coord
                cl3_right = cls.cl3_max_coord

            curr_cl = r1_seq[cl3_left:cl3_right]
            matched_cl = cls.cl_match_data(curr_cl, cls.cell_labels_3, search_rule, fuzzy_only=fuzzy_only)
        else:
            raise ValueError("Allowed values for the argument 'cl_num' are: 1, 2 or 3.")

        return matched_cl
