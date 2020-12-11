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
    linker1 = 'ACTGGCCTGCGA'
    linker2 = 'GGTAGCGGTGACA'
    umi_length = 8
    # cl_padding = 2
    # cl1_min_coord = 0
    # cl1_max_coord = 9
    # cl2_min_coord = 21
    # cl2_max_coord = 30
    # cl3_min_coord = 43
    # cl3_max_coord = 52

    def __init__(self, read_length: int):
        self._read_length = read_length

    @classmethod
    def get_linkers_coordinates(cls, read_seq: str) -> tuple or None:
        linker1_match = linker1_coords = linker2_match = linker2_coords = None
        search_rule = "i<=1,d<=1,s<=2,2i+2d+1s<=2:[ATGC]"

        linker1_match = regex.search(rf"(?be)({cls.linker1}){{{search_rule}}}", read_seq)
        if linker1_match:
            linker1_coords = linker1_match.span()
            linker2_match = regex.search(rf"(?be)({cls.linker2}){{{search_rule}}}", read_seq)
            if linker2_match:
                linker2_coords = linker2_match.span()
                return linker1_coords, linker2_coords
            else:
                return None
        else:
            return None

    @classmethod
    def search_umi(cls, read_seq: str, linker2_max_coord: int, last_matched_pos: int) -> UmiData:
        umi_left = linker2_max_coord + last_matched_pos
        return UmiData(umi_seq=read_seq[umi_left:umi_left+cls.umi_length], last_matched_pos=umi_left+cls.umi_length)

    # @classmethod
    # def search_umi(cls, read_seq: str, altered_read: bool, last_matched_pos: int) -> UmiData:
        # if not altered_read:
            # return UmiData(umi_seq=read_seq[52:60], last_matched_pos=60)
        # else:
            # cl3_left = cls.cl3_min_coord - (cls.cl_padding * 4)
            # umi_left = cl3_left+last_matched_pos
            # return UmiData(umi_seq=read_seq[umi_left:umi_left+8], last_matched_pos=umi_left+8)

    @staticmethod
    def cl_match_data(cl_curr: str, cl_num: int, cell_labels: list, search_rule: str,
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
            if cl_num == 1 or cl_num == 2:
                res = None
                for cl_label in cell_labels:
                    res = regex.search(rf"(?be)({cl_label}){{{search_rule}}}", cl_curr)
                    if res:
                        return ClMatchData(orig_seq=cl_label,
                                           matched_seq=cl_curr,
                                           last_matched_pos=res.span()[1])
                return res

            elif cl_num == 3:
                res = None
                for cl_label in cell_labels:
                    res = regex.search(rf"(?be)({cl_curr}){{{search_rule}}}", cl_label)
                    if res:
                        if res.fuzzy_counts[1] == 0 and res.fuzzy_counts[2] == 0:  # substitutions or end indels in cl_curr
                            return ClMatchData(orig_seq=cl_label,
                                               matched_seq=cl_curr,
                                               last_matched_pos=res.span()[1])
                        elif (res.fuzzy_counts[1] == 1 and res.fuzzy_counts[2] == 0) or res.span()[0] == 1:  # mid or start deletion in cl_curr
                            return ClMatchData(orig_seq=cl_label,
                                               matched_seq=cl_curr,
                                               last_matched_pos=res.span()[1]-1)
                        elif res.fuzzy_counts[1] == 0 and res.fuzzy_counts[2] == 1:  # mid or start insertion in cl_curr
                            return ClMatchData(orig_seq=cl_label,
                                               matched_seq=res[0],
                                               last_matched_pos=res.span()[1]+2)
                        else:
                            return None
                return res
            else:
                raise ValueError("Allowed values for the argument 'cl_num' are: 1, 2 or 3.")
                
    # @classmethod
    # def cl_match_data(cls, cl_curr: str, cell_labels: list, matching_rule: str,
                      # fuzzy_only: bool = False) -> ClMatchData or None:
        # """

        # :param fuzzy_only:
        # :param cl_curr:
        # :param cell_labels:
        # :param matching_rule:
        # :return:
        # """
        # if not fuzzy_only and any(cl_curr in s for s in cell_labels):
            # return ClMatchData(orig_seq=cl_curr,
                               # matched_seq=cl_curr,
                               # last_matched_pos=len(cl_curr))
        # else:
            # # Fuzzy search
            # res = None
            # for cl_label in cell_labels:
                # res = regex.search(rf"(?:{cl_label}){{{matching_rule}}}", cl_curr)
                # if res:
                    # return ClMatchData(orig_seq=cl_label,
                                       # matched_seq=cl_curr,
                                       # last_matched_pos=res.regs[0][1])

            # return res
            
    # @classmethod
    # def search_cl(cls, r1_seq: str, cl_num: int, altered_read: bool) -> ClMatchData or None:
        # if altered_read:
            # search_rule = "i<=1,d<=1,s<=3,2i+2d+1s<=4"
            # fuzzy_only = True
        # else:
            # search_rule = "s<=3"
            # fuzzy_only = False

        # if cl_num == 1:
            # cl1_left = cls.cl1_min_coord
            # if altered_read:
                # cl1_right = cls.cl1_max_coord + cls.cl_padding
            # else:
                # cl1_right = cls.cl1_max_coord

            # curr_cl = r1_seq[cl1_left:cl1_right]
            # matched_cl = cls.cl_match_data(curr_cl, cls.cell_labels_1, search_rule, fuzzy_only=fuzzy_only)
        # elif cl_num == 2:
            # if altered_read:
                # cl2_left = cls.cl2_min_coord - (cls.cl_padding * 2)
                # cl2_right = cls.cl2_max_coord + (cls.cl_padding * 2)
            # else:
                # cl2_left = cls.cl2_min_coord
                # cl2_right = cls.cl2_max_coord

            # curr_cl = r1_seq[cl2_left:cl2_right]
            # matched_cl = cls.cl_match_data(curr_cl, cls.cell_labels_2, search_rule, fuzzy_only=fuzzy_only)
        # elif cl_num == 3:
            # if altered_read:
                # cl3_left = cls.cl3_min_coord - (cls.cl_padding * 4)
                # cl3_right = cls.cl3_max_coord + (cls.cl_padding * 4)
            # else:
                # cl3_left = cls.cl3_min_coord
                # cl3_right = cls.cl3_max_coord

            # curr_cl = r1_seq[cl3_left:cl3_right]
            # matched_cl = cls.cl_match_data(curr_cl, cls.cell_labels_3, search_rule, fuzzy_only=fuzzy_only)
        # else:
            # raise ValueError("Allowed values for the argument 'cl_num' are: 1, 2 or 3.")

        # return matched_cl
