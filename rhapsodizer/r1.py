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
    min_read_length = 66
    max_snf = 0.55
    min_qual = 20
    min_t_in_polyt = 6

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
        if numt >= cls.min_t_in_polyt:
            return True
        else:
            return False

    @staticmethod
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
                        pbar.update(len(temp) + 1)  # +1 accounts for + of the previous line
                    else:
                        # skip "+" separator line
                        r1_f.readline()

                        # Check mean base quality score
                        read_qual = r1_f.readline().strip()
                        pbar.update(len(read_qual) + 1)  # +1 accounts for + of the previous line
                        if not R1.has_minimum_quality_value(read_qual):
                            r1_dropped.append((header, "min_qual"))
                        else:
                            cl1, cl2, cl3, umi, polyt = R1.analyze_read(read_seq, unaltered_read_length)

                            if (cl1 and cl2 and cl3 and umi and polyt) and R1.is_valid_polyt(polyt) and "N" not in umi:
                                r1_passed[header] = (cl1, cl2, cl3, umi)
                            else:
                                r1_dropped.append((header, cl1, cl2, cl3, umi))

                pbar.update(uncompressed_size)
                pbar.set_description(f"Processed {uncompressed_size} reads")
        return r1_passed, r1_dropped
