""""
File description
"""

from abc import ABC, abstractmethod
from collections import Counter

__author__ = "Tommaso Mazza"
__copyright__ = "Copyright 2020, The rhapsodizer_project Project"
__version__ = "0.0.1"
__maintainer__ = "Tommaso Mazza"
__email__ = "bioinformatics@css-mendel.it"
__status__ = "Development"
__date__ = "30/10/2020"
__creator__ = "t.mazza"
__license__ = u"""
  Copyright (C) 2020  Tommaso Mazza <t,mazza@css-mendel.it>
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


class Read(ABC):
    min_read_length = 66
    max_snf = 0.55
    min_qual = 20

    @classmethod
    def has_minimum_read_length(cls, read_seq: str) -> bool:
        """
        Read length: If the length of R1 read is <66 or R2 read is <64, the R1/R2 read pair is dropped.
        :param read_seq:
        :return:
        """
        if len(read_seq) < cls.min_read_length:
            return False
        else:
            return True

    @classmethod
    def check_snf(cls, read_seq: str) -> bool:
        """
        Highest Single Nucleotide Frequency (SNF) observed across the bases of the read: If SNF is ≥0.55 for the R1 read
        or SNF ≥0.80 for the R2 read, the read pair is dropped.
        This criterion removes reads with low complexity such as strings of identical bases and tandem repeats.
        :param read_seq:
        :return:
        """
        ncounts = Counter(read_seq)
        ncounts_freq = [x[1] / len(read_seq) for x in ncounts.most_common()]
        if any(i >= cls.max_snf for i in ncounts_freq):
            return False
        else:
            return True

    @classmethod
    def has_minimum_quality_value(cls, read_qual: str) -> bool:
        """
        Mean base quality score of the read: If the mean base quality score of either R1 read or R2 read is <20,
        the read pair is dropped
        :param read_qual:
        :param qual_limit:
        :return:
        """
        read_qual_decoded = [ord(c) - 33 for c in read_qual]
        if sum(read_qual_decoded) / len(read_qual_decoded) < cls.min_qual:
            return False
        else:
            return True

    @staticmethod
    # @abstractmethod
    def parse(read_seq: str):
        pass
