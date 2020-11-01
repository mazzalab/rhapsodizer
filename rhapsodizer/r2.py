""""
File description
"""

import os
from tqdm import tqdm
from gzip import open as gzopen
from rhapsodizer.read import Read

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


class R2(Read):
    min_read_length = 64
    max_snf = 0.80
    min_qual = 20

    @staticmethod
    def read_st(st_f_name: str) -> dict:
        st = {}
        with open(st_f_name, 'r') as st_f:
            for tag in st_f:
                tag_token = tag.split("\t")
                # {SAMPLE_TAG} -> SAMPLE_TAG_NAME
                st[tag_token[0]] = tag_token[1].strip()

        st_f.close()
        return st

    @staticmethod
    def read_cartridge_index(ind_f_name: str) -> dict:
        ind = {}
        with open(ind_f_name, 'r') as ind_f:
            for index in ind_f:
                index_token = index.split("\t")
                # {INDEX} -> {INDEX_NAME}
                ind[index_token[0]] = index_token[1].strip()

        ind_f.close()
        return ind

    @staticmethod
    def get_file_name(file_name: str) -> str:
        # get file name without extension
        base_name = os.path.basename(file_name)
        if base_name.endswith(".fastq.gz"):
            base_out_file_name = base_name.replace('.fastq.gz', '')
        else:
            raise TypeError("The fastq file must be gzipped and the file extension must be '.fastq.gz'")

        return base_out_file_name

    @staticmethod
    def parse(r2_file: str, stags: dict, cart_idx: dict) -> tuple:
        # get file name without extension
        try:
            base_out_file_name = R2.get_file_name(r2_file)
        except TypeError as te:
            raise

        r2_passed = {}
        r2_dropped = []

        with gzopen(r2_file, 'rt') as r2_f:
            r2_f.seek(0, 2)
            uncompressed_size = r2_f.tell()
            r2_f.seek(0)

            with tqdm(total=uncompressed_size) as pbar:
                for header in r2_f:
                    header = header.strip()
                    pbar.update(len(header))

                    # the read index is located in the last 8 character of the header
                    read_idx = header[-8:]
                    if read_idx in cart_idx:
                        read_seq = r2_f.readline().strip()
                        pbar.update(len(read_seq))

                        # Check Read length and Highest Single Nucleotide Frequency (SNF)
                        if not R2.has_minimum_read_length(read_seq) or not R2.check_snf(read_seq):
                            r2_dropped.append((header, "min_length or snf"))
                            r2_f.readline()  # skip "+"
                            temp = r2_f.readline()  # skip qual line
                            pbar.update(len(temp) + 1)
                        else:
                            temp_st_name = [st_names for st, st_names in stags.items() if st in read_seq]

                            r2_f.readline()  # skip "+"
                            read_qual = r2_f.readline().strip()  # process qual line
                            pbar.update(len(read_qual) + 1)
                            if len(temp_st_name) == 1:
                                if not R2.has_minimum_quality_value(read_qual):
                                    r2_dropped.append((header, "min_qual"))
                                else:
                                    r2_passed[header] = base_out_file_name + "_" + temp_st_name[0] + "_" + cart_idx[
                                        read_idx]
                            elif len(temp_st_name) > 1:
                                r2_dropped.append((header, "multiple_st"))
                    else:
                        r2_dropped.append((header, "cart_idx"))
                        temp = len(r2_f.readline())  # skip fasta line
                        r2_f.readline()  # skip "+"
                        temp = temp + len(r2_f.readline().strip())  # skip qual line
                        pbar.update(temp)

                pbar.update(uncompressed_size)
                pbar.set_description(f"Processed {uncompressed_size} reads")
        return r2_passed, r2_dropped
