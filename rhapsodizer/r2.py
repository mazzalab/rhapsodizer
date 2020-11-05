""""
File description
"""

import os
import bamnostic as bs
from tqdm import tqdm
from gzip import open as gzopen
from rhapsodizer import log
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


class R2(Read):
    min_read_length = 64
    max_snf = 0.80
    min_qual = 20
    min_mapping_qual = 4
    priming_window = 5
    total_cigar_m = 60

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

    @classmethod
    def parse_bam(cls, r2_bam: str, bed: str) -> tuple:
        """
        An R2 read is a valid gene alignment if all of these criteria are met:
        The read aligns uniquely to a transcript sequence in the reference.
        The R2 alignment begins within the first five nucleotides. This criterion
        ensures that the R2 read originates from an actual PCR priming event.
        The length of the alignment that can be a match or mismatch in the CIGAR
        string is >60.
        The read does not align to phiX174.
        :param r2_bam:
        :param bed:
        :param min_mapping_qual:
        :param priming_window:
        :param total_cigar_m:
        :return:
        """
        r2_map_passed = dict()
        r2_map_dropped = set()
        
        # read bam file
        log.info('Processing BAM file')
        bam = bs.AlignmentFile(r2_bam, 'rb')
        
        for read in bam:
            
            # check if read is uniquely mapped
            if read.mapping_quality >= cls.min_mapping_qual:
            
                # check if priming occurs in the first n nucleotides
                nt = cls.priming_window
                while nt > 0:
                    for operator in read.cigar:
                        if operator[0] == 0:
                            priming = True
                            break  # skip remaining cigar operators
                        else:
                            nt = nt - operator[1]
                    break  # exit while loop
                else:
                    r2_map_dropped.add(read.read_name)
                    priming = False
                
                # check if the total CIGAR M-operation is > m
                if priming:
                    cigar_dict = dict()
                    for n,m in read.cigar:
                        cigar_dict.setdefault(n, []).append(m)
                    if sum(cigar_dict[0]) > cls.total_cigar_m:
                    
                        # the read is a valid gene alignment if at least 1 nt is overlapping
                        # should be using query_length, query_alignment_length or reference_length?
                        
                        # read bed file
                        with open(bed, 'r') as bedf:
                            read_start = read.pos + 1  # 1-based transcription start
                            read_end = read.pos + 1 + read.query_length  # 1-based transcription end
                            for line in bedf:
                                gene_pos, gene_start, gene_end, gene_symbol = line.split('\t')
                                
                                # check if read maps in the chromosome of the current gene coordinates
                                if read.reference_name == gene_pos:
                                    
                                    # check in which gene the read aligns
                                    if (int(gene_start) <= read_start < int(gene_end)) \
                                    or (read_start < int(gene_start) <= read_end) \
                                    or (read_start <= int(gene_end) < read_end):
                                        r2_map_passed[read.read_name] = gene_symbol.rstrip('\n')
                                        break  # skip remaining lines of bed file
                                    else:
                                        continue  # read next gene coordinates
                                    
                                else:
                                    continue  # read next gene coordinates
                            
                            # drop good reads not mapping in any of the given genes
                            if read.read_name not in r2_map_passed.keys():
                                r2_map_dropped.add(read.read_name)
                    
                    else:
                        r2_map_dropped.add(read.read_name)
                        
                else:
                    pass
            
            else:
                r2_map_dropped.add(read.read_name)

        bam.close()

        return r2_map_passed, r2_map_dropped

    @staticmethod
    def parse_fastq(r2_file: str, stags_file_name: str, index_file_name: str) -> tuple:
        # get file name without extension
        try:
            base_out_file_name = R2.get_file_name(r2_file)
        except TypeError as te:
            raise

        # read sample tags
        log.info('Reading sample tags file')
        stags: dict = R2.read_st(stags_file_name)

        # read cartridge index
        log.info('Reading cartridge index list')
        cart_idx: dict = R2.read_cartridge_index(index_file_name)

        r2_passed = {}
        r2_dropped = set()
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
                            r2_dropped.add((header, "min_length or snf"))
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
                                    r2_dropped.add((header, "min_qual"))
                                else:
                                    r2_passed[header] = base_out_file_name + "_" + temp_st_name[0] + "_" + cart_idx[
                                        read_idx]
                            elif len(temp_st_name) > 1:
                                r2_dropped.add((header, "multiple_st"))
                    else:
                        r2_dropped.add((header, "cart_idx"))
                        temp = len(r2_f.readline())  # skip fasta line
                        r2_f.readline()  # skip "+"
                        temp = temp + len(r2_f.readline().strip()) + 1  # skip qual line
                        pbar.update(temp)

                pbar.update(uncompressed_size)
                pbar.set_description(f"Processed {uncompressed_size} reads")
        return r2_passed, r2_dropped
