import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import argparse
from rhapsodizer import log
from rhapsodizer.r2 import R2
from rhapsodizer.r1 import R1
from rhapsodizer.matrix import Matrix


def main(r1: str, r2: str, r2_map: str, bed: str, read_length: int, stags_file_name: str, index_file_name: str):
    log.info('Processing R2')
    # read alignment passed/dropped read headers
    r2_map_passed: dict
    r2_map_dropped: set
    r2_map_passed, r2_map_dropped = R2.parse_bam(r2_map, bed)
    
    ###testing
    #r2_map_passed['NS500299:176:HH2GTBGXG:1:11101:7200:1063'] =  'PTPRC'
    #r2_map_passed['NS500299:176:HH2GTBGXG:1:11101:7201:1063'] =  'PIK3AP1'
    #r2_map_passed['NS500299:176:HH2GTBGXG:1:11101:7202:1063'] =  'PIK3AP1'
    #r2_map_passed['NS500299:176:HH2GTBGXG:1:11101:7203:1063'] =  'PIK3AP1'
    #r2_map_passed['NS500299:176:HH2GTBGXG:1:11101:7204:1063'] =  'SOD2'
    #r2_map_passed['NS500299:176:HH2GTBGXG:1:11101:7205:1063'] =  'SOD2'
    ###
    
    #print('r2_map_passed')
    #for k,v in r2_map_passed.items():
    #    print(k,v)
    
    #print('r2_map_dropped')
    #for x in r2_map_dropped:
    #    print(x)

    r2_passed: dict
    r2_dropped: set
    r2_passed, r2_dropped = R2.parse_fastq(r2, stags_file_name, index_file_name)
    
    #print('r2_passed')
    #for k,v in r2_passed.items():
    #    print(k,v)
    
    #print('r2_dropped')
    #for x in r2_dropped:
    #    print(x)

    # merge dropped reads in R2 and dispose variables
    r2_all_dropped: set = r2_map_dropped.union(r2_dropped)
    del r2_map_dropped
    del r2_dropped
    
    # r2_all_dropped.clear()  # needed for testing r1.py

    log.info('Processing R1')
    r1_passed: dict
    r1_dropped: list
    r1_passed, r1_dropped = R1.parse_fastq(r1, read_length, r2_all_dropped)
    
    ###testing
    #r1_passed['@NS500299:176:HH2GTBGXG:1:11101:7200:1063 1:N:0:CGAGGCTG'] = ('ATCACGTTA', 'CGATGTTTA', 'TTAGGCATA', 'CCCCAGGG')
    #r1_passed['@NS500299:176:HH2GTBGXG:1:11101:7201:1063 1:N:0:CGAGGCTG'] = ('GAGTACATT', 'TACTAGTCA', 'ATCGAGTCT', 'CCCCGGGG')
    #r1_passed['@NS500299:176:HH2GTBGXG:1:11101:7202:1063 1:N:0:CGAGGCTG'] = ('GAGTACATT', 'TACTAGTCA', 'ATCGAGTCT', 'CCCCGGGG')
    #r1_passed['@NS500299:176:HH2GTBGXG:1:11101:7203:1063 1:N:0:CGAGGCTG'] = ('GAGTACATT', 'TACTAGTCA', 'ATCGAGTCT', 'CCCCAGGG')
    #r1_passed['@NS500299:176:HH2GTBGXG:1:11101:7204:1063 1:N:0:CGAGGCTG'] = ('GAGTACATT', 'TACTAGTCA', 'ATCGAGTCT', 'CCCCGGGG')
    #r1_passed['@NS500299:176:HH2GTBGXG:1:11101:7205:1063 1:N:0:CGAGGCTG'] = ('GAGTATTAG', 'TACTAGTCA', 'ATCGAGTCT', 'CCCCGGGG')
    ###
    
    #print('r1_passed')
    #for k,v in r1_passed.items():
    #    print(k,v)
        
    #print('r1_dropped')
    #for x in r1_dropped:
    #    print(x)
    
    log.info('Generating expression matrix')
    r1_compiled: dict  = Matrix.compile_r1_passed(r1_passed)
    
    #print('r1_compiled')
    #for k,v in r1_compiled.items():
    #    print(k,v)
    
    r1_map_info: dict = Matrix.merge_r1_map_info(r1_compiled, r2_map_passed)
    
    #print('r1_map_info')
    #for k,v in r1_map_info.items():
    #    print(k,v)
    
    gene_symbols: list = Matrix.get_gene_symbols(bed)
    
    rasd_df: pd.DataFrame = Matrix.generate_rasd_matrix(r1_compiled, r1_map_info, gene_symbols)
    
    #rasd_df.to_csv('rasd_df')
    #print('rasd_df')
    
    rsec_df: pd.DataFrame = Matrix().generate_rsec_matrix(rasd_df, gene_symbols)
    
    rsec_df.to_csv('rsec_df')
    #print('rsec_df')
    
    # dbec_df: pd.DataFrame = Matrix.generate_dbec_matrix(rsec_df)
    
    # print('dbec_df')
    # dbec_df.to_csv('dbec_df')


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
