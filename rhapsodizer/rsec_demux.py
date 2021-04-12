import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import argparse
import pandas as pd
from rhapsodizer import log
from rhapsodizer.r2 import R2
from rhapsodizer.r1 import R1
from rhapsodizer.matrix import Matrix


def main(r1: str, r2: str, r2_map: str, bed: str, read_length: int, stags_file_name: str, index_file_name: str, matrix_name: str):
    '''
    log.info('Processing R2')
    # read alignment passed/dropped read headers
    r2_map_passed: dict
    r2_map_dropped: set
    r2_map_passed, r2_map_dropped = R2.parse_bam(r2_map, bed)

    r2_passed: dict
    r2_dropped: set
    r2_passed, r2_dropped = R2.parse_fastq(r2, stags_file_name, index_file_name)

    # merge dropped reads in R2 and dispose variables
    r2_all_dropped = r2_map_dropped.union(r2_dropped)
    del r2_map_dropped
    del r2_dropped
    
    log.info('Processing R1')
    r1_passed: dict
    r1_dropped: list
    r1_passed, r1_dropped = R1.parse_fastq(r1, read_length, r2_all_dropped)
    
    # merge dropped reads
    # r2_passed_good = R2.purge_reads(r2_passed, r1_dropped)
    # r1_passed_good = R1.purge_reads(r1_passed, r2_dropped)
    '''
    
    log.info('Generating expression matrix')
    #r1_compiled: dict  = Matrix.compile_r1_passed(r1_passed)
    
    #del r1_passed
    #del r1_dropped
    
    #r1_map_info: dict = Matrix.merge_r1_map_info(r1_compiled, r2_map_passed)
    
    #del r2_map_passed
    
    gene_symbols: list = Matrix.get_gene_symbols(bed)
    
    #rasd_df: pd.DataFrame = Matrix.generate_rasd_matrix(r1_compiled, r1_map_info, gene_symbols)
    # edit file name
    rasd_df: pd.DataFrame = pd.read_csv('../files/cov1-200723_rasd_rebuild.csv', header=0, index_col=0, \
        dtype={'ABCB1': 'uint16', 'ABCB4': 'uint16', 'ADA': 'uint16', 'ADGRE1': 'uint16', 'ADGRG3': 'uint16', 'ADORA1': 'uint16', \
        'AIM2': 'uint16', 'ALAS2': 'uint16', 'ANXA5': 'uint16', 'AOC3': 'uint16', 'APOBEC3G': 'uint16', 'APOE': 'uint16', \
        'AQP4': 'uint16', 'AQP9': 'uint16', 'ARG1': 'uint16', 'ARL4C': 'uint16', 'ARNTL': 'uint16', 'ARNTL2': 'uint16', \
        'ATF6B': 'uint16', 'AURKB': 'uint16', 'AZU1': 'uint16', 'B3GAT1': 'uint16', 'BACH2': 'uint16', 'BAX': 'uint16', \
        'BCL11B': 'uint16', 'BCL2': 'uint16', 'BCL2A1': 'uint16', 'BCL2L11': 'uint16', 'BCL6': 'uint16', 'BHLHE41': 'uint16', \
        'BIN2': 'uint16', 'BIRC3': 'uint16', 'BLK': 'uint16', 'BLNK': 'uint16', 'BPI': 'uint16', 'BTG1': 'uint16', 'BTLA': 'uint16', \
        'C10orf54': 'uint16', 'C1QA': 'uint16', 'C1QB': 'uint16', 'CASP5': 'uint16', 'CBLB': 'uint16', 'CCL1': 'uint16', \
        'CCL13': 'uint16', 'CCL17': 'uint16', 'CCL19': 'uint16', 'CCL2': 'uint16', 'CCL20': 'uint16', 'CCL22': 'uint16', \
        'CCL3': 'uint16', 'CCL4': 'uint16', 'CCL5': 'uint16', 'CCND2': 'uint16', 'CCR1': 'uint16', 'CCR10': 'uint16', 'CCR2': 'uint16', \
        'CCR3': 'uint16', 'CCR4': 'uint16', 'CCR5': 'uint16', 'CCR7': 'uint16', 'CCR8': 'uint16', 'CCR9': 'uint16', 'CD14': 'uint16', \
        'CD160': 'uint16', 'CD163': 'uint16', 'CD1A': 'uint16', 'CD1B': 'uint16', 'CD1C': 'uint16', 'CD2': 'uint16', 'CD200': 'uint16', \
        'CD209': 'uint16', 'CD22': 'uint16', 'CD24': 'uint16', 'CD244': 'uint16', 'CD247': 'uint16', 'CD27': 'uint16', 'CD274': 'uint16', \
        'CD28': 'uint16', 'CD300A': 'uint16', 'CD33': 'uint16', 'CD34': 'uint16', 'CD36': 'uint16', 'CD37': 'uint16', 'CD38': 'uint16', \
        'CD3D': 'uint16', 'CD3E': 'uint16', 'CD3G': 'uint16', 'CD4': 'uint16', 'CD40': 'uint16', 'CD44': 'uint16', 'CD48': 'uint16', \
        'CD5': 'uint16', 'CD52': 'uint16', 'CD6': 'uint16', 'CD63': 'uint16', 'CD69': 'uint16', 'CD7': 'uint16', 'CD70': 'uint16', \
        'CD72': 'uint16', 'CD74': 'uint16', 'CD79A': 'uint16', 'CD79B': 'uint16', 'CD80': 'uint16', 'CD86': 'uint16', 'CD8A': 'uint16', 'CD8B': 'uint16', 'CD9': 'uint16', 'CEACAM8': 'uint16', 'CEBPD': 'uint16', 'CHI3L1': 'uint16', 'CHI3L2': 'uint16', 'CLC': 'uint16', 'CLEC10A': 'uint16', 'CLEC4D': 'uint16', 'CLEC4E': 'uint16', 'CLOCK': 'uint16', 'CMKLR1': 'uint16', 'CMTM2': 'uint16', 'CNOT2': 'uint16', 'CNTNAP3': 'uint16', 'COL1A2': 'uint16', 'CPA3': 'uint16', 'CR2': 'uint16', 'CRY1': 'uint16', 'CRY2': 'uint16', 'CSF2': 'uint16', 'CSF3': 'uint16', 'CST7': 'uint16', 'CTLA4': 'uint16', 'CTSD': 'uint16', 'CTSG': 'uint16', 'CTSW': 'uint16', 'CX3CR1': 'uint16', 'CXCL1': 'uint16', 'CXCL10': 'uint16', 'CXCL11': 'uint16', 'CXCL13': 'uint16', 'CXCL16': 'uint16', 'CXCL2': 'uint16', 'CXCL3': 'uint16', 'CXCL5': 'uint16', 'CXCL8': 'uint16', 'CXCL9': 'uint16', 'CXCR1': 'uint16', 'CXCR2': 'uint16', 'CXCR3': 'uint16', 'CXCR4': 'uint16', 'CXCR5': 'uint16', 'CXCR6': 'uint16', 'DBP': 'uint16', 'DEC1': 'uint16', 'DEFA3': 'uint16', 'DEFA4': 'uint16', 'DOCK8': 'uint16', 'DPP4': 'uint16', 'DUSP1': 'uint16', 'DUSP2': 'uint16', 'DUSP4': 'uint16', 'E': 'uint16', 'EBF1': 'uint16', 'EGR1': 'uint16', 'EGR3': 'uint16', 'ELANE': 'uint16', 'ENTPD1': 'uint16', 'EOMES': 'uint16', 'EPHA1': 'uint16', 'EPX': 'uint16', 'F13A1': 'uint16', 'F5': 'uint16', 'FAM129C': 'uint16', 'FAM65B': 'uint16', 'FAS': 'uint16', 'FASLG': 'uint16', 'FCER1A': 'uint16', 'FCER1G': 'uint16', 'FCER2': 'uint16', 'FCGR3A': 'uint16', 'FCN1': 'uint16', 'FLT3': 'uint16', 'FN1': 'uint16', 'FOSB': 'uint16', 'FOSL1': 'uint16', 'FOXO1': 'uint16', 'FOXP1': 'uint16', 'FOXP3': 'uint16', 'FTH1': 'uint16', 'FUT4': 'uint16', 'FYB': 'uint16', 'FYN': 'uint16', 'GAB2': 'uint16', 'GAPDH': 'uint16', 'GIMAP2': 'uint16', 'GIMAP5': 'uint16', 'GNAI2': 'uint16', 'GNLY': 'uint16', 'GZMA': 'uint16', 'GZMB': 'uint16', 'GZMH': 'uint16', 'GZMK': 'uint16', 'HAVCR2': 'uint16', 'HLA-A': 'uint16', 'HLA-DMA': 'uint16', 'HLA-DPA1': 'uint16', 'HLA-DQB1': 'uint16', 'HLA-DRA': 'uint16', 'HLF': 'uint16', 'HMMR': 'uint16', 'HSP90AA1': 'uint16', 'ICAM1': 'uint16', 'ICOS': 'uint16', 'IER3': 'uint16', 'IFITM2': 'uint16', 'IFITM3': 'uint16', 'IFNA1': 'uint16', 'IFNG': 'uint16', 'IFNGR1': 'uint16', 'IGBP1': 'uint16', 'IGFBP2': 'uint16', 'IGHA1': 'uint16', 'IGHD': 'uint16', 'IGHE': 'uint16', 'IGHG1': 'uint16', 'IGHG2': 'uint16', 'IGHG3': 'uint16', 'IGHG4': 'uint16', 'IGHM': 'uint16', 'IGKC': 'uint16', 'IGLC3': 'uint16', 'IKZF1': 'uint16', 'IKZF2': 'uint16', 'IL12A': 'uint16', 'IL12RB1': 'uint16', 'IL12RB2': 'uint16', 'IL13': 'uint16', 'IL15': 'uint16', 'IL15RA': 'uint16', 'IL17A': 'uint16', 'IL17F': 'uint16', 'IL18': 'uint16', 'IL18R1': 'uint16', 'IL18RAP': 'uint16', 'IL1B': 'uint16', 'IL1R2': 'uint16', 'IL1RL1': 'uint16', 'IL1RN': 'uint16', 'IL2': 'uint16', 'IL21': 'uint16', 'IL22': 'uint16', 'IL23R': 'uint16', 'IL25': 'uint16', 'IL2RA': 'uint16', 'IL2RB': 'uint16', 'IL3': 'uint16', 'IL31': 'uint16', 'IL32': 'uint16', 'IL33': 'uint16', 'IL3RA': 'uint16', 'IL4': 'uint16', 'IL4R': 'uint16', 'IL5': 'uint16', 'IL6': 'uint16', 'IL7R': 'uint16', 'IL9': 'uint16', 'IRF1': 'uint16', 'IRF4': 'uint16', 'IRF6': 'uint16', 'IRF8': 'uint16', 'IRF9': 'uint16', 'ITGA4': 'uint16', 'ITGAE': 'uint16', 'ITGAM': 'uint16', 'ITGAX': 'uint16', 'ITGB2': 'uint16', 'JCHAIN': 'uint16', 'JUN': 'uint16', 'JUNB': 'uint16', 'KCNE3': 'uint16', 'KCNK5': 'uint16', 'KDELR1': 'uint16', 'KIAA0101': 'uint16', 'KIR2DL1': 'uint16', 'KIT': 'uint16', 'KITLG': 'uint16', 'KLRB1': 'uint16', 'KLRC1': 'uint16', 'KLRC3': 'uint16', 'KLRC4': 'uint16', 'KLRF1': 'uint16', 'KLRG1': 'uint16', 'KLRK1': 'uint16', 'LAG3': 'uint16', 'LAIR2': 'uint16', 'LAMP1': 'uint16', 'LAMP3': 'uint16', 'LAP3': 'uint16', 'LAT': 'uint16', 'LAT2': 'uint16', 'LCK': 'uint16', 'LEF1': 'uint16', 'LGALS1': 'uint16', 'LGALS3': 'uint16', 'LGALS9': 'uint16', 'LIF': 'uint16', 'LILRB4': 'uint16', 'LIPA': 'uint16', 'LRRC32': 'uint16', 'LTA': 'uint16', 'LTB': 'uint16', 'LY86': 'uint16', 'LYN': 'uint16', 'M': 'uint16', 'MCM2': 'uint16', 'MCM4': 'uint16', 'MGST1': 'uint16', 'MITF': 'uint16', 'MKI67': 'uint16', 'MME': 'uint16', 'MMP12': 'uint16', 'MMP9': 'uint16', 'MS4A1': 'uint16', 'MTHFR': 'uint16', 'MYC': 'uint16', 'MZB1': 'uint16', 'N': 'uint16', 'NAMPT': 'uint16', 'NCAM1': 'uint16', 'NCR3': 'uint16', 'NFIL3': 'uint16', 'NFKBIA': 'uint16', 'NINJ2': 'uint16', 'NKG7': 'uint16', 'NPAS2': 'uint16', 'NR1D1': 'uint16', 'NR1D2': 'uint16', 'NRP1': 'uint16', 'NT5E': 'uint16', 'ORF10': 'uint16', 'ORF1ab': 'uint16', 'ORF3a': 'uint16', 'ORF6': 'uint16', 'ORF7a': 'uint16', 'ORF7b': 'uint16', 'ORF8': 'uint16', 'PASK': 'uint16', 'PAX5': 'uint16', 'PCNA': 'uint16', 'PDCD1': 'uint16', 'PDIA4': 'uint16', 'PDIA6': 'uint16', 'PER1': 'uint16', 'PER2': 'uint16', 'PER3': 'uint16', 'PI3': 'uint16', 'PIGF': 'uint16', 'PIK3AP1': 'uint16', 'PIK3IP1': 'uint16', 'PMCH': 'uint16', 'POU2AF1': 'uint16', 'PPARGC1B': 'uint16', 'PRDM1': 'uint16', 'PRF1': 'uint16', 'PSEN1': 'uint16', 'PSME2': 'uint16', 'PTGDR2': 'uint16', 'PTPRC': 'uint16', 'PTTG2': 'uint16', 'QPCT': 'uint16', 'RGS1': 'uint16', 'RNASE2': 'uint16', 'RNASE6': 'uint16', 'RORA': 'uint16', 'RORC': 'uint16', 'RPN2': 'uint16', 'RUNX3': 'uint16', 'S': 'uint16', 'S100A10': 'uint16', 'S100A12': 'uint16', 'S100A9': 'uint16', 'SDC4': 'uint16', 'SELL': 'uint16', 'SELPLG': 'uint16', 'SERPINB1': 'uint16', 'SERPINE2': 'uint16', 'SLC25A37': 'uint16', 'SLC7A7': 'uint16', 'SMAD3': 'uint16', 'SMAD4': 'uint16', 'SNCA': 'uint16', 'SOD2': 'uint16', 'SOX9': 'uint16', 'SPP1': 'uint16', 'STAT1': 'uint16', 'STAT3': 'uint16', 'STAT4': 'uint16', 'STAT5A': 'uint16', 'STAT6': 'uint16', 'TARP': 'uint16', 'TBX21': 'uint16', 'TCF4': 'uint16', 'TCF7': 'uint16', 'TCL1A': 'uint16', 'TEF': 'uint16', 'TGFB1': 'uint16', 'TGFB3': 'uint16', 'TGFBI': 'uint16', 'THBD': 'uint16', 'THBS1': 'uint16', 'TIAF1': 'uint16', 'TIGIT': 'uint16', 'TLR2': 'uint16', 'TLR7': 'uint16', 'TLR8': 'uint16', 'TLR9': 'uint16', 'TMEM97': 'uint16', 'TNF': 'uint16', 'TNFRSF13C': 'uint16', 'TNFRSF17': 'uint16', 'TNFRSF25': 'uint16', 'TNFRSF4': 'uint16', 'TNFRSF8': 'uint16', 'TNFRSF9': 'uint16', 'TNFSF10': 'uint16', 'TNFSF13': 'uint16', 'TNFSF13B': 'uint16', 'TNFSF14': 'uint16', 'TNFSF8': 'uint16', 'TOP2A': 'uint16', 'TPSAB1': 'uint16', 'TRAC': 'uint16', 'TRAT1': 'uint16', 'TRBC2': 'uint16', 'TRDC': 'uint16', 'TREM1': 'uint16', 'TRIB2': 'uint16', 'TSPAN32': 'uint16', 'TXK': 'uint16', 'TYMS': 'uint16', 'UBE2C': 'uint16', 'VCAM1': 'uint16', 'VEGFA': 'uint16', 'VMO1': 'uint16', 'VNN2': 'uint16', 'VPREB3': 'uint16', 'VPS28': 'uint16', 'VSIG4': 'uint16', 'XBP1': 'uint16', 'YBX3': 'uint16', 'YY1': 'uint16', 'ZAP70': 'uint16', 'ZBED2': 'uint16', 'ZBTB16': 'uint16', 'ZNF683': 'uint16'})
    
    #del r1_compiled
    #del r1_map_info
    
    rsec_df: pd.DataFrame = Matrix().generate_rsec_matrix(rasd_df, gene_symbols)
    
    del rasd_df
    
    matrix_path = '../files/' + matrix_name + '.csv'
    
    rsec_df.to_csv(matrix_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Input and output file arguments')
    parser.add_argument("--r1", help="R1 fastq file")
    parser.add_argument("--r2", help="R2 fastq file")
    parser.add_argument("--r2_bam", help="R2 bam file")
    parser.add_argument("--bed", help="Target genes coordinates")
    parser.add_argument("--tag", help="Sample tags file")
    parser.add_argument("--index", help="Cartridge indices file")
    parser.add_argument("--read_length", help="Sequencing full read length", type=int)
    parser.add_argument("--matrix", help="Output expression matrix file")
    args = parser.parse_args()

    r1 = args.r1
    r2 = args.r2
    r2_bam = args.r2_bam
    bed = args.bed
    st_f_name = args.tag
    ind_f_name = args.index
    read_length = args.read_length
    matrix_name = args.matrix

    log.info('Start processing BD Rhapsody reads')
    main(r1, r2, r2_bam, bed, read_length, st_f_name, ind_f_name, matrix_name)
