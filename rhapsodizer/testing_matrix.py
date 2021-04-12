import pandas as pd
import regex


class Matrix:

    @staticmethod
    def compile_r1_passed(r1_passed: dict) -> dict:
        r1_compiled = {header[1:41] : labels[0] + labels[1] + labels[2] + labels[3] \
                      for header, labels in r1_passed.items()}
                      
        return r1_compiled
        
    @staticmethod
    def merge_r1_map_info(r1_compiled: dict, r2_map_passed: dict) -> dict:
        r1_map_info = {header : [r2_map_passed[header], labels] \
                       for header, labels in r1_compiled.items() \
                       if header in r1_compiled and header in r2_map_passed}
        
        return r1_map_info
    
    @staticmethod
    def get_gene_symbols(bed: str) -> list:
        gene_symbols = set()
        
        with open(bed, 'r') as bedf:
            for line in bedf:
                *_, symbol = line.split('\t')
                gene_symbols.add(symbol.rstrip('\n'))
                
        gene_symbols = sorted(gene_symbols)
        
        return gene_symbols
        
    @staticmethod
    def generate_rasd_matrix(r1_compiled: dict, r1_map_info: dict, gene_symbols: list) -> pd.DataFrame:
        cell_index = set(r1_compiled.values())
        rasd_df = pd.DataFrame(index=cell_index, columns=gene_symbols)
        rasd_df.fillna(0, inplace=True)
        
        for header, info in r1_map_info.items():
            rasd_df.at[info[1], info[0]] += 1
        
        return rasd_df
    
    def collapse_umi(self, cls_umi_dict: dict) -> dict:
        parents_dict = {}  # parent labels : [cumulative frequency, [list of children]]
        matching_rule = 's<=1:[ATGC]'
        
        for cls_umi_1, frequency_1 in cls_umi_dict.items():
            cumulative_frequency = frequency_1
            parent_values = [cumulative_frequency, []]
            res = None
            
            for cls_umi_2, frequency_2 in cls_umi_dict.items():
                if cls_umi_1[:-8] == cls_umi_2[:-8]:
                    res = regex.fullmatch(rf"(?be)({cls_umi_1[-8:]}){{{matching_rule}}}", cls_umi_2[-8:])
                    if res and cls_umi_1[-8:] != cls_umi_2[-8:]:
                        # check parenthood criteria
                        if parent_values[0] >= frequency_2:
                            parent_values[0] += frequency_2
                            parent_values[1].append(cls_umi_2[-8:])
            
            parents_dict[cls_umi_1] = parent_values
        
        purged_parents_dict = {}
        
        for parent in parents_dict.keys():
            if parent not in purged_parents_dict.keys():
                for child in parents_dict.values():
                    is_parent = True
                    if parent[-8:] in child[1]:
                        is_parent = False
                        break
                if is_parent:
                    purged_parents_dict[parent] = parents_dict[parent]
        
        del parents_dict
        
        return purged_parents_dict
    
    def generate_rsec_matrix(self, rasd_df: pd.DataFrame, gene_symbols: list) -> pd.DataFrame:
        rsec_list = []
        rsec_labels = []
        gene_count = len(gene_symbols)
        
        for gene in gene_symbols:
            cls_umi_dict = {}  # labels : frequencies
            gene_pos = gene_symbols.index(gene)
            
            for cls_umi, row in rasd_df.iterrows():
                if row[gene] > 0:
                    cls_umi_dict[cls_umi] = row[gene]
                    
            purged_parents_dict = self.collapse_umi(cls_umi_dict)
            
            for label, frequency in purged_parents_dict.items():
                temp_rsec_list = [0] * gene_count
                temp_rsec_list[gene_pos] = frequency[0]
                rsec_labels.append(label)
                rsec_list.append(temp_rsec_list)
            
        rsec_df = pd.DataFrame(rsec_list, index=rsec_labels, columns=gene_symbols).groupby(level=0).sum()
        
        return rsec_df
        
    @staticmethod
    def generate_dbec_matrix(rsec_df: pd.DataFrame) -> pd.DataFrame:
        pass
