import numpy as np
import pandas as pd
import logomaker as lm
from Bio import AlignIO
from Bio import SeqIO
from Bio import Align
from Bio.Align import MultipleSeqAlignment
import matplotlib.pyplot as plt

class ShannonEntropy(object):

    """
        sh = ShannonEntropy()

        sh.reaMSA(MSA_FILE = "alignment.fasta")
        
        sh.calculate()

    """
        
    def __init__(self):
        pass

    def readMSA(self, MSA_FILE: str =None, alnformat="fasta", verbose=0):
    
        alignment = Align.read(MSA_FILE, alnformat)

        # Do a little sanity checking:
        seq_lengths_list = [len(record) for record in alignment]
        
        seq_lengths = set(seq_lengths_list)

        if verbose > 0: print("Alignment length is:" + str(list(seq_lengths)))

        if len(seq_lengths) != 1:
            raise ValueError("Your alignment lengths aren't equal. Check your alignment file.")

        index = range(1, list(seq_lengths)[0]+1)

        self.alignment = alignment

        return alignment, list(seq_lengths), index
    
    def __shannon_entropy(self, list_input):
      
        import math
        
        unique_base = set(list_input)
        
        M   =  len(list_input)
        entropy_list = []
        # Number of residues in column
        for base in unique_base:
            n_i = list_input.count(base) # Number of residues of type i
            P_i = n_i/float(M) # n_i(Number of residues of type i) / M(Number of residues in column)
            entropy_i = P_i*(math.log(P_i,2))
            entropy_list.append(entropy_i)

        sh_entropy = -(sum(entropy_list))

        return sh_entropy
    
    def calculate(self):
        
        self.shannon_entropy_list = []
        
        for col_no in range(len(list(self.alignment[0]))):
            list_input = list(self.alignment[:, col_no])
            self.shannon_entropy_list.append(self.__shannon_entropy(list_input))

        return self.shannon_entropy_list 
    
    def __calculate_frequency_matrix(self):
        
        if hasattr(self, "frequency"):
            return self.frequency
        
        counts = {letter: [] for letter in self.alphabet}
        for i in range(len(self.alignment[0])):
            column = [seq[i] for seq in self.alignment]
            for letter in self.alphabet:
                counts[letter].append(column.count(letter) / len(column))
        
        self.counts = counts
        self.frequency = pd.DataFrame(counts)
        
        return self.frequency
    

    def logo_plot(self, alphabet: str = "ACDEFGHIKLMNPQRSTVWY-", cutoff: float = 0.01, top_positions: int = None):
        
        self.alphabet = alphabet

        freq_matrix = self.__calculate_frequency_matrix()
        if top_positions:
            top_positions_indexes = np.argsort(self.shannon_entropy_list)[-top_positions:]
            top_positions_indexes = np.sort(top_positions_indexes)
        else:
            top_positions_indexes = list(np.where(np.array(self.shannon_entropy_list) >= cutoff)[0])
        
        top_freq_matrix = freq_matrix.iloc[top_positions_indexes]
        indexes = list(top_freq_matrix.index)
        top_freq_matrix.reset_index(inplace=True, drop=True)

        lm.Logo(top_freq_matrix, color_scheme="chemistry", figsize=(14, 7), )
        plt.title("Sequence Logo for High-Entropy Positions")
        plt.xlabel("Position")
        plt.ylabel("Frequency")
        plt.xticks(range(len(indexes)), labels=indexes, rotation=45, ha='right')
        plt.savefig("Sequence_Logo_for_High-Entropy_Positions.pdf", format="pdf")
        plt.show()
    
    def bar_plot(self, cutoff: float = 0.01, top_positions: int = None):
        
        if top_positions:
            top_positions_indexes = np.argsort(self.shannon_entropy_list)[-top_positions:]
            top_positions_indexes = np.sort(top_positions_indexes)
        else:
            top_positions_indexes = list(np.where(np.array(self.shannon_entropy_list) >= cutoff)[0])
        
        sh_values = np.array(self.shannon_entropy_list)[top_positions_indexes]
        plt.bar(top_positions_indexes, sh_values, color='blue')
        plt.title("High-Entropy Positions")
        plt.xlabel("Position")
        plt.ylabel("Entropy")
        #plt.xticks(range(len(top_positions_indexes)), labels=top_positions_indexes, rotation=45, ha='right')
        plt.xticks(rotation=45, ha='right')
        plt.savefig("High-Entropy_Positions.pdf", format="pdf")
        plt.show()

