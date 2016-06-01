import random as rnd
import string

import translation
import molecules as mol
import re


class ModelData:
    """
    class to process data sources for usage in the model
    """

    code = dict([('UCA', 'S'), ('UCG', 'S'), ('UCC', 'S'), ('UCU', 'S'),
             ('UUU', 'F'), ('UUC', 'F'), ('UUA', 'L'), ('UUG', 'L'),
             ('UAU', 'Y'), ('UAC', 'Y'), ('UAA', '*'), ('UAG', '*'),
             ('UGU', 'C'), ('UGC', 'C'), ('UGA', '*'), ('UGG', 'W'),
             ('CUA', 'L'), ('CUG', 'L'), ('CUC', 'L'), ('CUU', 'L'),
             ('CCA', 'P'), ('CCG', 'P'), ('CCC', 'P'), ('CCU', 'P'),
             ('CAU', 'H'), ('CAC', 'H'), ('CAA', 'Q'), ('CAG', 'Q'),
             ('CGA', 'R'), ('CGG', 'R'), ('CGC', 'R'), ('CGU', 'R'),
             ('AUU', 'I'), ('AUC', 'I'), ('AUA', 'I'), ('AUG', 'M'),
             ('ACA', 'T'), ('ACG', 'T'), ('ACC', 'T'), ('ACU', 'T'),
             ('AAU', 'N'), ('AAC', 'N'), ('AAA', 'K'), ('AAG', 'K'),
             ('AGU', 'S'), ('AGC', 'S'), ('AGA', 'R'), ('AGG', 'R'),
             ('GUA', 'V'), ('GUG', 'V'), ('GUC', 'V'), ('GUU', 'V'),
             ('GCA', 'A'), ('GCG', 'A'), ('GCC', 'A'), ('GCU', 'A'),
             ('GAU', 'D'), ('GAC', 'D'), ('GAA', 'E'), ('GAG', 'E'),
             ('GGA', 'G'), ('GGG', 'G'), ('GGC', 'G'), ('GGU', 'G')])

    def __init__(self):
        pass

    #do we need this method generally?
    def get_states(self, molecule_class):
        """
        retrieves the information required to construct the different model molecules
        @param molecule_class: BioMolecule class
        @return: list
        """

        if molecule_class == mol.MRNA:
            alphabet = list(self.code.keys())
            mrnas = []
            genes = {}
            for i in range(10):
                sequence = ''.join([rnd.choice(alphabet) for i in range(rnd.randint(50, 500))])
                genes[''.join([rnd.choice(string.ascii_uppercase) for i in range(3)])] = sequence


            for gene in genes:
                for i in range(rnd.randint(1, 10)):
                    mrnas.append(("MRNA_{}_{}".format(gene, i), gene, genes[gene]))
            return mrnas

    def createchromosomes():
    
        chr1=mol.Chromosome(1,"fsa_sequences/S288C_Chromosome I.fsa")
        chr2=mol.Chromosome(2,"fsa_sequences/S288C_Chromosome II.fsa")
        chr3=mol.Chromosome(3,"fsa_sequences/S288C_Chromosome III.fsa")
        chr4=mol.Chromosome(4,"fsa_sequences/S288C_Chromosome IV.fsa")
        chr5=mol.Chromosome(5,"fsa_sequences/S288C_Chromosome V.fsa")
        chr6=mol.Chromosome(6,"fsa_sequences/S288C_Chromosome VI.fsa")
        chr7=mol.Chromosome(7,"fsa_sequences/S288C_Chromosome VII.fsa")
        chr8=mol.Chromosome(8,"fsa_sequences/S288C_Chromosome VIII.fsa")
        chr9=mol.Chromosome(9,"fsa_sequences/S288C_Chromosome IX.fsa")
        chr10=mol.Chromosome(10,"fsa_sequences/S288C_Chromosome X.fsa")
        chr11=mol.Chromosome(11,"fsa_sequences/S288C_Chromosome XI.fsa")
        chr12=mol.Chromosome(12,"fsa_sequences/S288C_Chromosome XII.fsa")
        chr13=mol.Chromosome(13,"fsa_sequences/S288C_Chromosome XIII.fsa")
        chr14=mol.Chromosome(14,"fsa_sequences/S288C_Chromosome XIV.fsa")
        chr15=mol.Chromosome(15,"fsa_sequences/S288C_Chromosome XV.fsa")
        chr16=mol.Chromosome(16,"fsa_sequences/S288C_Chromosome XVI.fsa")
        chrmito=mol.Chromosome(17,"fsa_sequences/S288C_Chromosome Mito.fsa")

        chr_list=[chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chrmito]
    
        return chr_list

    
    def createwholegenome(chr_list):

        whole_genome = ""
        for i in range(len(chr_list)):
            whole_genome += whole_genome + chr_list[i].sequence

        return whole_genome

    def creategenes():

        with open("fsa_sequences/orf_coding.fasta") as orf_fasta:       #open the file, read it and create a list, where each element is a gene with header+sequence
           orf_list = orf_fasta.read().replace("i>", "").replace("sub>", "").replace("->", "").split(">")
    
        orf_list = orf_list[1:] #entfernen des ersten nicht vorhandenden elements
        orf_splitlist = [""]*len(orf_list)  #initialise the new list[gen][header=0 or gen=1]
        
        for i in range(len(orf_list)):  #Trennen von header und sequenz
            orf_splitlist[i] = orf_list[i].split("\n", 1)   

        """
        Gen Sequenz
        """

        for i in range(len(orf_list)):  #replace "\n" 
            orf_splitlist[i][1] = orf_splitlist[i][1].replace("\n", "")


        gene_seq = [x[1] for x in orf_splitlist]
    

        """
        Gen ID
        """

        header_list = [x[0] for x in orf_splitlist]
        header_split = [""]*len(header_list)

        for i in range(len(header_list)):
            header_split[i] = header_list[i].split(" ", 1)


        gene_id = [x[0] for x in header_split]

        """
        Which Chr
        """

        which_chr_list = re.findall(",\s(Chr\s([IVX]+|Mito)|2-micron plasmid)", "".join(header_list)) # find out which chr with regular expressions and put them in a nested list

        which_chr = [x[0] for x in which_chr_list]

        converter_dictionary = {"Chr I" : 1, "Chr II": 2, "Chr III": 3, "Chr IV": 4, "Chr V": 5, "Chr VI": 6, "Chr VII": 7, "Chr VIII": 8, "Chr IX": 9, "Chr X": 10, "Chr XI": 11, "Chr XII": 12, "Chr XIII": 13, "Chr XIV": 14, "Chr XV": 15, "Chr XVI": 16, "Chr Mito": 17, "2-micron plasmid": 18}

        for i in range(len(which_chr)):
            which_chr[i] = converter_dictionary[which_chr[i]]


        """
        gene_name
        """

        gene_name_list = [""]*len(header_list)
   
        for i in range(len(header_split)):
            gene_name_list[i] = header_split[i][1].split(" ", 1)

        gene_name = [x[0] for x in gene_name_list]



        """
        Location
        """

        loc_list = [0]*len(header_list)

        for i in range(len(header_list)):   #creates a list where loc_list[0] = all the tuples of locations for gene 0 in str form
            loc_list[i] = re.findall("([^-][\d]+?)-([\d]+?),", header_list[i])

    
        gene = {}
        for i in range(len(gene_id)):
            gene[gene_id[i]] = mol.Gene(gene_id[i], gene_name[i], which_chr[i], gene_seq[i], loc_list[i])

        return gene
