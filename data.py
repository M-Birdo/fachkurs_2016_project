import re
import sys



class Chromosome:
    """
    An object for the chromosome with the attributes:

    id -> gives out the chromosome in integer form 1-16, 17 = mitochondrial chromosome
    fastaname -> is only important for the initiation of the object and used by the function createchromosomes
    sequence -> sequence of the whole chromosome
    revsequence -> the reverse sequence of the whole chromosome

    possible operations are:
    gene.id -> gives the chromosome number 1-16, 17 = mitochondrial chromosome
    gene.sequence -> gives the sequence of the chromosome
    gene.revsequence -> gives out the reverse sequence of the chromosome
    """
    def __init__(self, id ,fastaname):
        self._id=id
        self._fastaname=fastaname

        
        #read in the file, delete the headers, concatenate the single lines of the sequence and store them as a single string in sequence
        with open(self._fastaname) as chr_fasta:
            chr_list = chr_fasta.read().splitlines() 
        chr_list[0] = ""
        #self._sequence = "".join(chr_list)
        


        sequence_str = "".join(chr_list)
        self._sequence = list(sequence_str)
        


        #generate the reverse sequence

        converter_dictionary = {"A" : "T", "T": "A", "C": "G", "G": "C"}

        sequence_str = list(sequence_str)
        for i in range(len(sequence_str)):
            sequence_str[i] = converter_dictionary[sequence_str[i]]
        self._revsequence = sequence_str

    #add befehl   
    def __add__(self,chromosome):
        if not isinstance(chromosome,Chromosome):
            raise TypeError
        self.sequence=self.sequence+ chromosome.sequence
        self.revsequence=self.revsequence+ chromosome.revsequence
        return self
    
    
    #getter für id, sequence & revsequence

    @property
    def id(self):
        return self._id
    @id.setter
    def id(self, value):
        if not isinstance(value, int):
            raise TypeError("ID must be an Integer.")
        self._id = value
        
    @property
    def sequence(self):
        return self._sequence
    #@sequence.setter
    #def sequence(self, value):
    #    if not isinstance(value, str):
    #        raise TypeError("Sequence must be a String.")
    #    self._sequence = value
        
    @property
    def revsequence(self):
        return self._revsequence
    @revsequence.setter
    def revsequence(self, value):
        if not isinstance(value, str):
            raise TypeError("RevSequence must be a String.")
        self._revsequence = value
        
    @property
    def fastaname(self):
        return self._fastaname


class Gene:
    """
    an object for the gene sequence with the attributes:
    id -> gene id
    name -> gene name
    chr -> which chromosome the gene is on in integer form, 1-15 und 17 = mitochondrial chromosome and 18 = 2-micron plasmid
    sequence -> sequence of the dna which is transcribed
    count ->
    sequence_binding -> records if the gene is bound (=1) or currently not bound (=0), important for transcription vs replikation (i think)

    possible operations are:
    gene.mid -> provides the gene name
    gene.chr -> provides the chromosome on which the gene is on
    gene.sequence -> provides the sequence of the gene
    gene.name -> provides the name of the gene
    """
    def __init__(self, mid, name, chr, sequence, location, count=0):
        
        self.__name = name
        self.__mid = mid
        self.__location = location
        self.__chr  = chr
        self.__sequence = sequence 
        self.sequence_binding=[0]*len(sequence)
        self.rnas_transcribed=0 
        self.pol_on_gen = []


    ###### COMMENT for DATA GROUP #######
    #feel free to replace 'sequence'-information by start-, end-positions and strand (+/-)

    ###### COMMENT FOR REPLICATION_GROUP #######
    #count: 1 for unreplicated gene, 2 for copied gene 

    #@property
    #def location(self):
    #    return self.__location

    @property
    def mid(self):
        return self.__mid

    @property
    def chr(self):
        return self.__chr

    @property
    def sequence(self):
        return self.__sequence
        
    @property
    def name(self):
        return self.__name

    @property
    def location(self):
        return self.__location
    
    
    @sequence.setter
    def sequence(self, value):
        if not isinstance(value, str):
            raise Exception("sequence must be a string")
            # TODO: check for valid nucleotides here
        self.__sequence = value.upper()



def createwholegenome(chr_list):

    whole_genome = ""
    for i in range(len(chr_list)):
        whole_genome += chr_list[i].sequence

    return whole_genome

def createchromosomes():
    
    chr1=Chromosome("Chr 1","fsa_sequences/S288C_Chromosome I.fsa")
    chr2=Chromosome("Chr 2","fsa_sequences/S288C_Chromosome II.fsa")
    chr3=Chromosome("Chr 3","fsa_sequences/S288C_Chromosome III.fsa")
    chr4=Chromosome("Chr 4","fsa_sequences/S288C_Chromosome IV.fsa")
    chr5=Chromosome("Chr 5","fsa_sequences/S288C_Chromosome V.fsa")
    chr6=Chromosome("Chr 6","fsa_sequences/S288C_Chromosome VI.fsa")
    chr7=Chromosome("Chr 7","fsa_sequences/S288C_Chromosome VII.fsa")
    chr8=Chromosome("Chr 8","fsa_sequences/S288C_Chromosome VIII.fsa")
    chr9=Chromosome("Chr 9","fsa_sequences/S288C_Chromosome IX.fsa")
    chr10=Chromosome("Chr 10","fsa_sequences/S288C_Chromosome X.fsa")
    chr11=Chromosome("Chr 11","fsa_sequences/S288C_Chromosome XI.fsa")
    chr12=Chromosome("Chr 12","fsa_sequences/S288C_Chromosome XII.fsa")
    chr13=Chromosome("Chr 13","fsa_sequences/S288C_Chromosome XIII.fsa")
    chr14=Chromosome("Chr 14","fsa_sequences/S288C_Chromosome XIV.fsa")
    chr15=Chromosome("Chr 15","fsa_sequences/S288C_Chromosome XV.fsa")
    chr16=Chromosome("Chr 16","fsa_sequences/S288C_Chromosome XVI.fsa")
    chrmito=Chromosome("Chr 17","fsa_sequences/S288C_Chromosome Mito.fsa")

    chr_list=[chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chrmito]
        
    return chr_list

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
        gene[gene_id[i]] = Gene(gene_id[i], gene_name[i], which_chr[i], gene_seq[i], loc_list[i])

    return gene