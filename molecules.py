__author__ = 'max'


class BioMolecule:
    """
    A generic molecule that has basic attributes like name and
    mass.

    @type name: str
    @type mass: float
    """

    def __init__(self, mid, name, mass=0):
        self.mid = mid                      #id ?
        self.name = name
        self.mass = mass

    @property
    def mid(self):
        return self.__mid

    @mid.setter
    def mid(self, value):
        self.__mid = value

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, value):
        self.__name = value

    @property
    def mass(self):
        return self.__mass

    @mass.setter
    def mass(self, value):
        if not (isinstance(value, float) or isinstance(value, int)):
            raise Exception("mass must be numeric")
        else:
            self.__mass = value

    def __repr__(self): #string "self.name,type"		#print(list(object))
        return ','.join([self.name, str(type(self))])

   # def __str__(self):	#print(object)
        # todo: each class should have something like this
 #       pass


class Polymer(BioMolecule):
    """
    A polymer molecule that has a sequence attribute which is
    accessible via indexing the object.

    @type name: str
    @type sequence: str
    @type mass: float
    """

    def __init__(self, mid, name, sequence, mass=0):
        super().__init__(mid, name, mass)
        self.__sequence = sequence

    def __getitem__(self, value):
        return self.sequence[value]

    def __setitem__(self, key, value):
        self.sequence[key] = value

    @property
    def sequence(self):
        return self.__sequence

    @sequence.setter
    def sequence(self, value):
        if not isinstance(value, str):
            raise Exception("sequence must be a string")
            # TODO: check for valid nucleotides here
        self.__sequence = value.upper()


class BioMoleculeCount(BioMolecule):        #new variable: number of molecules of this type
    def __init__(self, mid, name, count=0):
        super().__init__(mid, name)
        self.count = count

    @property
    def count(self):
        return self.__count

    @count.setter
    def count(self, value):
        self.__count = value


class NucleotidPool(BioMolecule):
    def __init__(self, count):
        super().__init__(mid, name, count)

        self.count_nuc={'A': count, 'C': count, 'G': count, 'T': count, 'U': count}


class MRNA(Polymer):
    def __init__(self, mid, name, sequence, mass=0):
        super().__init__(mid, name, sequence, mass)
        self.sequence_triplet_binding = [0] * (len(sequence) // 3)

    def calculate_mass(self):
        self.mass = 0
        NA_mass = {'A': 1.0, 'U': 2.2, 'G': 2.1, 'C': 1.3}
        for na in self.sequence:
            self.mass += NA_mass[na]


class Protein(Polymer):
    """
    Protein with Polymer features and mass calculation. A global class
    attribute counts the number of proteins that have been instantiated.

    A protein can be elongated using the "+" operator:

    >> protein = Protein(1, "prot", "MVFT")
    >> protein + "A"
    >> protein.sequence
    MVFTA



    """
    number_of_proteins = 0

    def __init__(self, mid, name, sequence, mass=0):
        super().__init__(mid, name, sequence, mass)

    def __iadd__(self, AS):
        self.sequence = self.sequence + AS
        return self

    def calculate_mass(self):
        AA_mass = dict(A=89.0929, R=175.208, N=132.118, D=132.094, C=121.158, Q=146.144, E=146.121, G=75.0664,
                       H=155.154, I=131.172, L=131.172, K=147.195, M=149.211, F=165.189, P=115.13, S=105.092, T=119.119,
                       W=204.225, Y=181.188, V=117.146)
        for aa in self.sequence:
            self.mass += AA_mass[aa]


class Ribosome(BioMoleculeCount):
    """
    A ribosome can bind MRNA and translate it. After translation is
    finished it produces a protein.

    During initiation the ribosome checks if a given MRNA is bound
    by another ribosome and binds only if position 0 is empty.

    Elongation checks if the next codon is unbound and only elongates
    if the ribosome can move on. If the ribosome encounters a stop
    codon ("*") translation terminates. The MRNA is detached from the
    ribosome and the finished protein is returned.
    """

    def __init__(self, mid, name, count=0):
        super().__init__(mid, name, count)


class Polymerase(BioMoleculeCount):
    """
    A polymerase that can elongate nucleotide molecules. This could be used to derive special
    RNA and DNA polymerases.
    """

    def __init__(self, mid, name, count=0):
    	super().__init__(mid, name, count)
   


class RNAPolymeraseI(Polymerase):
    """
    A polymerase that generates rRNAs from DNA sequences.
    """
    def __init__(self, mid, name, count=0):
        super().__init__(mid, name, count)

class RNAPolymeraseII(Polymerase):
    """
    A polymerase that generates mRNAs from DNA sequences.
    """
    def __init__(self, mid, name, count=0):
        super().__init__(mid, name, count)

class RNAPolymeraseIII(Polymerase):
    """
    A polymerase that generates tRNAs from DNA sequences.
    """
    def __init__(self, mid, name, count=0):
        super().__init__(mid, name, count)





class Gene(BioMoleculeCount):

    """
    an object for the gene sequence with the attributes:
    id -> gene id
    name -> gene name
    chr -> which chromosome the gene is on in integer form, 1-15 und 17 = mitochondrial chromosome and 18 = 2-micron plasmid
    sequence -> sequence of the dna which is transcribed
    count -> number of copies of the gene
    sequence_binding -> records if the gene is bound (=1) or currently not bound (=0), important for transcription vs replikation 
    pol_on_gene -> positions of polymerases on the gene

    possible operations are:
    gene.mid -> gives the gene name
    gene.chr -> gives the chromosome on which the gene is on
    gene.sequence -> gives the sequence of the gene
    """

    def __init__(self, mid, name, chr, sequence, location, count=0):
        
        super().__init__(mid, name, count)
        self.__location = location
        self.__chr  = chr
        self.__sequence = sequence 
        self.sequence_binding=[0]*len(sequence)
        self.rnas_transcribed=0 							#number of transcribed RNAs during one transcription-process
        self.pol_on_gen = []								#nucleotides transcribed by polymerases on the gene
        self.rate=32
						
        
    ###### COMMENT for DATA GROUP #######
    #feel free to replace 'sequence'-information by start-, end-positions and strand (+/-)

    ###### COMMENT FOR REPLICATION_GROUP #######
    #count: 1 for unreplicated gene, 2 for copied gene 

    @property
    def chr(self):
        return self.__chr

    @property
    def sequence(self):
        return self.__sequence

    @property
    def location(self):
        return self.__location


        
    @sequence.setter
    def sequence(self, value):
        if not isinstance(value, str):
            raise Exception("sequence must be a string")
            # TODO: check for valid nucleotides here
        self.__sequence = value.upper()

    


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
        self.binding_molecules=[[],[]]	#list with tuples of start and end positions of occupied regions in [0] and the binding molecule in [1]

        
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


#####this method will check for bound regions in the chromosome######

    def chromosome_bound(self, range):

        if isinstance(range, list):
            #range[0]: startposition of test, range[1]: endposition of test
            iter=0
            bound_stuff=[]
            bound_found=False
            for tuples in self.binding_molecules[0]:

                if bound_found==False:
                    if tuples[0]>range[1] and iter==0:
                        return None
                    elif tuples[0]>range[0] and tuples[0]>=range[1]:
                        return None
                    elif tuples[1]<range[0] and tuples[1]<range[0]:
                        pass
                    elif tuples[0]<=range[0] and tuples[1]>=range[1]:
                        bound_stuff.append(self.binding_molecules[1][self.binding_molecules[0].index(tuples)])
                        return bound_stuff
                    elif tuples[1]<range[1]:
                        bound_stuff.append(self.binding_molecules[1][self.binding_molecules[0].index(tuples)])
                        bound_found=True
                    
                else:
                    if tuples[0]>range[1]:
                        return bound_stuff
                    elif tuples[1]>range[1]:
                        bound_stuff.append(self.binding_molecules[1][self.binding_molecules[0].index(tuples)])
                        return bound_stuff
                    else:
                        bound_stuff.append(self.binding_molecules[1][self.binding_molecules[0].index(tuples)])
                iter+=1
            return bound_stuff



        elif isinstance(range, int):
            for tuples in self.binding_molecules[0]:
                if tuples[0]<=range and tuples[1]>=range:
                    return self.binding_molecules[1][self.binding_molecules[0].index(tuples)]
                elif tuples[0]>range:
                    break
           
        else:
            print("Argument type not expected (list or int)")

    #### method needed which stores a tuple of start&end and bound molecule in the ordered (!!!) list bindnig_molecules
    def bind_to_chrom(start, end):
        pass
