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

    def __str__(self):	#print(object)
        # todo: each class should have something like this
        pass


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
   


class RNAPolymeraseII(Polymerase):
    """
    A polymerase that generates mRNAs from DNA sequences.
    """
    def __init__(self, mid, name, count=0):
    	super().__init__(mid, name, count)



class Gene(BioMoleculeCount):
    def __init__(self, mid, name, strand, count=0):
        super().__init__(mid, name, count)
        self.__location = None 
        self.__chr  = None
        self.__sequence = None 
        self.sequence_binding=[0]*(end-start)


    ###### COMMENT for DATA GROUP #######
    #feel free to replace 'sequence'-information by start-, end-positions and strand (+/-)

    ###### COMMENT FOR REPLICATION_GROUP #######
    #count: 1 for unreplicated gene, 2 for copied gene 

    @property
    def location(self):
        return self.__location

    @property
    def chr(self):
        return self.__chr

    @property
    def sequence(self):
        return self.__sequence
        
    @sequence.setter
    def sequence(self, value):
        if not isinstance(value, str):
            raise Exception("sequence must be a string")
            # TODO: check for valid nucleotides here
        self.__sequence = value.upper()
