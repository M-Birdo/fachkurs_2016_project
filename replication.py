from molecules import BioMoleculeCount
from molecules import BioMolecule, NucleotidPool
from processes import Process
import random as random

# generell:
    # start replication? > abhängig von Proteinzahl
    # replication fork
    # origins of Replication? ARS , Origin recognition complex bindet an A-rich DNA
    # Wo bindet Helikase, wo Polymerase
    # koordination mit Transkription
    # output
   

class Chromosome:
    def __init__(self, name, sequence):
        self.sequence = sequence
        self.name = name
        self.replication_ori_bound = False # beginnt ungebunden, 

class Helicase(BioMoleculeCount):
    def __init__(self, mid, name, count=0):
        super().__init__(mid, name, count)

    def __str__(self):
        return 'helicases: {0}'.format(self.count) 
    
class Polymerase(BioMoleculeCount):
    def __init__(self, mid, name, count=0):
        super().__init__(mid, name, count)

    def __str__(self):
        return 'polymerase: {0}'.format(self.count, self.name)

class Replication(Process):
    def __init__(self, id, name):
        super().__init__(id, name)
        self.chromosomes = {} # dict für das neue Chromosom
        self.duplication = {} # dict für duplizierte Chromosomen (Ausschlusskriterium)

    def set_states(self, substrate_ids, enzyme_ids):
        super().set_states(substrate_ids, enzyme_ids)
        if self.duplication == {}:
            self.duplication = {chrom: None for chrom in self.substrate_ids} # nur für den Anfang, weil dict noch leer

        self.chromosome_names = self.substrate_ids 
    
    def update(self, model):
    
        self.polymerase =  model.states['Polymerase3']
        self.helicase = model.states['DnaB']
        self.old_chromosomes = [model.states[chromsome_name] for chromsome_name in self.chromosome_names]
        
        for i, old_chromosome in enumerate(self.old_chromosomes[0:16]): # wenn das aufgerufene Chromosom folgende Bedingungen erfüllt wird es zur entprechenden Phase weitergeleitet
            if not old_chromosome.replication_ori_bound and not self.duplication[self.chromosome_names[i]] and self.polymerase.count > 0 and self.helicase.count > 0 and sum(pool.count_nuc.values()) > len(old_chromosome.sequence)*2:
            #and not transcription an aktueller Stelle:
                self.initiate(old_chromosome) 
            elif old_chromosome.replication_ori_bound: #and not transcription an der Stelle:
                self.elongate(old_chromosome)
            elif (self.polymerase.count == 0 or self.helicase.count == 0) and not old_chromosome.replication_ori_bound and not self.duplication[self.chromosome_names[i]]:
                continue
            else:
                raise NotImplementedError

        print(sum(pool.count_nuc.values()))
        print(self.helicase.count)
        print(self.polymerase.count)
        print(self.chromosomes["Chr 1"].sequence)
        print(len(self.chromosomes["Chr 1"].sequence))


        #for c_name in self.chromosomes:
            #print(self.chromosomes[c_name].sequence)
            
    def initiate(self, chrom: Chromosome):
        print('initiation')
                    # legt das dict 'Chromosomes' an, Name bleibt erhalten, Sequenz ändert sich   
        self.chromosomes[chrom.name]=Chromosome(chrom.name,[])
        self.polymerase.count -= 1 
        self.helicase.count -= 1
        chrom.replication_ori_bound = True 
        self.duplication[chrom.name] = False 
  
    def elongate(self, old_chrom: Chromosome):
        print('elongation')
        
        """
        generiert jeweils eine 100-Nukleotid-Liste (für 100bp/s) aus der alten Sequenz . 
        Jedem Nukleotid wird eine Zahl zw. 0-100 zugeordnet (= p(Mutation)) und ensprechend der 
        Wahrscheinlichkeit wird das neue Nukleotid in der replizierten Sequenz angehängt 

        Sequenz_to_replicate > prob > new_chrom

        """

        new_chrom = self.chromosomes[old_chrom.name]
        sequence_to_replicate = []
        prob=[]
        
        # dict für Mutation
        nucleotids = {'A':'T','T':'A','G':'C','C':'G'}
        transition = {'A':'C','T':'G','G':'T','C':'A'}
        transversion = {'A':'A','T':'T','G':'G','C':'C'}

        #100er Schritte
        if len(new_chrom.sequence) == 0:
            sequence_to_replicate = old_chrom.sequence[0:100]
        elif len(old_chrom.sequence) - len(new_chrom.sequence) >= 100:
            sequence_to_replicate = old_chrom.sequence[len(new_chrom.sequence):len(new_chrom.sequence)+100]
        elif len(old_chrom.sequence) - len(new_chrom.sequence) == 0:
            return self.terminate(new_chrom.sequence)
        else:
            sequence_to_replicate = old_chrom.sequence[len(new_chrom.sequence):]
        
        # Mutation
        for i in range(len(sequence_to_replicate)):
            x = random.randint(0,10000)
            if x == 0:                             #transversion p=0.0001
                new_chrom.sequence.append(transversion[sequence_to_replicate[i]])     
                prob.append(x)
                pool.count_nuc[sequence_to_replicate[i]] -=1
                pool.count_nuc[transversion[sequence_to_replicate[i]]] -=1  
            elif x == 1:                            #transition p=0.0001
                new_chrom.sequence.append(transition[sequence_to_replicate[i]])
                prob.append(x)
                pool.count_nuc[sequence_to_replicate[i]] -=1
                pool.count_nuc[transition[sequence_to_replicate[i]]] -=1
            else:                                   #complementary strand p=0.9998
                new_chrom.sequence.append(nucleotids[sequence_to_replicate[i]])   
                prob.append(x)
                pool.count_nuc[sequence_to_replicate[i]] -=1
                pool.count_nuc[nucleotids[sequence_to_replicate[i]]] -=1
                
        print(pool.count_nuc)
        #print("original")
        #print (sequence_to_replicate)#
        #print("probalilities")
        #print (prob)
        #print ("new_sequence")
        #print (new_chrom.sequence)
        #print ("mutation")
        #print (mut)
   
    def terminate(self, chrom, i):  
        
        self.polymerase += 1
        self.helicase += 1
        self.duplication[chrom.name] = True
        chrom.replication_ori_bound = False

        return new_chrom
        
pool = NucleotidPool(1,'pool',1000000000000)


if __name__ == '__main__':
    rep = Replication('replication', 'replication')
    #print(rep.helicase)
    #print (rep.polymerase)
    #chrom1 = Chromosome('chrom1', ['A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A''A','G','C','T','T','G','A','C','T','A','A','G','C','T','T','G','A','C','T','A'])
        
    #chrom1 = chr_list[0]
    #chrom2 = chr_list[1]
    #print (len(chrom1.sequence))

    #rep.initiate(chrom1)
    #rep.elongate(chrom1)
   

