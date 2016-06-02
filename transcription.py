import processes 
import molecules
import random
import numpy as np
from matplotlib import pyplot as plt


################### TO DO !!!!!! ########################
#look to data: elongation-processes per 1s and transcripion-rates (Huyen)

#include real data and global data-construction

#get information how many mRNAs are transcribed at the moment

#at the moment: only one polymerase transcribes an RNA ! (Paula)

#include rRNA and tRNA - transcription

#######################################################

#mRNA half lifes

########################################################

class Transcription(processes.Process):
	"Transcription process needs to be described."


	def __init__(self, id, name):
		""" All new Transcription variables need to be implemented."""
		
		super().__init__(id, name)
		#initialized through set_states in model: 
		#id_enzymes= one single RNA Polymerase -> not necessary as polymerase is an argument of Transcription 
		# -> multiple transcription processes with different polymerases possible
		#id_substrates= a dictionary with all genes of the genome: {geneid1: Geneobject, geneid2: GeneObject, ...} 

		#self.mypolymerase=polymerase
		#self.my_nucleotides=nucleotide_pool

				#polymerase=RNAPolymeraseII(needs to be specified) -> initialization can be done in transcription or model
				#expect polymerase to be unbound -> check in model.py

		#length of occupied sequence by RNA-Polymerase II on DNA: already validated with data!
		self.polymerase_size=17

		#####for visualization of selected genes
		#self.allgenes=[[],[]]
		#########################################


	def __repr__(self):
		# todo: each process class should have something like this
		pass

	def __str__(self):
		# todo: each process class should have something like this
		pass



	def update(self, model):
	#def update(self, genedic, rna_pool):

		genedic=model.genes
		rna_pool=model.mrnas
		self.mypolymerase=self.enzyme_ids
		self.my_nucleotides=model.nucleotides

		#number still needed: empirically! 
		update_per_s=2000

		gene_id, weights=self.make_weights(genedic)

		for steps in range(update_per_s):
			#rna = self.onestep(genedic)
			rna=self.onestep(genedic, weights, gene_id)
			if isinstance(rna, molecules.MRNA):
				#rna_pool.append(rna)
				if rna.name in model.states:
					model.states[rna.name].append(rna)
				else:
					model.states[rna.name] = [rna]

    	####visualization of selected genes ######
		#plt.plot(range(len(self.allgenes[0])),self.allgenes[1])
		#plt.xlabel(self.allgenes[0])
		#plt.savefig('tests/count_histogram_genes.pdf')
		###########################################

		#return rna_pool


			
	#def onestep(self,genedic):			#expect a dictionary of genes 
	def onestep(self, genedic, weights, gene_ids):

		transc_gene=self.select_gene(genedic)

		#testprints
		#print(transc_gene.name)
		#print(transc_gene.sequence)	
		#print(transc_gene.pol_on_gene)

		######for visualization of selected genes ################
		#if transc_gene.name in self.allgenes[0]:
		#	self.allgenes[1][self.allgenes[0].index(transc_gene.name)]+=1
		#else:
		#	self.allgenes[0].append(transc_gene.name)
		#	self.allgenes[1].append(1)
		#########################################################

		if not transc_gene.pol_on_gene:
			if transc_gene.sequence_binding[0]==0 and self.mypolymerase.count>0:
				self.initiate(transc_gene)
		else:
			new_pol = self.rand_distr(transc_gene, gene_ids, weights) 
			#function that will compare the probability of ocurrence of the transcription rate of the gene with a random number

			#if transc_gene.sequence_binding[self.polymerase_size-1]==0 and self.mypolymerase.count>0 and random.randint(1,200)==1:
			if transc_gene.sequence_binding[self.polymerase_size-1]==0 and self.mypolymerase.count>0 and new_pol == 1:

				#check whether the whole region the polymerase will occupy is free
				self.initiate(transc_gene)
			else:
				pol_position=random.choice(transc_gene.pol_on_gene)
				mrna = self.transcribe(transc_gene, pol_position)
				if isinstance(mrna, molecules.MRNA):
					#if mrna.name in model.states:
					#	model.states[mrna.name].append(mrna)
					#else:
					#	model.states[mrna.name] = [mrna]
					return mrna

	def make_weights(self, genedic):
		transc_rate=[]
		gene_ids=list(genedic.keys())
		for g in gene_ids:
			transc_rate.append(genedic[g].transrate)
		transc_rate=np.array(transc_rate)

		weights=transc_rate
		weights=weights/sum(weights)
		return (gene_ids, weights)



	def select_gene(self, genedic):

		#transcribed gene identified: no random selection, but weighted by copies of each gene

		copies=[]
		transc_rate=[]
		weight_for_binding=[]
		gene_ids=list(genedic.keys())
		for g in gene_ids:
			copies.append(genedic[g].count)
			transc_rate.append(genedic[g].transrate)
			if genedic[g].pol_on_gene:
				weight_for_binding.append(100000000)
			else:
				weight_for_binding.append(1)
		copies=np.array(copies)
		transc_rate=np.array(transc_rate)

		weights=copies*transc_rate*weight_for_binding
		weights=weights/sum(weights)

		rand_index=np.random.choice(range(len(weights)),p=weights)

			
		transc_gene=genedic[gene_ids[rand_index]]

		return transc_gene


	def initiate(self,gene):

		""" for one Polymerase: look for a gene to transcribe. when ORF is found: initiate """

		#### from the data group: we expect an self.genes-dictionary in model.py containing all genes

		#print(gene.sequence_binding)
		if isinstance(self.mypolymerase, molecules.RNAPolymeraseII): 
			gene.pol_on_gene.append(0)
			mrna=molecules.MRNA("mRNA_{}".format(gene.mid), "mRNA_{0}".format(gene.name.split("_")[-1]), '',)

		#if isinstance(self.mypolymerase, molecules.RNAPolymeraseI):
			#gene.pol_on_gene[0].append(0)
			#rna=molecules.RNA("RNA_{}".format(gene.mid), "mRNA_{0}".format(gene.name.split("_")[-1]), '',)

		#if isinstance(self.mypolymerase, molecules.RNAPolymeraseIII):
			#gene.pol_on_gene[0].append(0)
			#rna=molecules.RNA("RNA_{}".format(gene.mid), "mRNA_{0}".format(gene.name.split("_")[-1]), '',)

		
		#bigger polymerase: has size of polymerase_size in both directions
		if self.polymerase_size>len(gene.sequence_binding):
			for i in range(len(gene.sequence_binding)):
				gene.sequence_binding[i]=mrna
		else:
			for i in range(self.polymerase_size):
				gene.sequence_binding[i]=mrna		
			
		self.mypolymerase.count+=-1


	def transcribe(self, gene, position):

		""" elongate mRNA for given ORF for only one step. if ORF ends: call terminate-function """
		pos=position
		index=gene.pol_on_gene.index(pos)

		mrna=gene.sequence_binding[pos]

		#appending the correct codon from the coding strand
		nuc=gene.sequence[pos]
		if nuc=='T':
			mrna.sequence+='U'
			self.my_nucleotides.count_nuc['U']+=-1
		else:
			mrna.sequence+=nuc
			self.my_nucleotides.count_nuc[nuc]+=-1

		#if we are not on the end of the ORF-string
		if pos+1<len(gene.sequence_binding):
			#!!!gene.pol_num +=1

			#if coding region is ending: no macromolecule is thought to disturb the elongation process
			if pos+self.polymerase_size>=len(gene.sequence_binding):
				gene.pol_on_gene[index]+=1
				if pos >= self.polymerase_size:
					gene.sequence_binding[pos-self.polymerase_size]=0
				return 0

			#if position space in front of polymerase is empty
			elif gene.sequence_binding[pos+self.polymerase_size]==0:
				gene.pol_on_gene[index]+=1
				gene.sequence_binding[pos+self.polymerase_size]=mrna
				if pos >= self.polymerase_size:
					gene.sequence_binding[pos-self.polymerase_size]=0
				return 0

		else:
			return self.terminate(gene, position)
		


	def terminate(self, gene, position):
		""" separate mRNA-DNA-Polymerase-complex. release and store new mRNA """
		
		self.mypolymerase.count+=1
		rna = gene.sequence_binding[position]
		for p in range(self.polymerase_size+1):
			gene.sequence_binding[position-p]=0
		del gene.pol_on_gene[gene.pol_on_gene.index(position)]

		return rna


	def rand_distr(self, gene, gene_ids, weights): #if gene_ids and weight are global we don't need to send them in the function 
		ran=random.uniform(0,5)
		index_id = gene_ids.index(gene.mid)
		if weights[index_id]>= ran:
			return 1
		else:
			return 0

