import random as random
import numpy as np



sequence_to_replicate = ['A','G','C','T','T','G','A','C','T','A',
                         'A','G','C','T','T','G','A','C','T','A',
                         'A','G','C','T','T','G','A','C','T','A',
                         'A','G','C','T','T','G','A','C','T','A',
                         'A','G','C','T','T','G','A','C','T','A',
                         'A','G','C','T','T','G','A','C','T','A',
                         'A','G','C','T','T','G','A','C','T','A',
                         'A','G','C','T','T','G','A','C','T','A',
                         'A','G','C','T','T','G','A','C','T','A',
                         'A','G','C','T','T','G','A','C','T','A',
                         'A','G','C','T','T','G','A','C','T','A']

prob= []
new_chrom= []
mut = []
ts_position = []
tv_position = []
ts_nucleotids =[]
tv_nucleotids = []


def rep_mutation():
    for i in range(len(sequence_to_replicate)):
            x = random.randint(0,100)
            #dictionary
            nucleotids = {'A':'T','T':'A','G':'C','C':'G'}
            transition = {'A':'C','T':'G','G':'T','C':'A'}
            transversion = {'A':'A','T':'T','G':'G','C':'C'}


            if x == 0:                             #transversion p=0.01
                new_chrom.append(transversion[sequence_to_replicate[i]])     
                prob.append(x)
            elif x == 1:                            #transition p=0.01
                new_chrom.append(transition[sequence_to_replicate[i]])
                prob.append(x)
            else:                                    #complementary strand p=0.98
                new_chrom.append(nucleotids[sequence_to_replicate[i]])   
                prob.append(x)

                
def proof_reading():
    
    for i in prob:
        if i == 1:
            mut.append('ts')
        elif i == 0:
            mut.append('tv')
        else:
            mut.append(0)
    
    
    if 'ts' or 'tv' in mut:
        position_ts = [i for i,j in enumerate(mut) if 'ts' == j]
        ts_position.append(position_ts)
        position_tv = [i for i,j in enumerate(mut) if 'tv' ==j]
        tv_position.append(position_tv)
    return mut

    
def show_mutation():
    
    if 'ts' or 'tv' in mut:
    	#3 LISTEN
    	#[3,4]
    	#[A, C, G, A, C]
    	#[T, G, c, T, G]

        i = 0
        j = 0   
        while i in range(len(ts_position)):
            for x in ts_position[i]:
                a = sequence_to_replicate[x]
                b = new_chrom[x]
                c = [a,b]
                ts_nucleotids.append(c)
                i +=1
        while j in range(len(tv_position)):
            for y in tv_position[j]:
                d = sequence_to_replicate[y]
                e = new_chrom[y]
                f = [d,e]
                tv_nucleotids.append(f)
                j +=1

    else:
        pass




    
    
rep_mutation()
proof_reading()

    

print("original")
print(sequence_to_replicate)
print("probabilities")
print(prob)    
print('replica')
print(new_chrom)
print("mutation")
print(mut)
print("number of mutations :", len(ts_nucleotids)+len(tv_nucleotids))
print("transition at:", ts_position,"total:", len(ts_nucleotids))
print("nucleotids ori > repl:",ts_nucleotids)
print("transversion at:", tv_position, "total:", len(tv_nucleotids))
print("nucleotids ori > repl:",tv_nucleotids)
