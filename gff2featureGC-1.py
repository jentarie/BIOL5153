#!/usr/bin/python

import sys

def clean_seq(input_seq):
    clean = input_seq.upper()
    clean = clean.replace('N', '')
    return clean

def nuc_freq(sequence, base, sig_figs=2):
    #calculate the length of the sequence
    length = len(sequence)


    # count the number of this nucleotide
    base_count = sequence.count(base)

    # calculate base frequency
    base_freq = base_count/length

    # return the frequency and length of nucleotide
    return(length, round(base_freq, sig_figs))

gff_file = 'watermelon.gff'
fsa_file = 'watermelon.fsa'
gff_in = open(gff_file, 'r')
fsa_in = open(fsa_file, 'r')

genome = ""

line_number = 0

for line in fsa_in:
#    print(str(line_number) + ": " + line )
    line = line.rstrip('\n')
    if line_number > 0:
        genome = genome + line
    line_number += 1
#print(len(genome))

fsa_in.close()

cds = ""
intron = ""
tRNA = ""
rRNA = ""
repeats = ""
misc = ""

for line in gff_in:
#    line = line.rstrip('\n')
#    types = line.split('; type ')
#    other_type = types[len(types)-1]
#    print(other_type)

    fields = line.split('\t')
    type = fields[2]
    start = int(fields[3])
    end = int(fields[4])
    #print(type, "\t", start, "\t", end)
    fragment = genome[start - 1:end]

#    g_count = fragment.count('G')
#    c_count = fragment.count('C')
#    gc_count = g_count + c_count

    fragment = clean_seq(fragment)


    if type == 'CDS':
        cds += fragment
    if type == 'intron':
        intron += fragment
    if type == 'tRNA':
        tRNA += fragment
    if type == 'rRNA':
        rRNA += fragment
    if type == 'repeat_region':
        repeats += fragment
    if type == 'misc_feature':
        misc += fragment
#print(cds.count('G'))
types = ['cds', 'intron', 'tRNA', 'rRNA', 'repeats', 'misc']
#print(len(cds))
#sys.exit()



# loop over features:
i = 0
#for feature_type in [CDS, intron, tRNA, rRNA, repeat_region, misc_feature]:
for feature_type in [cds, intron, tRNA, rRNA, repeats, misc]:
    #loop over 4 nucleotides
    for nucleotide in ['A', 'C', 'G', 'T']:
        #calculate the nucleotide composition for each feature
        (feature_length, feature_comp) = nuc_freq(sequence=feature_type, base=nucleotide, sig_figs=2)
        print(types[i] + "\t\t" + str(feature_length) + ' nucleotides\t' + str(feature_comp) + "% " + nucleotide)
    i = i + 1
