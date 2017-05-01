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


#reverse complement dna
def rev_compl(fragment):
    compl = fragment.replace("ACTG", "tgac")
    return compl.upper()

#function for G+C
### Try to extract all the gene sequences from the gcc file for exons only.
#fine if on the coding strand, but could be on the negative strang (function for reverse complement)
#print/store/build a fasta sequence of each gene (coding sequence) eg:
#>cox1
#TGATTT...


#key = feeature_type, value
feature_sequences = {}

#key = gene, value
gene_sequences = {}

gff_file = 'watermelon.gff'
fsa_file = 'watermelon.fsa'
gff_in = open(gff_file, 'r')
fsa_in = open(fsa_file, 'r')


#extract names of genes
#def gene_extract(file):
#    for line in file:
#        if "CDS" in line:
#            gene = line.split("Gene ")[1].split(";")[0]
#            return gene
#print(gene)
#print(gene_extract(gff_in))

gene_list = []
for line in gff_in:
    if "CDS" in line:
        gene = line.split("Gene ")[1].split(" ;")[0]
        if gene[0] not in gene_list:
            gene_list.append(gene)
print(gene_list)

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

for line in gff_in:
#    line = line.rstrip('\n')
#    types = line.split('; type ')
#    other_type = types[len(types)-1]
#    print(other_type)
    gene = line.split("Gene ")[1].split(" ;")[0]
    fields = line.split('\t')
    type = fields[2]
    start = int(fields[3])
    end = int(fields[4])
    #print(type, "\t", start, "\t", end)
    fragment = genome[start - 1:end]
    fragment = clean_seq(fragment)
    

  
    
    if gene in gene_sequences:
        gene_sequences[gene] += fragment
    else:
        gene_sequences[gene] = fragment

  
for gene, sequence in gene_sequences.items():
    print(gene + '\t' + sequence)

#    g_count = fragment.count('G')
#    c_count = fragment.count('C')
#    gc_count = g_count + c_count   

    #populate our dictionary
#    if type in feature_sequences:
#        feature_sequences[type] += fragment
#    else:
#        feature_sequences[type] = fragment
    
#gff_in.close()

#for feature, sequence in feature_sequences.items():
#    print(feature + '\t' + str(len(sequence)))
#    #print(sequence)
#    #call function here to calculate GC content



types = ['cds', 'intron', 'tRNA', 'rRNA', 'repeats', 'misc']       
#print(len(cds))
#sys.exit()