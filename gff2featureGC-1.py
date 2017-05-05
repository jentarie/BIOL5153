#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 12:56:29 2017

@author: jen
"""

#!/usr/bin/python
 
 
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
 
def GC_frequency(sequence, sig_figs=2):
    length = len(sequence)
    GC_count = sequence.count('G') + sequence.count('C')
    GC_freq = GC_count/length
    return(round(GC_freq, sig_figs))


 
#reverse complement dna
def rev_compl(fragment):
    compl = fragment.replace("ACTG", "tgac")
    return compl.upper()
 

#key = feeature_type, value
feature_sequences = {}

 
#tuple method - almost works!
gene_sequences = {}
 

 
 
gff_file = '/home/jen/Desktop/BIOL5153/watermelon_files/watermelon.gff'
fsa_file = '/home/jen/Desktop/BIOL5153/watermelon_files/watermelon.fsa'
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
 
#gene_list = []
#for line in gff_in:
#    if "CDS" in line:
#        gene = line.split("Gene ")[1].split(" ;")[0]
#        if gene[0] not in gene_list:
#            gene_list.append(gene)
#print(gene_list)
 
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
    #gene = line.split("Gene ")[1].split(" ;")[0]
    fields = line.split('\t')
    type = fields[2]
    start = int(fields[3])
    end = int(fields[4])
    #print(type, "\t", start, "\t", end)
   
    #extract and clean sequence
    fragment = genome[start - 1:end]
    fragment = clean_seq(fragment)   
 
    if (fields[6] == '-'):
        fragment = rev_compl(fragment)
 
#    if gene in gene_sequences:
#        gene_sequences[gene] += fragment
#    else:
#        gene_sequences[gene] = fragment
       
    
    
    #get the gene name
    if(type == 'CDS'):
        attributes = fields[8].split(' ; ')
        #print(attributes[0])
       
        gene_fields = attributes[0].split(' ')
        gene_name = gene_fields[1]   
    #determine if there are multiple nmber of exons
    #get the exon number
        if('exon' in gene_fields):
            exon_num = gene_fields[-1]
            #exon_num = gene_fields[3]
            #print(gene_name, exon_num)
        else:
            exon_num = '0'
 

#build dictionary using tuple as key
        gene_sequences[(gene_name, exon_num)] = fragment
#print(gene_sequences)   
 
gff_in.close()        


genes_compiled={}

for gene_name, exon_num in sorted(gene_sequences):
    #print(">" + gene_name + ' ' + str(exon_num) + '\n' + gene_sequences[(gene_name, exon_num)])
    if gene_name in genes_compiled:
        genes_compiled[gene_name] += gene_sequences[(gene_name, exon_num)]
    else:
            genes_compiled[gene_name] = gene_sequences[(gene_name, exon_num)]
#print(genes_compiled)


for gene, sequence in genes_compiled.items():
    print(">" + gene + '\n' + sequence)

        
for gene, sequence in genes_compiled.items():
    GC_content = GC_frequency(sequence)
    print("GC content of " + gene + ':\t ' + str(GC_content) + '%')

  

 
#for feature, sequence in feature_sequences.items():
#    print(feature + '\t' + str(len(sequence)))
#    #print(sequence)
#    #call function here to calculate GC content
 
 
 
#types = ['cds', 'intron', 'tRNA', 'rRNA', 'repeats', 'misc']
#print(len(cds))
#sys.exit()