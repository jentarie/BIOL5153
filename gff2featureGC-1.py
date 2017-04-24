#!/usr/bin/python


                 
gff_in = "watermelon.gff"
fsa_in = "watermelon.fsa"
gff_file = open(gff_in, 'r')
fsa_file = open(fsa_in, 'r')

 

genome = ""

line_number = 0

for line in fsa_file:
#    print(str(line_number) + ": " + line )
    line = line.rstrip('\n')
    if line_number > 0:
        genome = genome + line
    line_number += 1
len_genome = len(genome)
#print(len_genome)

#write each feature class to separate files



exon_out = open("exon.txt", "w")
CDS_out = open("CDS.txt", "w")
intron_out = open("intron.txt", "w")
misc_out = open("misc_feature.txt", "w")
rRNA_out = open("rRNA.txt", "w")
repeat_out = open("repeat_region.txt", "w")
tRNA_out = open("tRNA.txt", "w")
for line in gff_file:
    if "exon" in line:
        exon = line
        exon_out.write(exon)
    if "CDS" in line:
        CDS = line
        CDS_out.write(CDS)
    if "intron" in line:
        intron = line
        intron_out.write(intron)
    if "misc_feature" in line:
        misc_feature = line
        misc_out.write(misc_feature)
    if "rRNA" in line:
        rRNA = line
        rRNA_out.write(rRNA)
    if "repeat_region" in line:
        repeat_region = line
        repeat_out.write(repeat_region)
    if "tRNA" in line:
        tRNA = line
        tRNA_out.write(tRNA)

gff_file.close()
fsa_file.close()
exon_out.close()
CDS_out.close()
intron_out.close()
misc_out.close()
rRNA_out.close()
repeat_out.close()
tRNA_out.close()

#trim first line from fasta file and save to new file, fsa_seq.txt
fsa_seq = open("fsa_seq.txt", "w")
fsa_file = open("watermelon.fsa", 'r')
watermelon_fsa = fsa_file.readlines()[1:]
fsa_seq.write(str(watermelon_fsa))
fsa_seq.close()
fsa_file.close()

fsa_seq = open("fsa_seq.txt").read()

#find exon length and GC content
all_exons = ""

exon_locations = open("exon.txt")

for line in exon_locations:
    positions = line.split('\t')
    exon_start = int(positions[3])
    exon_stop = int(positions[4])
    exon = fsa_seq[exon_start - 1:exon_stop]
    exon_length = len(exon)
#    print("exon_start: " + str(start))
#    print("exon_stop: " + str(stop))
#    print("exon_length: " + str(exon_length))
#    print("exon: " + str(exon))
    all_exons = all_exons + exon
    exon_percent = ((len(all_exons) / len_genome ) * 100)
    exon_g_count = all_exons.count('G')
    exon_c_count = all_exons.count('C')
    exon_gc_count = exon_g_count + exon_c_count
    exon_gc = round(((exon_gc_count/len(all_exons))*100), 2)
#print("length of exon coding sequence is: " + str(len(all_exons)))
#print(exon_percent)
#print(exon_gc)

#find intron length and GC content
CDS_coding_sequence = ""

CDS_locations = open("CDS.txt")

for line in CDS_locations:
    positions = line.split('\t')
    CDS_start = int(positions[3])
    CDS_stop = int(positions[4])
    CDS = fsa_seq[CDS_start -1:CDS_stop]
    CDS_coding_sequence = CDS_coding_sequence + CDS
    CDS_length = len(CDS_coding_sequence)
    CDS_percent = ((len(CDS_coding_sequence) / len_genome ) * 100)
    CDS_g_count = CDS_coding_sequence.count('G')
    CDS_c_count = CDS_coding_sequence.count('C')
    CDS_gc_count = CDS_g_count + CDS_c_count
    CDS_gc = round(((CDS_gc_count/len(CDS_coding_sequence))*100), 2)
#print("length of CDS coding sequence is: " + str(len(CDS_coding_sequence)))

#find intron length and GC content
intron_coding_sequence = ""

intron_locations = open("intron.txt")

for line in intron_locations:
    positions = line.split('\t')
    intron_start = int(positions[3])
    intron_stop = int(positions[4])
    intron = fsa_seq[intron_start - 1:intron_stop]
    intron_coding_sequence = intron_coding_sequence + intron
    intron_length = len(intron_coding_sequence)
    intron_percent = ((len(intron_coding_sequence) / len_genome ) * 100)
    intron_g_count = intron_coding_sequence.count('G')
    intron_c_count = intron_coding_sequence.count('C')
    intron_gc_count = intron_g_count + intron_c_count
    intron_gc = round(((intron_gc_count/len(intron_coding_sequence))*100), 2)
#print("length of intron coding sequence is: " + str(len(intron_coding_sequence)))

#find misc_feature length and GC content
misc_coding_sequence = ""

misc_locations = open("misc_feature.txt")

for line in misc_locations:
    positions = line.split('\t')
    misc_start = int(positions[3])
    misc_stop = int(positions[4])
    misc = fsa_seq[misc_start -1:misc_stop]
    misc_coding_sequence = misc_coding_sequence + misc
    misc_length = len(misc_coding_sequence)
    misc_percent = ((len(misc_coding_sequence) / len_genome ) * 100)
    misc_g_count = misc_coding_sequence.count('G')
    misc_c_count = misc_coding_sequence.count('C')
    misc_gc_count = misc_g_count + misc_c_count
    misc_gc = round(((misc_gc_count/len(misc_coding_sequence))*100), 2)
#print("length of misc_feature coding sequence is: " + str(len(misc_coding_sequence)))

#find rRNA length and GC content
rRNA_coding_sequence = ""

rRNA_locations = open("rRNA.txt")

for line in rRNA_locations:
    positions = line.split('\t')
    rRNA_start = int(positions[3])
    rRNA_stop = int(positions[4])
    rRNA = fsa_seq[rRNA_start -1:rRNA_stop]
    rRNA_length = len(rRNA)
    rRNA_coding_sequence = rRNA_coding_sequence + rRNA
    rRNA_length = len(rRNA_coding_sequence)
    rRNA_percent = ((len(rRNA_coding_sequence) / len_genome ) * 100)
    rRNA_g_count = rRNA_coding_sequence.count('G')
    rRNA_c_count = rRNA_coding_sequence.count('C')
    rRNA_gc_count = rRNA_g_count + rRNA_c_count
    rRNA_gc = round(((rRNA_gc_count/len(rRNA_coding_sequence))*100), 2)
#print("length of rRNA coding sequence is: " + str(len(rRNA_coding_sequence)))

#find repeat_region length and GC content
repeat_coding_sequence = ""

repeat_locations = open("repeat_region.txt")

for line in repeat_locations:
    positions = line.split('\t')
    repeat_start = int(positions[3])
    repeat_stop = int(positions[4])
    repeat = fsa_seq[repeat_start -1:repeat_stop]
    repeat_length = len(repeat)
    repeat_coding_sequence = repeat_coding_sequence + repeat
    repeat_percent = ((len(repeat_coding_sequence) / len_genome ) * 100)
    repeat_g_count = repeat_coding_sequence.count('G')
    repeat_c_count = repeat_coding_sequence.count('C')
    repeat_gc_count = repeat_g_count + repeat_c_count
    repeat_gc = round(((repeat_gc_count/len(repeat_coding_sequence))*100), 2)
#print("length of repeat_region coding sequence is: " + str(len(repeat_coding_sequence)))

#find tRNA length and GC content
tRNA_coding_sequence = ""

tRNA_locations = open("tRNA.txt")

for line in tRNA_locations:
    positions = line.split('\t')
    tRNA_start = int(positions[3])
    tRNA_stop = int(positions[4])
    tRNA = fsa_seq[tRNA_start - 1:tRNA_stop]
    tRNA_length = len(tRNA)
    tRNA_coding_sequence = tRNA_coding_sequence + tRNA
    tRNA_length = len(tRNA_coding_sequence)
    tRNA_percent = ((len(tRNA_coding_sequence) / len_genome ) * 100)
    tRNA_g_count = tRNA_coding_sequence.count('G')
    tRNA_c_count = tRNA_coding_sequence.count('C')
    tRNA_gc_count = tRNA_g_count + tRNA_c_count
    tRNA_gc = round(((tRNA_gc_count/len(tRNA_coding_sequence))*100), 2)
#print("length of tRNA coding sequence is: " + str(len(tRNA_coding_sequence)))



print("exon",'\t\t', CDS_length, '\t', "(" + str(round(CDS_percent, 1)) + "%)",'\t', CDS_gc)
print("intron", '\t\t', intron_length, '\t', "(" + str(round(intron_percent, 1)) + "%)",'\t', intron_gc)
print("misc_feature", '\t', misc_length, '\t', "(" + str(round(misc_percent, 1)) + "%)",'\t', misc_gc)
print("repeat_region", '\t', repeat_length, '\t', "(" + str(round(repeat_percent, 1)) + "%)",'\t', repeat_gc)
print("rRNA", '\t\t', rRNA_length, '\t', "(" + str(round(rRNA_percent, 1)) + "%)",'\t', rRNA_gc)
print("tRNA", '\t\t', tRNA_length, '\t', "(" + str(round(tRNA_percent, 1)) + "%)",'\t', tRNA_gc)