# -*- coding: utf-8 -*-
"""
Created on Thu Dec 04 14:32:21 2014

@author: marbec
"""
from datetime import datetime
import random
from operator import itemgetter
from math import fabs
import scipy.stats as ss

print(datetime.now())

#The GWAS data should be in the format chromosomenumber:baseposition [tab] p-value (3:476672   0.034)
#example  specify_GWAS_file = "test_GWAS.tabular"

specify_bed_file = "HumanForebrainVISTA+500bp_noOverlappingenhancers.bed"
specify_GWAS_file = "test_GWAS.tabular"
name_of_resultfile_for_tested_elements = "SNP_associations_in_bed_file_elements.txt"
name_of_resultfile_for_permutation = "Permutationresults.txt"
zscore_of_tested_elements = "TotalZscore_of_tested_elements.txt"
zscore_of_permutations = "TotalZscores_of_permutations.txt"



#Open bed file and read in regions 
bed_file = open(specify_bed_file, "r")

regions = []
region_index = 0
for line in bed_file:
    this_element = line.split("\t")
    output = []
    for value in this_element:
        try:
            output.append(int(value))
        except ValueError:
            output.append(value)
    output.append(region_index)
    region_index += 1        
    regions.append(output)

bed_file.close()
#print(regions)


#Enable sorting of elements and SNPs per chromosome
def chr_index(Chr):
    x = 0    
    if (Chr == "1" or Chr == "chr1"):
        x = 0
    if (Chr == "2" or Chr == "chr2"):
        x = 1
    if (Chr == "3" or Chr == "chr3"):
        x = 2
    if (Chr == "4" or Chr == "chr4"):
        x = 3
    if (Chr == "5" or Chr == "chr5"):
        x = 4
    if (Chr == "6" or Chr == "chr6"):
        x = 5
    if (Chr == "7" or Chr == "chr7"):
        x = 6		
    if (Chr == "8" or Chr == "chr8"):
        x = 7
    if (Chr == "9" or Chr == "chr9"):
        x = 8
    if (Chr == "10" or Chr ==  "chr10"):
        x = 9
    if (Chr == "11" or Chr ==  "chr11"):
        x = 10
    if (Chr == "12" or Chr ==  "chr12"):
        x = 11
    if (Chr == "13" or Chr ==  "chr13"):
        x = 12
    if (Chr == "14" or Chr ==  "chr14"):
        x = 13
    if (Chr == "15" or Chr ==  "chr15"):
        x = 14
    if (Chr == "16" or Chr ==  "chr16"):
        x = 15
    if (Chr == "17" or Chr ==  "chr17"):
        x = 16
    if (Chr == "18" or Chr ==  "chr18"):
        x = 17
    if (Chr == "19" or Chr ==  "chr19"):
        x = 18
    if (Chr == "20" or Chr ==  "chr20"):
        x = 19
    if (Chr == "21" or Chr ==  "chr21"):
        x = 20
    if (Chr == "22" or Chr == "chr22"):
        x = 21
    if (Chr == "X" or Chr == "chrX"):
        x = 22
    if (Chr == "Y" or Chr == "chrY"):
        x = 23
    return int(x)
    


#Open GWAS file and read in SNPs
GWAS_file = open(specify_GWAS_file, "r")
#next(GWAS_file)        # removes first line of GWAS_inactive for HeightAndBMI GWAs data

chromosomes = [[] for i in range(23)]
for line in GWAS_file:
    
    this_SNP = line.split()
    output= []
    chr_position= this_SNP[0].split(":")
    output.append(chr_position[0])

    try:
        chr_position[1] = int(chr_position[1])
    except ValueError:
        chr_position[1] = chr_position[1][:-1]
        chr_position[1] = int(chr_position[1])

    output.append(chr_position[1])
    output.append(this_SNP[1])
    
    chr = chr_index(output[0])
    chromosomes[chr].append(output)
    
GWAS_file.close()
print("finished reading GWAS file")
print(datetime.now())


#sort SNPs stored in each chromosomes by their nucleotide position
for chromosome in chromosomes:
    chromosome.sort(key=itemgetter(1))

print("finished sorting GWAS file")
print(datetime.now())


#obtain SNPs and their p-values for each element in the bed file
global_stats = []

for region in regions:
    #print(region)
    
    region_i = region[4]    #index of region
    chromosome = chr_index(region[0])   #get chromosome
    start = region[1]-1000   #get starting nucleotide of element
    end = region[2]+1000 #get end nucleotide of element
    length_region = region[2] - region[1] +2000  #get length of element
    i = 0   #new index
    start_SNP = 0   #start at first SNP
    end_SNP = 0 #new_end
    try:
        SNP = chromosomes[chromosome][0]    #retrieve first SNP from chromosome
    except IndexError:
        break
    SNP_position = SNP[1]   #positionOfFirstSNP in chromosome

    while end >= SNP_position:    #while end position of element is upstream of SNP
        i +=10000   #... jump 1000 SNPs further ...
        try:    
            SNP = chromosomes[chromosome][i]   #... and test this SNP ...
        except IndexError:
            SNP = chromosomes[chromosome][-1]   #if last SNP is reached take last one
            i = len(chromosomes[chromosome])-1
            break
        SNP_position = SNP[1]   #pass on updated SNP position
        end_SNP = i     #this is the last SNP in element

    while start <= SNP_position:
        i -=100
        try:
            SNP = chromosomes[chromosome][i]    
        except IndexError:
            break
        SNP_position = SNP[1]
        start_SNP = i
    
    N_SNPs = 0
    for SNP in chromosomes[chromosome][start_SNP:end_SNP]:
        SNP_position = int(SNP[1])
        if (SNP_position >= start and SNP_position <= end):
                global_stats.append(SNP)
                N_SNPs +=1
    #print(datetime.now())
    regions[region_i].append(N_SNPs)

print("finished with global stats")
print(datetime.now())

#Write result file for SNP associations within bed file elements
result_file = open(name_of_resultfile_for_tested_elements, "w")

for line in global_stats:
    for i in line:
        result_file.write(str(i) + "\t")
    result_file.write("\n")
result_file.close()

print("finished writing global stats")
print(datetime.now())



#define SNP serch function for permutation test
def call_SNPs(region):
    local_stats = []
    chromosome = chr_index(region[0])
    i = 0
    start_SNP = 0
    end_SNP = 0
    SNP = chromosomes[chromosome][0]    
    SNP_position = SNP[1]

    while region[2] >= SNP_position:
        i +=10000
        try:
            SNP = chromosomes[chromosome][i]    
        except IndexError:
            SNP = chromosomes[chromosome][-1]   
            i = len(chromosomes[chromosome])-1
            break
        SNP_position = SNP[1]
        end_SNP = i
        
    #print(str(SNP_position) + "define end")

    while region[1] <= SNP_position:
        i -=100
        try:
            SNP = chromosomes[chromosome][i]    
        except IndexError:
            break
        SNP_position = SNP[1]
        start_SNP = i
    #print(str(SNP_position) + "define start")
    
    N_SNPs = 0
    for SNP in chromosomes[chromosome][start_SNP:end_SNP]:
        SNP_position = int(SNP[1])
        if (SNP_position >= region[1] and SNP_position <= region[2]):
                local_stats.append(SNP)
                N_SNPs +=1
                
    #print(N_SNPs)
    #print(local_stats)
  
    return local_stats




#permutation
#Number of permutations should be > 1 000 (!), ideally 10 000
permut_stats = []
flex_region = 1000
number_of_permutations = 200

for i in range(number_of_permutations):
    new_regions = []
    
    for region in regions:
        #print("region" + str(region))
        chromosome = chr_index(region[0])
        try:
            maximum = chromosomes[chromosome][-1][1]
        except IndexError:
            break
        len_region = region[2]-region[1]
        end = region[2]
        N_SNPs = region[5]
        
        """
        new_end = random.randrange(0, maximum)
        new_start = new_end -len_region - flex_region
        new_region = [chromosome, new_start, new_end]
        new_SNPs = call_SNPs(new_region)
        """
        
        chromo_length = len(chromosomes[chromosome])
        new_start_SNP = random.randrange(0,chromo_length)
        new_start_SNP = chromosomes[chromosome][new_start_SNP]
        new_end = new_start_SNP[1]+ (len_region + flex_region / 2)
        new_start = new_end - len_region - flex_region
        new_region = [chromosome, new_start, new_end]
        new_SNPs = call_SNPs(new_region)
        
        while len(new_SNPs) <= N_SNPs:
                    new_start_SNP = random.randrange(0,chromo_length)
                    new_start_SNP = chromosomes[chromosome][new_start_SNP]
                    new_end = new_start_SNP[1]+ (len_region + flex_region / 2)
                    new_start = new_end - len_region - flex_region
                    new_region = [chromosome, new_start, new_end]
                    new_SNPs = call_SNPs(new_region)

        index = 0
        while index < N_SNPs:            
            new_regions.append(new_SNPs[index][2])
            index +=1

                
    permut_stats.append(new_regions)
    print("Permutation " + str(i))
    print(datetime.now())
           
 


#print permutation file
permut_result = open(name_of_resultfile_for_permutation, "w")

for permutation in permut_stats:
    for i in permutation:
        permut_result.write(i + "\t")
    permut_result.write("\n")
    
permut_result.close()

#calculate z-score for original global stats
z_scores_original = []
for SNP in global_stats:    
    p_value = SNP[2]
    z_score = fabs(ss.norm.ppf(float(p_value)/2))
    z_scores_original.append(z_score)
    
sum_of_z_scores = sum(z_scores_original)

zscore_result_original = open(zscore_of_tested_elements, "w")
zscore_result_original.write(str(sum_of_z_scores) + "\n")
zscore_result_original.close()    


#calculate z-score for each permutation
z_scores = []
for i, permutation in enumerate(permut_stats):    
    
    permutation_stats = []    
    for j, p_value in enumerate(permutation):
        
        z_score = fabs(ss.norm.ppf(float(permut_stats[i][j])/2))
        permutation_stats.append(z_score)
        
    sum_of_z_scores = sum(permutation_stats)
    z_scores.append(sum_of_z_scores)

zscore_result = open(zscore_of_permutations, "w")
for i in z_scores:
    zscore_result.write(str(i) + "\n")
zscore_result.close()    
    
print("Finished")