#GIA TREKSIMO
#python argtest.py -controls_file data/mikrakicontrols.txt -cases_file data/mikrakicases.txt -output testaki -allele_frequency
import os
import sys
import argparse

#einai apo to arxeio argtest. einai mono oi grammes pou xreiazontai gia na treksw to allelefreq
parser = argparse.ArgumentParser(description='This is our projectara')
parser.add_argument('-controls_file', help='File containing control genotypes',required=True)
parser.add_argument('-cases_file', help='File containing cases genotypes',required=True)
parser.add_argument('-output', help='Output file name',required=True)
parser.add_argument('-allele_frequency', action='store_true')           #!ftiakse ta help
args = parser.parse_args()
#%%
def genotype_counts(datasetLINE):    #prepei na to taiso lines    
    splittedline= datasetLINE.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
    
    snp= splittedline[0]
    homozygous_refrence=0
    homozygous_alternative=0
    heterozygous=0
    for i in range(5, len(splittedline), 3): #!tsekare to range
                individual= splittedline[i] + splittedline[i+1] + splittedline[i+2]
                if individual=='100':
                    homozygous_refrence+=1
                elif individual=='010':
                    heterozygous+=1
                elif individual=='001':
                    homozygous_alternative+=1
    return snp , homozygous_refrence, homozygous_alternative, heterozygous 
    #+' '+  str(homozygous_alternative) +' '+  str(homozygous_alternative) 
    #me tabs eixe thema        
#%%
def allele_freq(datasetLINE):
    #gia upologismo atomwn dataset (kathe grammi exei ton idio arithmo atomwn me tis alles vevaia)
    splittedline= datasetLINE.split(' ') #ena line twn 500: exei len 1505(=3*N+5)
    N=(len(splittedline)-5)/3
    
    p= round((genotype_counts(datasetLINE)[1]*2 + genotype_counts(datasetLINE)[3])/(2*N), 3)  #einai /2N opou N=500   
    q= round((genotype_counts(datasetLINE)[2]*2 + genotype_counts(datasetLINE)[3])/(2*N), 3)
    
    return splittedline[0], p, q
#%%             ΤΕΣΤ ΜΕ 1 ΑΡΧΕΙΟ
#==============================================================================
# with open('/home/rantaplan/master/togamatoproject/data/mikrakicases.txt') as file:
#     for line in file:
#         print(allele_freq(line))
#==============================================================================

#%%            TEST ME 2 ARXEIA
#==============================================================================
#  
# with open('/home/rantaplan/master/togamatoproject/data/mikrakicases.txt') as cases, open('/home/rantaplan/master/togamatoproject/data/mikrakicontrols.txt') as controls:
#     
#     for line_cases, line_controls in zip(cases, controls):
#                 line_cases=line_cases.rstrip('\n')
#                 line_controls=line_controls.rstrip('\n')        #!to suxnotita(line_cases) exei type: NoneType kaii otan to kanw str() mou typwnei ena None san deuteri grammi
#                 
#                 #to print einai anapiro giati to thelw se string   
#                 print(allele_freq(line_controls)[0], allele_freq(line_controls)[1], allele_freq(line_controls)[2], allele_freq(line_cases)[1], allele_freq(line_cases)[2], round(allele_freq(line_cases)[1]+allele_freq(line_controls)[1], 3), round(allele_freq(line_cases)[2]+ allele_freq(line_cases)[2], 3))             
#     
#==============================================================================
#%%
if args.allele_frequency:
    
    
#Kanonika tha prepei na elegxw an to arxeio uparxei idi!!! if not os.path.exists('{'):
    output= open('{}.frequency'.format(args.output), 'w')
    with open(args.controls_file) as cases, open(args.controls_file) as controls:
        for line_cases, line_controls in zip(cases, controls):
                 line_cases=line_cases.rstrip('\n')
                 line_controls=line_controls.rstrip('\n')        
                 
                 #to print einai anapiro giati to thelw se string   
                 print(allele_freq(line_controls)[0], allele_freq(line_controls)[1], allele_freq(line_controls)[2], allele_freq(line_cases)[1], allele_freq(line_cases)[2], round(allele_freq(line_cases)[1]+allele_freq(line_controls)[1], 3), round(allele_freq(line_cases)[2]+ allele_freq(line_cases)[2], 3), file=output)             
     








