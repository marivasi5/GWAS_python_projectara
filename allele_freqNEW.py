#GIA TREKSIMO
#python allele_freqNEW.py -controls_file data/controlsmikraki.txt -cases_file data/casesmikraki.txt -output testaki -allele_frequency

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
def genotype_counts(datasetLINE):        
    '''Ypologismos arithmou gonotupwn
    Input: grammi arxeiou se Genotype File Format 
    Output: tuple tis morfis (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn)'''
    
    splittedline= datasetLINE.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
    snp= splittedline[0]
    N=(len(splittedline)-5)/3   #Arithmos atomwn. px ena line twn 500 atomwn: exei len=1505(3*N+5)
    
    homozygous_refrence=0 ;homozygous_alternative=0 ; heterozygous=0
    for i in range(5, len(splittedline) -len(splittedline)%3, 3): #!tsekare to range
                individual= splittedline[i] + splittedline[i+1] + splittedline[i+2]
                if individual=='100':
                    homozygous_refrence+=1
                elif individual=='010':
                    heterozygous+=1
                elif individual=='001':
                    homozygous_alternative+=1
    return snp , homozygous_refrence, homozygous_alternative, heterozygous, N 
            
#%%                     TO SWSTO
def allele_freqLINE(datasetLINE):
    '''Ypologismos suxnotitas allilomorfwn (p=suxnotita refrence allilomorfou, q=suxnotita alternative allilomorfou)
    Input: grammi arxeiou se Genotype File Format 
    Output: tuple (snp_ID, p, q) '''
    
    snp, R, A, het, N = genotype_counts(datasetLINE)
    p= round((R*2 + het)/(2*N), 3)  
    q= round((A*2 + het)/(2*N), 3)
#===================TO LATHOS===========================================================
#     p= round((genotype_counts(datasetLINE)[1]*2 + genotype_counts(datasetLINE)[3])/(2*N), 3)  #einai /2N opou N=500   
#     q= round((genotype_counts(datasetLINE)[2]*2 + genotype_counts(datasetLINE)[3])/(2*N), 3)
#==============================================================================
    return snp, p, q
#%%                
def allele_freq(x):  #isws prepei na to treksoume etsi gia na vroume kai gia to merged xwris na exoume kanei merging 
    '''Ypologismos suxnotitas allilomorfwn (p=suxnotita refrence allilomorfou, q=suxnotita alternative allilomorfou)
    *Prepei prwta na treksi i genotype_counts*
    Input: tuple  (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn)
    Output: tuple (snp_ID, p, q) '''
    
    snp, R, A, het, N = x
    p= round((R*2 + het)/(2*N), 3)  
    q= round((A*2 + het)/(2*N), 3)
    return snp, p, q
#%%            TEST ME 2 ARXEIA
#==============================================================================
# import time
# start=time.time() 
# j = 0
#     
# with open('/home/rantaplan/master/projectara/data/gwas.cases.gen') as cases, open('/home/rantaplan/master/projectara/data/gwas.controls.gen') as controls:
#     
#     for line_cases, line_controls in zip(cases, controls):
#                 line_cases=line_cases.rstrip('\n')
#                 line_controls=line_controls.rstrip('\n')        #!to suxnotita(line_cases) exei type: NoneType kaii otan to kanw str() mou typwnei ena None san deuteri grammi
#                 if j % 2000 == 0:
#                      print(j)   
#                 j += 1
#                 
#                 counts_cases=genotype_counts(line_cases)                
#                 counts_controls=genotype_counts(line_controls)
#                 
#                 snp, p_controls, q_controls = allele_freq(counts_controls)
#                 snp, p_cases, q_cases = allele_freq(counts_cases)
# 
#     print(time.time()-start)     
# 
#==============================================================================
#%%
if args.allele_frequency:

    import time
    start = time.time()     
#Kanonika tha prepei na elegxw an to arxeio uparxei idi!!! if not os.path.exists('{'):
    j = 0
    
    output= open('{}.frequency'.format(args.output), 'w')
    with open(args.cases_file) as cases, open(args.controls_file) as controls:
        for line_cases,line_controls in zip(cases, controls):
                line_cases=line_cases.rstrip('\n')
                line_controls=line_controls.rstrip('\n')        
                 
                 #gia elegxo taxutitas, mou deixnei poses fores exei treksei(ara poses grammes)
                if j % 2000 == 0:
                     print(j)   
                j += 1
                 
#=======================ME LINES=======================================================
#                 snp, p_controls, q_controls = allele_freqLINE(line_controls)
#                 snp1, p_cases, q_cases = allele_freqLINE(line_cases)
#                  
#                 print(snp1, p_controls, q_controls, p_cases, q_cases, file=output)
#==============================================================================
                counts_cases=genotype_counts(line_cases)                
                counts_controls=genotype_counts(line_controls)
                counts_merged= (counts_controls[0], counts_cases[1]+counts_controls[1], counts_cases[2]+counts_controls[2], counts_cases[3]+counts_controls[3], counts_cases[4]+counts_controls[4]) 
                
                snp, p_controls, q_controls = allele_freq(counts_controls)
                snp, p_cases, q_cases = allele_freq(counts_cases)
                snp, p_merged, q_merged = allele_freq(counts_merged)
            
                print(snp, p_controls, q_controls, p_cases, q_cases, p_merged, q_merged, file=output)
    print(time.time()-start)     
#%%

