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
#NEW FLAG!!!
parser.add_argument('-HWE', action='store_true')
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
    
#%%
##HWE:

def HWE(line_cases,line_controls):
    from scipy import stats
    counts_cases=genotype_counts(line_cases)                
    counts_controls=genotype_counts(line_controls)
    counts_merged= (counts_controls[0], counts_cases[1]+counts_controls[1], counts_cases[2]+counts_controls[2], counts_cases[3]+counts_controls[3], counts_cases[4]+counts_controls[4]) 
    if counts_merged[1] > 5 and counts_merged[2] > 5  and counts_merged[3] > 5: #This test is invalid when the observed or expected frequencies in each category are too small. A typical rule is that all of the observed and expected frequencies should be at least 5.
        snp, p_merged, q_merged = allele_freq(counts_merged)
        N = counts_merged[4]
        expect_homozygous_refrence = (p_merged**2)*N #pairnoume to plh8os twn anamenomenw reference
        expect_homozygous_alternative = (q_merged**2)*N
        expect_heterozygous = p_merged*q_merged*2*N
        x2test_hw = stats.chisquare([counts_merged[1], counts_merged[3],counts_merged[2]], [expect_homozygous_refrence, expect_heterozygous, expect_homozygous_alternative])
        return snp, x2test_hw.pvalue, round(expect_homozygous_refrence,3), round(expect_homozygous_alternative,3), round(expect_heterozygous,3)

#%%

if args.HWE:

    import time
    start = time.time()     
#Kanonika tha prepei na elegxw an to arxeio uparxei idi!!! if not os.path.exists('{'):
    j = 0
    
    output= open('{}.hwe'.format(args.output), 'w')
    with open(args.cases_file) as cases, open(args.controls_file) as controls:
        for line_cases,line_controls in zip(cases, controls):
                line_cases=line_cases.rstrip('\n')
                line_controls=line_controls.rstrip('\n')
                #gia elegxo taxutitas, mou deixnei poses fores exei treksei(ara poses grammes)
                if j % 2000 == 0:
                     print(j)   
                j += 1
                
                snp, x2test_hw,ehr,eha,ehet = HWE(line_cases,line_controls)
                print(snp, x2test_hw, file=output)
    print(time.time()-start)
                

