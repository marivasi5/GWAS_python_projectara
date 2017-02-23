#%%
import os

os.chdir('gwas')
#%%
#%%
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
parser.add_argument('-association_test', action='store_value')
args = parser.parse_args()




#%%
def genotype_counts(datasetLINE):        
    '''Ypologismos arithmou gonotupwn
    Input: grammi arxeiou se Genotype File Format 
    Output: tuple tis morfis (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn, thesi SNP)'''
    
    splittedline= datasetLINE.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
    snp= splittedline[0]
    N=(len(splittedline)-5)/3   #Arithmos atomwn. px ena line twn 500 atomwn: exei len=1505(3*N+5)
    locus=splittedline[2]
    
    homozygous_refrence=0 ;homozygous_alternative=0 ; heterozygous=0
    for i in range(5, len(splittedline) -len(splittedline)%3, 3): #!tsekare to range
                individual= splittedline[i] + splittedline[i+1] + splittedline[i+2]
                if individual=='100':
                    homozygous_refrence+=1
                elif individual=='010':
                    heterozygous+=1
                elif individual=='001':
                    homozygous_alternative+=1
    return snp , homozygous_refrence, homozygous_alternative, heterozygous, N , locus
            
#%%                
def allele_freq(x):  #isws prepei na to treksoume etsi gia na vroume kai gia to merged xwris na exoume kanei merging 
    '''Ypologismos suxnotitas allilomorfwn (p=suxnotita refrence allilomorfou, q=suxnotita alternative allilomorfou)
    *Prepei prwta na treksi i genotype_counts*
    Input: tuple  (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn)
    Output: tuple (snp_ID, p, q) '''
    
    snp, R, A, het, N, loci = x
    p= round((R*2 + het)/(2*N), 3)  
    q= round((A*2 + het)/(2*N), 3)
    return snp, p, q
    
#%%


def HWE(x,y):# x=counts_cases, y=counts_controls
    from scipy import stats
    
    counts_merged= (y[0], x[1]+y[1], x[2]+y[2], x[3]+y[3], x[4]+y[4], x[5]) 
    
    if counts_merged[1] > 0 and counts_merged[2] > 0  and counts_merged[3] > 0: #This test is invalid when the observed or expected frequencies in each category are too small. A typical rule is that all of the observed and expected frequencies should be at least 5.
        snp, p_merged, q_merged = allele_freq(counts_merged)
        N = counts_merged[4]
        expect_homozygous_refrence = (p_merged**2)*N #pairnoume to plh8os twn anamenomenw reference
        expect_homozygous_alternative = (q_merged**2)*N
        expect_heterozygous = p_merged*q_merged*2*N
        x2test_hw = stats.chisquare([counts_merged[1], counts_merged[3],counts_merged[2]], [expect_homozygous_refrence, expect_heterozygous, expect_homozygous_alternative])
    
        return snp, x2test_hw.pvalue, round(expect_homozygous_refrence,3), round(expect_homozygous_alternative,3), round(expect_heterozygous,3)
    else:
        snp = counts_merged[0]
        return snp, "cannot compute pvalue", "ehr", "eha", "ehet"
#%%
#%%
#Genotypic Association --> http://www.gwaspi.org/?page_id=332

def association_test(x,y): # taizw tuples tou genotype_counts # x=counts_cases, y=counts_controls
            from scipy import stats            
            snp = x[0]
            loci = x[5]
            if x[1] > 0 and y[1] > 0 and x[2] >0 and y[2] > 0 and x[3] > 0 and y[3] > 0:
                snp,x2test_hw,ehr,eha,ehet = HWE(x,y)
                snp, p_cases, q_cases= allele_freq(x)
                snp, p_controls, q_controls= allele_freq(y)
                N=x[4]
                x2test_ga = stats.chisquare([x[1],y[1],x[2],y[2],x[3],y[3]], [(p_cases**2)*N,(p_controls**2)*N, (q_cases**2)*N, (q_controls**2)*N ,2*p_cases*q_cases*N, 2*p_controls*q_controls*N]) 
                x2test_ga_final= x2test_ga.pvalue
            else:
                x2test_ga_final = "cannot compute pvalue"
            if x[2]!= 0 and y[1]!=0 and y[3]!=0:
    
                OR_RRAA = (x[1]*y[2])/(x[2]*y[1])
                OR_RAAA = (x[3]*y[2])/(x[2]*y[3])
            else:
                OR_RRAA = "none"
                OR_RAAA = "none"

            return snp, loci, x2test_ga_final, OR_RRAA, OR_RAAA
            
#%%                 HWE treksimo

with open('gwas.cases.gen') as cases, open('gwas.controls.gen') as controls:
    j=0
    
    for line_cases, line_controls in zip(cases, controls):
        
                line_cases=line_cases.rstrip('\n')
                line_controls=line_controls.rstrip('\n')
                #gia elegxo taxutitas, mou deixnei poses fores exei treksei(ara poses grammes)
                
                if j % 2000 == 0:
                     print(j)   
                j += 1
                
                counts_cases=genotype_counts(line_cases)                
                counts_controls=genotype_counts(line_controls)
                
                snp, x2test_hw,ehr,eha,ehet = HWE(counts_cases,counts_controls)
                print(snp, x2test_hw)

#%%                 ASSOCIATION TREKSIMO
with open('gwas.cases.gen') as cases, open('gwas.controls.gen') as controls:
    j=0
    loci_list = []
    pvalue_list = []
    
    for line_cases, line_controls in zip(cases, controls):
        
                line_cases=line_cases.rstrip('\n')
                line_controls=line_controls.rstrip('\n')
                #gia elegxo taxutitas, mou deixnei poses fores exei treksei(ara poses grammes)
                
                if j % 2000 == 0:
                     print(j)   
                j += 1
                
                counts_cases=genotype_counts(line_cases)                
                counts_controls=genotype_counts(line_controls)
                
                snp, loci, pvalue, OR_RRAA, OR_RAAA  = association_test(counts_cases,counts_controls)
                loci_list.append(loci)
                pvalue_list.append(pvalue)
                #print(snp, loci, pvalue, OR_RRAA, OR_RAAA)






#%% 
#==============================================================================
#Manhattan plot

import matplotlib.pyplot as plt
import math

i=0
length=len(pvalue_list)
while i <= length:
    
    if pvalue_list[i] == 'cannot compute pvalue':
        del pvalue_list[i]
        del loci_list[i]
        i-=1 #meiwnw mesa sto if to i kata 1 wste na elegxw to epomeno stoixeio p mpainei sth 8esh tou string
        
    i+=1# auksanw kata 1 wste na proxwraw ston elegxo ths listas
    length=len(pvalue_list)
    print (length)
    
    
#pvaluelist_new = list(map(lambda x:x if x!='cannot compute pvalue',pvalue_list))

pvalues = list(map(lambda x:(-math.log(x)),pvalue_list))

#fig,ax = plt.subplots()
plt.plot(loci_list, pvalues, ls='', marker='.')
plt.savefig("manhattan")
plt.show()
#%%

        



