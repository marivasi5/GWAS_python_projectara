import os

os.chdir('gwas')


    
#%%
def genotype_counts(datasetLINE):    #prepei na to taiso lines    
    splittedline= datasetLINE.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
    
    snp= splittedline[0]
    homozygous_refrence=0
    homozygous_alternative=0
    heterozygous=0
    for i in range(5, len(splittedline) - len(splittedline)%3, 3): #!tsekare to range
                individual= splittedline[i] + splittedline[i+1] + splittedline[i+2]
                if individual=='100':
                    homozygous_refrence+=1
                elif individual=='010':
                    heterozygous+=1
                elif individual=='001':
                    homozygous_alternative+=1
    return snp , homozygous_refrence, homozygous_alternative, heterozygous 
    
  

#%%
def allele_freq(datasetLINE):
    #gia upologismo atomwn dataset (kathe grammi exei ton idio arithmo atomwn me tis alles vevaia)
    splittedline= datasetLINE.split(' ') #ena line twn 500: exei len 1505(=3*N+5)
    N =(len(splittedline)-5)/3
    
    p = round((2*genotype_counts(datasetLINE)[1] + genotype_counts(datasetLINE)[3])/(2*N), 3)  #einai /2N opou N=500   
    q = round(1-p,3)    
    return splittedline[0], p, q


#Kalw thn allele_freq gia na ftiakxw 2 listes,pou 8 ginoun oi times twn aksonwn g t plot minor_allele_freq

with open('gwas.cases.gen') as cases, open('gwas.controls.gen') as controls:
    minor_allele_freq_cases=[]
    minor_allele_freq_controls=[] 
    for line_cases, line_controls in zip(cases, controls):
         line_cases=line_cases.rstrip('\n')
         line_controls=line_controls.rstrip('\n')  
         minor_allele_freq_cases.append(allele_freq(line_cases)[2])
         minor_allele_freq_controls.append(allele_freq(line_controls)[2])#to print einai anapiro giati to thelw se string   
         #snps.append(allele_freq(line_controls)[0])
         #print(allele_freq(line_controls)[0], allele_freq(line_controls)[1], allele_freq(line_controls)[2], allele_freq(line_cases)[1], allele_freq(line_cases)[2], round(allele_freq(line_cases)[1]+allele_freq(line_controls)[1], 3), round(allele_freq(line_cases)[2]+ allele_freq(line_cases)[2], 3))             

#Plot minor allele freq - Gia ta mikra data set vgazei ok graph,sta megala gamietai!
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
#ax.plot(list(range(0,195145)), minor_allele_freq_cases,   mew=5, color='magenta',linestyle= "-",lw=1, alpha=1)#lw=line width, alpha=opacity
#ax.plot(list(range(0,195145)), minor_allele_freq_controls,  mew=5, color='cyan',linestyle= "-",lw=1, alpha=1)#lw=line width, alpha=opacity
X=list(range(0,195145))
Ycases=minor_allele_freq_cases
Ycontrols=minor_allele_freq_controls
legends = ax.plot(X,Ycases, 'b', X, Ycontrols, 'r')
print (legends)
plt.legend(legends, ["Cases", "Controls"], loc=1) # loc=2 shmainei panw aristera 
ax.set_xlabel("SNPs", fontsize=15)
ax.set_ylabel("Minor Allele Frequencies", fontsize=15)
ax.set_title("graaaaaaph", fontsize=25)
plt.show

   
#==============================================================================
# #%%

# 
# def suxnotita(datasetLINE):    #prepei na to taiso lines    
#     splittedline= datasetLINE.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
#     snp= splittedline[0]
#             
#     homozygous_refrence=0
#     homozygous_alternative=0
#     heterozygous=0
#     for i in range(5, len(splittedline), 3): #!tsekare to range
#                 individual= splittedline[i] + splittedline[i+1] + splittedline[i+2]
#                 if individual=='100':
#                     homozygous_refrence+=1
#                 elif individual=='010':
#                     heterozygous+=1
#                 elif individual=='001':
#                     homozygous_alternative+=1
#         #ypologismos suxnotitwn:
#                 p= round((2*homozygous_refrence + heterozygous)/2000 , 3)  #einai /2N opou N=500   
#                 q= round((2*homozygous_alternative + heterozygous)/2000, 3)     
#                     
#                 
#     return snp+'\t'+ str(p) +'\t'+  str(q) +'\t'+ str(homozygous_refrence) + '\t' + str(heterozygous) + '\t' + str(homozygous_alternative)       #na dw ta rounds
# 
#     #an to eixa se return snp, p, q to ekane tuple! /an to eixa se print meta i epomeni print mmou ekane nera (none)  ##to return to vazei se tuple(!)
# 
#==============================================================================
#%%
#==============================================================================
# #allele_frequency:
#     
# #def allele_frequency(x,y):
# 
#     with open('gwas.cases.gen') as cases, open('gwas.controls.gen') as controls: 
#             for line_cases, line_controls in zip(cases, controls):
#                 line_cases=line_cases.rstrip('\n')
#                 line_controls=line_controls.rstrip('\n') #!to suxnotita(line_cases) exei type: NoneType kaii otan to kanw str() mou typwnei ena None san deuteri grammi
#                 splitted_suxnotita_controls= suxnotita(line_controls).split('\t')
#                 splitted_suxnotita_cases= suxnotita(line_cases).split('\t')
#                 print(suxnotita(line_cases), splitted_suxnotita[1], '\t',splitted_suxnotita[2],'\t',splitted_suxnotita[2], '\t',splitted_suxnotita[3],splitted_suxnotita[4],'\t',splitted_suxnotita[5] )#einai lista
#      
# x=cases
# y=controls    
#==============================================================================
    

#%%

#Merged Dataset

#==============================================================================
#==============================================================================
# # def merger(x,y):
#     with open(x) as cases, open(y) as controls:
#         mergedlist=[] #to vzoume se lista gia na mh ftiaxnoume extra endiameso arxeio
#         for line_cases, line_controls in zip(cases, controls):
#             line_cases=line_cases.rstrip('\n')
#             line_controls=line_controls.rstrip('\n')
#             splitted_controls=line_controls.split(' ')
#             if line_cases.split(' ')[0]==splitted_controls[0]:
#                 line_all=line_cases
#                 for i in splitted_controls[5:]:
#                     line_all+= i+' '
#                 #print(line_all) # to len tou string  einai 6029 giati exoume 3029 xarakthres apo line_cases kai oi upoloipoi 3000 einai ta 500*3 atoma + ta kena                
#                 mergedlist.append(line_all)
#         return mergedlist
#==============================================================================
#==============================================================================
        
x='gwas.cases.gen'
y='gwas.controls.gen'

x='100cases.txt'
y='100controls.txt'


def merger(x,y):
    mergedlist=[] #to vαzoume se lista gia na mh ftiaxnoume extra endiameso arxeio
    with open(x) as cases, open(y) as controls:
        for line_cases, line_controls in zip(cases, controls):
            line_cases=line_cases.rstrip('\n')
            line_controls=line_controls.rstrip('\n')
            splitted_controls=line_controls.split(' ')
            if line_cases.split(' ')[0]==splitted_controls[0]:
                line_all=line_cases
                for i in splitted_controls[5:]:
                    line_all+= ' ' +i  #XRISTOS KAI PANAGIA 
                    #print(line_all) # to len tou string  einai 6029 giati exoume 3029 xarakthres apo line_cases kai oi upoloipoi 3000 einai ta 500*3 atoma + ta kena                
                mergedlist.append(line_all)
        return mergedlist

merger(x,y)
    
    
#%%
##HWE:

from scipy import stats


def HWE(dataline):    
    splittedline = dataline.split(' ') 
    N = (len(splittedline)-5)/3 #ena line twn 1000: exei len 3005(=3*N+5)
    expect_homozygous_refrence = (allele_freq(dataline)[1]**2)*N #pairnoume to plh8os twn anamenomenw reference
    expect_homozygous_alternative = (allele_freq(dataline)[2]**2)*N
    expect_heterozygous = allele_freq(dataline)[1]*allele_freq(dataline)[2]*2*N
    x2test_hw = stats.chisquare([genotype_counts(dataline)[1] , genotype_counts(dataline)[3], genotype_counts(dataline)[2]], [expect_homozygous_refrence, expect_heterozygous, expect_homozygous_alternative])
    return allele_freq(dataline)[0], x2test_hw, round(expect_homozygous_refrence,3), round(expect_homozygous_alternative,3), round(expect_heterozygous,3)


#==============================================================================
# test tou function HWE              
# for mergline in mergedlist:
#     print(HWE(mergline))   
#==============================================================================
#==============================================================================
#TEST gia HWE 
#from scipy import stats
# with open ('gwas.cases.gen') as q:
# 
#     for line in q:
#         expect_homozygous_refrence = (allele_freq(line)[1]**2)*500 #pairnoume to p
#         expect_homozygous_alternative = (allele_freq(line)[2]**2)*500
#         expect_heterozygous= allele_freq(line)[1]*allele_freq(line)[2]*2*500
#             #print(splitted_suxnotita[0], '\t', '\t',expect_homozygous_refrence, '\t',expect_homozygous_alternative, '\t',expect_heterozygous)
#             
#         x2test = stats.chisquare([genotype_counts(line)[1] , genotype_counts(line)[3], genotype_counts(line)[2]], [expect_homozygous_refrence, expect_heterozygous, expect_homozygous_alternative])
#             
#         print(allele_freq(line)[0],genotype_counts(line)[1],expect_homozygous_refrence,x2test)
#==============================================================================
    
    #==============================================================================
#         if expect_homozygous_refrence!=0 and expect_homozygous_alternative!=0 and expect_heterozygous!=0: 
#         # an k oi 3 ektimwmenes times einai mh midenikes painroume klassiko tupo x2
#             x2test=((int(hr)-expect_homozygous_refrence)**2)/expect_homozygous_refrence+((int(he)-expect_heterozygous)**2)/expect_heterozygous+((int(ha)-expect_homozygous_alternative)**2)/expect_homozygous_alternative
#        #an enas ap tous 2 ektimwmenous homozygous gonotypous einai mhdenikos tote kai o expected eterozygos mhdenizetai kai ara to x2 upologizetai mono ap ton allo omozygo  
#         elif expect_homozygous_refrence==0:
#             x2test=((int(ha)-expect_homozygous_alternative)**2)/expect_homozygous_alternative
#         elif expect_homozygous_alternative==0:
#             x2test=((int(hr)-expect_homozygous_refrence)**2)/expect_homozygous_refrence
#==============================================================================
#==============================================================================
#         
#         print(splitted_suxnotita[0],' ',x2test)
#==============================================================================

  

#%% 
#association_test:
x='gwas.cases.gen'
y='gwas.controls.gen'

x='100cases.txt'
y='100controls.txt'
#Genotypic Association --> http://www.gwaspi.org/?page_id=332


def association_test(x,y): # taizw datasets,gt einai genika argo!
    with open(x) as cases, open(y) as controls:
        pvalues=[]
        loci=[]
        for line_cases, line_controls in zip(cases, controls):
            line_cases=line_cases.rstrip('\n')
            line_controls=line_controls.rstrip('\n') #!to suxnotita(line_cases) exei type: NoneType kaii otan to kanw str() mou typwnei ena None san deuteri grammi
            splitted_controls= line_controls.split(' ')
            x2test_ga= stats.chisquare([genotype_counts(line_cases)[1],genotype_counts(line_controls)[1],genotype_counts(line_cases)[2],genotype_counts(line_controls)[2],genotype_counts(line_cases)[3],genotype_counts(line_controls)[3]], [HWE(line_cases)[2],HWE(line_controls)[2], HWE(line_cases)[3],HWE(line_controls)[3],HWE(line_cases)[4],HWE(line_controls)[4]]) 
            pvalues.append(x2test_ga[1])
            loci.append(splitted_controls[2])
            if genotype_counts(line_cases)[2]!= 0 and genotype_counts(line_controls)[1]!=0 and genotype_counts(line_controls)[3]!=0:
                OR_RRAA = (genotype_counts(line_cases)[1]*genotype_counts(line_controls)[2])/(genotype_counts(line_cases)[2]*genotype_counts(line_controls)[1])
                OR_RAAA = (genotype_counts(line_cases)[3]*genotype_counts(line_controls)[2])/(genotype_counts(line_cases)[2]*genotype_counts(line_controls)[3])
            else:
                OR_RRAA = "none"
                OR_RAAA = "none"
        return splitted_controls[0], splitted_controls[2], x2test_ga[1], OR_RRAA, OR_RAAA,pvalues,loci


        #print(association_test(line_cases,line_controls)[0:5])       

        
    
#==============================================================================
# ORAA-aa = (caseAA × ctrlaa) / (caseaa × ctrlAA)
# ORAa-aa = (caseAa × ctrlaa) / (caseaa × ctrlAa)

#%% 
#==============================================================================
#Manhattan plot

import matplotlib.pyplot as plt
import math


pvalues = list(map(lambda x:(-math.log(x)),association_test(x,y)[5])) #pvalues
loci = association_test(x,y)[6] #loci

plt.plot(loci, pvalues, ls='', marker='.')
plt.savefig("manhattan")
plt.show()

#sta kanonika datasets KANEI 125000 WRES NA TREKSEI!!
#%%
#qq plot

import numpy as np 
import pylab 

pv = np.asarray(pvalues)

stats.probplot(pvalues, dist = stats.exponnorm,sparams=(2.5,), plot=pylab)
#pylab.set_title("Probplot for exponential distr with shape parameter 2.5")
pylab.show()

#OR

import statsmodels.api as sm
import pylab
pv = np.asarray(pvalues)

sm.qqplot(pv, line='s')
pylab.show()


#%%Plotakiaaaa

#minor allele frequency plot












