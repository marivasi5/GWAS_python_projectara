import os

os.chdir('gwas')

f= open('100cases.txt')

for line in f: 

    line=line.rstrip('\n') 
         # if line: (gia an exw thema me tis kenes seires)
    splittedline= line.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
    snp= splittedline[0]
    
    homozygous_refrence=0
    homozygous_alternative=0
    heterozygous=0
    
    for i in range(5,len(splittedline) , 3): #!tsekare to range
        individual= splittedline[i] + splittedline[i+1] + splittedline[i+2] #epeidh einai lista to splittedline
        if individual=='100':
            homozygous_refrence+=1
        elif individual=='010':
            heterozygous+=1
        elif individual=='001':
            homozygous_alternative+=1
     
            
         #ypologismos suxnotitwn:
        p= (2*homozygous_refrence + heterozygous)/1000  #einai /2N opou N=500   
        q= (2*homozygous_alternative + heterozygous)/1000 #OR einai q=1-p  
    print(snp +' ' + str(p) +' '  + str(q))
    
#%%

def suxnotita(datasetLINE):    #prepei na to taiso lines    
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
        #ypologismos suxnotitwn:
                p= round((2*homozygous_refrence + heterozygous)/2000 , 3)  #einai /2N opou N=500   
                q= round((2*homozygous_alternative + heterozygous)/2000, 3)     
                    
                
    return snp+'\t'+ str(p) +'\t'+  str(q) +'\t'+ str(homozygous_refrence) + '\t' + str(heterozygous) + '\t' + str(homozygous_alternative)       #na dw ta rounds

    #an to eixa se return snp, p, q to ekane tuple! /an to eixa se print meta i epomeni print mmou ekane nera (none)  ##to return to vazei se tuple(!)

#%%

#Merged Dataset

with open('gwas.cases.gen') as cases, open('gwas.controls.gen') as controls:
    merged=open('merged.txt', 'w')
    for line_controls in controls:
        line_controls=line_controls.rstrip('\n')# string 
        for line_cases in cases:
            line_cases=line_cases.rstrip('\n')
            if line_cases[:5]==line_controls[:5]:
                print(line_cases,line_controls[-2999:], file=merged) #pros8esame sto line_cases ta line_controls ksekinwntas apo to telos ftasame sto stoixeio akrivws meta to alt allele (upologizontas kai ta spaces) kai phrame ola ta stoixeia tou string mexri to telos tou line
                break
            else:
                continue
            
    merged.close()
    
#%%
##Hardy-Weinberg


with open('merged.txt') as merged:
    for line in merged:
        mergedline=line.rstrip('\n')  # if line: (gia an exw thema me tis kenes seires)
        splitted_suxnotita= suxnotita(mergedline).split('\t') #dinei gia ka8e snp to p kai to q 
        p=splitted_suxnotita[1]
        q=splitted_suxnotita[2]
        expect_homozygous_refrence= (float(p)**2)*1000
        expect_homozygous_alternative= (float(q)**2)*1000
        expect_heterozygous= 2*float(p)*float(q)*1000
        hr=splitted_suxnotita[3]
        ha=splitted_suxnotita[5]
        he=splitted_suxnotita[4]
        if expect_homozygous_refrence!=0 and expect_homozygous_alternative!=0 and expect_heterozygous!=0: 
        # an k oi 3 ektimwmenes times einai mh midenikes painroume klassiko tupo x2
            x2test=((int(hr)-expect_homozygous_refrence)**2)/expect_homozygous_refrence+((int(he)-expect_heterozygous)**2)/expect_heterozygous+((int(ha)-expect_homozygous_alternative)**2)/expect_homozygous_alternative
       #an enas ap tous 2 ektimwmenous homozygous gonotypous einai mhdenikos tote kai o expected eterozygos mhdenizetai kai ara to x2 upologizetai mono ap ton allo omozygo  
        elif expect_homozygous_refrence==0:
            x2test=((int(ha)-expect_homozygous_alternative)**2)/expect_homozygous_alternative
        elif expect_homozygous_alternative==0:
            x2test=((int(hr)-expect_homozygous_refrence)**2)/expect_homozygous_refrence
        
        print(splitted_suxnotita[0],' ',x2test)

  

#%% 

#Genotypic Association

x2test

ORAA=
 
#==============================================================================
# 
#==============================================================================

