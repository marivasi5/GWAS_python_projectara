def genotype_counts(datasetLINE):    #prepei na to taiso lines    
    splittedline= datasetLINE.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
    
    snp= splittedline[0]
    homozygous_refrence=0
    homozygous_alternative=0
    heterozygous=0
    for i in range(5, len(splittedline) -len(splittedline)%3, 3): #!tsekare to range
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
#%%
def merger(x,y):
    
    with open(x) as cases, open(y) as controls:
        mergedlist=[] #to vzoume se lista gia na mh ftiaxnoume extra endiameso arxeio
        for line_cases, line_controls in zip(cases, controls):
            line_cases=line_cases.rstrip('\n')
            line_controls=line_controls.rstrip('\n')
            splitted_controls=line_controls.split(' ')
            if line_cases.split(' ')[0]==splitted_controls[0]:
                line_all=line_cases
                for i in splitted_controls[5:]:
                    line_all+= i+ ' '  #XRISTOS KAI PANAGIA 
                #print(line_all) # to len tou string  einai 6029 giati exoume 3029 xarakthres apo line_cases kai oi upoloipoi 3000 einai ta 500*3 atoma + ta kena                
                mergedlist.append(line_all)
        return mergedlist

#%%
         #   TREKSIMO
merger('/home/rantaplan/master/togamatoproject/data/mikrakicases.txt', '/home/rantaplan/master/togamatoproject/data/mikrakicontrols.txt')

#%%
from scipy import stats
with open ('/home/rantaplan/master/togamatoproject/data/mikrakicases.txt') as q:

    for line in q:          #SOS NA ALLAZOUME TA N
        splittedline= line.split(' ') #ena line twn 500: exei len 1505(=3*N+5)
        N=(len(splittedline)-5)/3
        
        expect_homozygous_refrence = (allele_freq(line)[1]**2)*N #pairnoume to p
        expect_homozygous_alternative = (allele_freq(line)[2]**2)*N
        expect_heterozygous= allele_freq(line)[1]*allele_freq(line)[2]*2*N
            #print(splitted_suxnotita[0], '\t', '\t',expect_homozygous_refrence, '\t',expect_homozygous_alternative, '\t',expect_heterozygous)
            
        x2test = stats.chisquare([genotype_counts(line)[1] , genotype_counts(line)[3], genotype_counts(line)[2]], [expect_homozygous_refrence, expect_heterozygous, expect_homozygous_alternative])
        print(x2test)



#%%

#Genotypic Association --> http://www.gwaspi.org/?page_id=332

with open('gwas.cases.gen') as cases, open('gwas.controls.gen') as controls: 
    for line_cases, line_controls in zip(cases, controls):
        line_cases=line_cases.rstrip('\n')
        line_controls=line_controls.rstrip('\n') #!to suxnotita(line_cases) exei type: NoneType kaii otan to kanw str() mou typwnei ena None san deuteri grammi
        splitted_suxnotita_controls= suxnotita(line_controls).split('\t')
        splitted_suxnotita_cases= suxnotita(line_cases).split('\t')
        #x2test= stats.chisquare([obsCaseAA+obsControlAA,obsCaseaa+obsControlaa,obsCaseAa+obsCaseAa], [expCaseAA+expControlAA,expCaseaa+expControlaa,expCaseAa+expControlAa])
        if float(splitted_suxnotita_cases[5])!=0 and float(splitted_suxnotita_controls[3])!=0 and float(splitted_suxnotita_cases[3])!=0 and float(splitted_suxnotita_controls[4])!=0:
            ORAACC = float(splitted_suxnotita_cases[3])*float(splitted_suxnotita_controls[5])/float(splitted_suxnotita_cases[5])*float(splitted_suxnotita_controls[3])
            ORACCC = float(splitted_suxnotita_cases[4])*float(splitted_suxnotita_controls[3])/float(splitted_suxnotita_cases[3])*float(splitted_suxnotita_controls[4])
        else:
            pass
            print(splitted_suxnotita_cases[0],ORAACC,ORACCC)
        


           

#==============================================================================
# ORAA-aa = (caseAA × ctrlaa) / (caseaa × ctrlAA)
# ORAa-aa = (caseAa × ctrlaa) / (caseaa × ctrlAa)

#==============================================================================


x2test = ((obsCaseAA-expCaseAA)² / expCaseAA) +
    ((obsCaseAa-expCaseAa)² / expCaseAa) +
    ((obsCaseAa-expCaseAa)² / expCaseaa) +
    ((obsControlAA-expControlAA)² / expControlAA) +
    ((obsCaseAa-expControlAa)² / expControlAa) +
    ((obsControlaa-expControlaa)² / expControlaa)
