with open('/home/rantaplan/master/projectara/data/controlsmikraki.txt') as file:
    snpA= 'snp_3'
    snpB= 'snp_8'
                    #vazw se lista ta 2 snp 
    snps=[]
    for line in file:
        if line[:5]==snpA or line[:5]==snpB:
            line=line.rstrip('\n')
            snps.append(line) #ta vazw arxika se lista gia na kanw tin anazitisi me mia mono if(ara me ena mono skanarisma tou arxeiou)
#%%            
snpA_splitted=snps[0].split(' ')        
snpB_splitted=snps[1].split(' ')

RR=0 ; Rh=0 ; RA=0
hR=0 ; hh=0 ; hA=0
AR=0 ; Ah=0 ; AA=0  

for i in range(1, len(snpA_splitted) -len(snpA_splitted)%3 ,3):
    
    a_individual= snpA_splitted[i] + snpA_splitted[i+1] + snpA_splitted[i+2]
    b_individual= snpB_splitted[i] + snpB_splitted[i+1] + snpB_splitted[i+2]
    
    #to elegxw ana grammi
    if a_individual == '100':
        if b_individual == '100':
            RR+=1
        elif b_individual == '010':
            Rh+=1
        elif b_individual == '001':
            RA+=1
    elif a_individual == '010':
        if b_individual == '100':
            hR+=1
        elif b_individual == '010':
            hh+=1
        elif b_individual == '001':
            hA+=1
    elif a_individual == '001':
        if b_individual == '100':
            AR+=1
        elif b_individual == '010':
            Ah+=1
        elif b_individual == '001':
            AA+=1
            
#%%            
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
#%%     gia na vrw ta pA pB
alleles=list(map(allele_freq, snps))
pA=alleles[0][1]           
pB=alleles[1][1]            
            
            
    