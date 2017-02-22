with open('/home/rantaplan/master/projectara/data/controlsmikraki.txt') as file:
    snpA= 'snp_0'#kanonika ta diavazei apo argsparse
    snpB= 'snp_2'        #!SOS MPOREI NA KANEI GREP KAI TO SN_8kati
                    #vazw se lista ta 2 snp 
    snps=[]
    for line in file:
        if line[:5]==snpA or line[:5]==snpB:
            line=line.rstrip('\n')
            snps.append(line) #ta vazw arxika se lista gia na kanw tin anazitisi me mia mono if(ara me ena mono skanarisma tou arxeiou)
#%%                YPOLOGISMOS DIPLOTUPWN
snpA_splitted=snps[0].split(' ')        
snpB_splitted=snps[1].split(' ')

RR=0 ; Rh=0 ; RA=0
hR=0 ; hh=0 ; hA=0
AR=0 ; Ah=0 ; AA=0  

for i in range(1, len(snpA_splitted) -len(snpA_splitted)%3 ,3):
    #oi gonotupoi tou kathe atomou stis 2 theseis 
    a_individual= snpA_splitted[i] + snpA_splitted[i+1] + snpA_splitted[i+2]
    b_individual= snpB_splitted[i] + snpB_splitted[i+1] + snpB_splitted[i+2]
    
    #diplotupoi: to elegxw ana grammi
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
#%%      YPOLOGISMOS pA pB
alleles=list(map(allele_freq, snps))
pA=alleles[0][1] 
pa=alleles[0][2]
pb=alleles[1][2]          
pB=alleles[1][1]  
#%%                     EM ALGORITHM
import math

pAB= pA * pB    #tsekare an einai swsti i arxiki timi

EMrun=True
count=0
while EMrun:
    count+=1
    E= 2*RR + Rh + hR + (pAB * (1 + pAB - pA - pB) * hh)/ (((pA - pAB) * (pB - pAB)) + pAB * (1 + pAB - pA - pB))
    
    pAB_MLE= E/1000    #kantw se sxesi me N   
    
    #pAB_new = (2*RR + Rh + hR + (pAB_MLE * (1 + pAB_MLE - pA - pB) * hh)/ (((pA - pAB_MLE) * (pB - pAB_MLE)) + pAB_MLE * (1 + pAB_MLE - pA - pB)))/1000
    
    print(count, abs(pAB_new - pAB))          
    if abs(pAB_new - pAB)<0.01:
        EMrun = False
    else:
        pAB=pAB_new
              
D=pAB_new- (pA*pB)
r2= D**2 /( pA * pa * pB * pb) #####EINAI MALLON LATHOS
if D <0:
    lista= [pA*pB, (1-pA) * (1-pB)]
    Dmax= min(lista)
    
else:
    lista=[pA*(1-pB), (1-pA)* pB]
    Dmax= min(lista)
          
Dtonos= D/Dmax
          
print(r2, D, Dtonos)
#%%
pAb_MLE= pA - pAB_MLE
    paB_MLE= pB - pAB_MLE
    pab_MLE= 1 - pAB_MLE -pA -pB
    


naB_new = 2*AR + hR + Ah + (hh * (pAb_new * paB_new) / ((pAb_new * paB_new) + (pAB_new* pab_new)))
nAB_new = 2*RR + Rh + hR + (hh * (pAB_new * pab_new) / ((pAb_new * paB_new) + (pAB_new * pab_new)))
nab_new = 2*AA + hA + Ah + (hh * (pAb_new * paB_new) / ((pAb_new * paB_new) + (pAB_new* pab_new)))
nAb_new = 2*RA + Rh + hA + (hh * (pAb_new * paB_new) / ((pAb_new * paB_new) + (pAB_new* pab_new)))

#==============================================================================
# L= (nAB * math.log(pAB)) + (naB * math.log(paB)) + (nAb * math.log(pAb)) + (nab * math.log(pab))
# L_new= (nAB_new * math.log(pAB_new)) + (naB_new * math.log(paB_new)) + (nAb_new * math.log(pAb_new)) + (nab_new * math.log(pab_new))
# 
#==============================================================================

#==============================================================================
# pAb= pA - pAB
# paB= pB - pAB
# pab= 1 - pAB -pA -pB
# naB = 2*AR + hR + Ah + (hh * (pAb * paB) / ((pAb * paB) + (pAB* pab)))
# nAB = 2*RR + Rh + hR + (hh * (pAB * pab) / ((pAb * paB) + (pAB* pab)))
# nab = 2*AA + hA + Ah + (hh * (pAb * paB) / ((pAb * paB) + (pAB* pab)))
# nAb = 2*RA + Rh + hA + (hh * (pAb * paB) / ((pAb * paB) + (pAB* pab)))
# 
#==============================================================================






























          
            
            
    