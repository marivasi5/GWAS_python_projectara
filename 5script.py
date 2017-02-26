#Exoume treksei to tool sto arxiko dataset gwas.cases.gen kai gwas.controls.gen
#COMMAND: python PROJECTARA.py -controls_file data/gwas.controls.gen -cases_file data/gwas.cases.gen -output STUDY -allele_frequency
#Apomonwthikeisame ta snps me maf>0.05 sta arxeia: STUDYremovedMAFS_cases STUDYremovedMAFS_controls
#Trexoume to tool me auta ta arxeia ws input kai flag -HWE wste na upologistei to pvalue gia kathe thesi
#python PROJECTARA.py -controls_file STUDYremovedMAFS_controls -cases_file STUDYremovedMAFS_cases -output STUDY -HWE
#Apothikeusame ta snps me pvalue > 0,001 apo to -HWE sta arxeia: STUDYremovedHWE_cases  STUDYremovedHWE_contros
#Trexoume to tool me ta parapanw arxeia ws input kai flag -association test -manhattan -qqplot
#COMMAND: python PROJECTARA.py -controls_file STUDYremovedHWE_controls -cases_file STUDYremovedHWE_cases -output STUDY -association_test -manhattan -qqplot
#Output file: STUDY.association

#%%             Apomonwsi statistical significant snps 
with open('STUDY.association')as file:
    remove=[]
    for line in file:
        line=line.rstrip('\n')        
        splittedline= line.split(' ')
        
        if float(splittedline[2]) < 1e-06:    
            remove.append(splittedline[0])

with open('STUDYsignificant_cases', 'w') as output_cases, open('STUDYremovedHWE_cases') as genotype_file_cases:
    for line in genotype_file_cases:
        if line.split(' ')[0] in remove:
            output_cases.write(line)
            
with open('STUDYsignificant_controls', 'w') as output_controls, open('STUDYremovedHWE_controls') as genotype_file_controls:
    for line in genotype_file_controls:
        if line.split(' ')[0] in remove:
            output_controls.write(line)
                       
#%%              Entopismos thesewn pou vriskontai se sundesi me ta statistical significant snps

#Gia tin analusi xrisimopoioume sunartiseis pou einai meros tou tool mas wstoso einai aparaitites k prosarmosmenes stis anagkes tou sugkekrimenou script
def genotype_counts(line):        
    '''Ypologismos arithmou gonotupwn
    Input: grammi arxeiou se Genotype File Format 
    Output: tuple tis morfis (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn, thesi SNP sumfwna me to NCBI build 36)'''
    splittedline= line.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
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
def allele_freq(x):   
    '''Ypologismos suxnotitas allilomorfwn (p=suxnotita refrence allilomorfou, q=suxnotita alternative allilomorfou)
    *Prepei prwta na treksi i genotype_counts*
    Input: tuple tis morfis (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn, thesi SNP sumfwna me to NCBI build 36)
    Output: tuple (snp_ID, p, q) '''
    
    snp, R, A, het, N, loci = x
    p= round((R*2 + het)/(2*N), 3)  
    q= round((A*2 + het)/(2*N), 3)
    return snp, p, q
def diplotupoi(snplist):
    '''Ypologismos twn diplotupwn pou xreiazontai gia na pragmatopoih8ei to haplotype estimation wste na upologistei to LD 2 thesewn
    Input: lista me 2 antikeimena: grammes twn 2 thesewn se Genotype File Format  
    Output: tuple tis morfis (arithmos atomwn omozugwn gia to refrence kai stis 2 theseis, arithmos atomwn omozugwn gia to refrence stin prwti thesi kai eterozugwn stin deuteri, arithmos atomwn eterozugwn stin prwti thesi kai omozugwn gia to refrence stin deuteri, arithmos atomwn eterozugwn kai stis duo theseis, arithmos atomwn'''    
    snpA, snpB=snplist    
    snpA_splitted=snpA.split(' ')        
    snpB_splitted=snpB.split(' ')
    N=(len(snpA_splitted)-5)/3
    RR=0 ; Rh=0 ; 
    hR=0 ; hh=0 ; 
    for i in range(5, len(snpA_splitted) -len(snpA_splitted)%3 ,3):
        #oi gonotupoi tou kathe atomou stis 2 theseis (theseis a kai p) 
        a_individual= snpA_splitted[i] + snpA_splitted[i+1] + snpA_splitted[i+2]
        b_individual= snpB_splitted[i] + snpB_splitted[i+1] + snpB_splitted[i+2]
        
        if a_individual == '100':
            if b_individual == '100':
                RR+=1
            elif b_individual == '010':
                Rh+=1
        elif a_individual == '010':
            if b_individual == '100':
                hR+=1
            elif b_individual == '010':
                hh+=1
    return RR, Rh, hR, hh, N
def EM(pAB,pA, pB, diplotype): #diplotupe: einai tuple apo tin sunartisi 
    '''Ypologismos tis ektimomenis aplotupikis suxnotitas ektelwntas ton EM algorithmo gia dyo snp
    *Prepei prwta na exei treksie i sunartisi diplotupoi gia ta snp pou theloume na meletisoume*
    Input: arxiki timi gia tin suxnotita tou aplotupou gia ton EM, suxnotita refrence allilomorfou stin mia thesi, suxnotita refrence allilomorfou stin alli thesi, tuple tis morfis (arithmos atomwn omozugwn gia to refrence kai stis 2 theseis, arithmos atomwn omozugwn gia to refrence stin prwti thesi kai eterozugwn stin deuteri, arithmos atomwn eterozugwn stin prwti thesi kai omozugwn gia to refrence stin deuteri, arithmos atomwn eterozugwn kai stis duo theseis    
    Output: ektimomeni aplotupiki suxnotita'''
    EMrun=True
    RR, Rh, hR, hh , N= diplotype
    
    while EMrun:
        E= 2*RR + Rh + hR + (pAB * (1 + pAB - pA - pB) * hh)/ (((pA - pAB) * (pB - pAB)) + pAB * (1 + pAB - pA - pB))
        pAB_new= E/(2*N)    
        
        if abs(pAB_new - pAB)< 10**(-6):
            EMrun = False
        else:
            pAB=pAB_new
    return pAB_new
   
#%%             Reverse lines (ara kai snps) gia tin downstream analysi

with open('/home/rantaplan/master/projectara/STUDYremovedHWE_cases') as file:
    listara=file.readlines()
with open('STUDYremovedHWE_cases_reversed', 'w') as output:
    listara.reverse()
    for line in listara:
        output.write(line)
     
#%%
output=open('significant_LD', 'w')

for i in remove:    
    snpA= i
    with open('/home/rantaplan/master/projectara/STUDYremovedHWE_cases') as file:
        for line in file:
            if line.split(' ')[0]==snpA:
                line=line.rstrip('\n')
                A=line
                B=next(file)
                B=B.rstrip('\n')
    
                counter=0
                saved=[]
                while counter < 1000:
                    snps= [A, B]
    
                    #Ypologismos pA pB
                    counts=list(map(genotype_counts, snps))
                    alleles=list(map(allele_freq, counts))
                    pA=alleles[0][1] ;pa=alleles[0][2]
                    pb=alleles[1][2] ;pB=alleles[1][1] 
    
                    import random
                    if pA>0.99 or pB>0.99 or pa>0.99 or pb>0.99:        
                        D=0
                        Dtonos=0
                    else:    
                        diplotype= diplotupoi(snps)
                        pABs=[random.uniform(0,1) for x in range(100)]
                        pAB_newS=[]
                        for i in pABs:
                            pAB_new= EM(i, pA, pB, diplotype)
                            pAB_newS.append(pAB_new)
                        pAB_mean= sum(pAB_newS)/len(pAB_newS)    
                        D= pAB_mean - (pA * pB)
                    
                        #Ypologismos Dtonos
                        if D <0:
                            lista= [pA*pB, (1-pA) * (1-pB)]
                            Dmax= min(lista)                        
                        else:
                            lista=[pA*(1-pB), (1-pA)* pB]
                            Dmax= min(lista)
                        Dtonos= D/Dmax                     
                    
                    if Dtonos > 0.9:
                        saved.append(B.split(' ')[0])
                        counter = 0
                    else:
                        counter += 1
                
                    B=next(file)
                    B=B.rstrip('\n')        

#%%                         DOWNSTREAM
    
    #trexoume tin idia diadikasia sto reversed (ara oustiastika kanoume downstream analusi)
    with open('/home/rantaplan/master/projectara/STUDYremovedHWE_cases_reversed') as file:
        for line in file:
            if line.split(' ')[0]==snpA:
                line=line.rstrip('\n')
                A=line
                B=next(file)
                B=B.rstrip('\n')
    
                counter=0
                while counter < 1000:
                    snps= [A, B]
    
                    #Ypologismos pA pB
                    counts=list(map(genotype_counts, snps))
                    alleles=list(map(allele_freq, counts))
                    pA=alleles[0][1] ;pa=alleles[0][2]
                    pb=alleles[1][2] ;pB=alleles[1][1] 
    
                    import random
                    if pA>0.99 or pB>0.99 or pa>0.99 or pb>0.99:        
                        D=0
                        Dtonos=0
                    else:    
                        diplotype= diplotupoi(snps)
                        pABs=[random.uniform(0,1) for x in range(100)]
                        pAB_newS=[]
                        for i in pABs:
                            pAB_new= EM(i, pA, pB, diplotype)
                            pAB_newS.append(pAB_new)
                        pAB_mean= sum(pAB_newS)/len(pAB_newS)    
                        D= pAB_mean - (pA * pB)
                    
                        #Ypologismos Dtonos
                        if D <0:
                            lista= [pA*pB, (1-pA) * (1-pB)]
                            Dmax= min(lista)                        
                        else:
                            lista=[pA*(1-pB), (1-pA)* pB]
                            Dmax= min(lista)
                        Dtonos= D/Dmax                     
                    
                    if Dtonos > 0.9:
                        saved.append(B.split(' ')[0])
                        counter = 0
                    else:
                        counter += 1
                
                    B=next(file)
                    B=B.rstrip('\n')         
                    
    print('The {} variants in linkage disequilibrium with {} are: {} \n'.format(len(saved), snpA, saved), file=output)
    print('euge kamari')
output.close()                
                
#%%                    INFO gia ta significant snps

import requests
from pyliftover import LiftOver     #gia to liftover
lo=LiftOver('hg18', 'hg38')         # downloads the hg18-to-hg38 coordinate conversion chain file from UCSC


#%%   Apothikevw se lista tis plirofories pou xreiazetai i get gia kathe significant snp

with open('Path to --> STUDYremovedHWE_cases') as cases:
    significant=[]  
    for line in cases:
              info=[]
              splittedline= line.split(' ')
              if splittedline[0] in remove:
                    info.append(splittedline[0])    #snp_ID                  
                    info.append(splittedline[2])    #position
                    info.append(splittedline[3])    #refrence
                    info.append(splittedline[4])    #alternative
                    significant.append(info) 

#%%
requestara=[]
for i in significant:
    info=[i[0]]    
    position=int(i[1])
    refrence=i[2]
    alternative=i[3]
        
    #LIFTOVER    
    conversion= lo.convert_coordinate('chr20', position)
    #The result is either None (if the source chromosome name is unrecognized) or a list of target positions in the new assembly. The list may be empty (locus is deleted in the new assembly), have a single element (locus matched uniquely)
    try:
        if len(conversion)==1:
            position_new=conversion[0][1]
            #GET REQUEST
            url= 'http://rest.ensembl.org/vep/human/hgvs/20:g.' + str(position_new) + refrence + '>' + alternative + '?' #to xromoswma twra einai to 20 alla tha eprepe na to vazw
            headers={ "Content-Type" : "application/json"}
            r = requests.get(url, headers=headers)
                        
            if r.ok:
                data = r.json()
                info.append(data)
                requestara.append(info)
            else:    
                print('The request for {} was not successful.'.format(i[0]))    
        else:
            print('The request for {} was not successful. Could not convert position to USCS hg38.'.format(i[0]))
    
    except TypeError:       #object of type 'NoneType' has no len()
        print('The request for {} was not successful. Could not convert position to USCS hg38.'.format(i[0]))
    

for i in requestara:
    print('The information acquired from Ensembl\'s Variant Effect Predictor for {} is: \n {}'.format(i[0], i[1]))    
    
           

            
            
            
            
            
