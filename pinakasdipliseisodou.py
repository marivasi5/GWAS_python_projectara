with open('/home/rantaplan/master/projectara/data/controlsmikraki.txt') as file:
    snpA= 'snp_3'
    snpB= 'snp_8'
                    #vazw se lista ta 2 snp 
    snps=[]
    for line in file:
        if line[:5]==snpA or line[:5]==snpB:
            line=line.rstrip('\n')
            snps.append(line)
#%%            
snpA_splitted=snps[0].split(' ')        #ta vazw arxika se lista gia na kanw tin anazitisi me mia mono if(ara me ena mono skanarisma tou arxeiou)
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
            
            
            
            
            
            
    