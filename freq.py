f=open('/home/rantaplan/master/togamatoproject/data/mikrakicases.txt')

for line in f: 
    line=line.rstrip('\n') 
#         # if line: (gia an exw thema me tis kenes seires)
    splittedline= line.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
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
        p= (2*homozygous_refrence + heterozygous)/1000  #einai /2N opou N=500   
        q= (2*homozygous_alternative + heterozygous)/1000 #OR einai q=1-p    
#  ----------------------------------------------------------   
#     #HW--------------------------------------------------------
#     #upologizw anamenomeno ARITHMO gonotupwn
#     expect_homozygous_refrence= (p**2)*500         
#     expect_homozygous_alternative= (q**2)*500 
#     expect_heterozygous= 2*p*q*500   
#     
#     xtest=((homozygous_refrence-expect_homozygous_refrence)**2)/expect_homozygous_refrence+((heterozygous-expect_heterozygous)**2)/expect_heterozygous+((homozygous_alternative-expect_heterozygous)**2)/expect_heterozygous
    print(snp +' ' + str(p) +' '  + str(q))
f.close()    
#%%
#SOS elegxos antistixisis tou x me tta cases kai tou y me ta control (kai oxi antistrofa)
def allele_frequencies(x,y,z):
    '''grapse ti kanei  x:cases y;Controls, z; onoma output'''
    
    output_name=open(z, 'w')
    with open(x) as cases, open(y) as controls: 
        for line_cases, line_controls in zip(cases, controls):
            splittedline_cases= line_cases.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
            splittedline_controls= line_controls.split(' ')          
            snp= splittedline_cases[0]
            
            #STA CASES
            CASEShomozygous_refrence=0
            CASEShomozygous_alternative=0
            CASESheterozygous=0
    
            for i in range(5, len(splittedline_cases), 3): #!tsekare to range
                individual= splittedline_cases[i] + splittedline_cases[i+1] + splittedline_cases[i+2]
                if individual=='100':
                    CASEShomozygous_refrence+=1
                elif individual=='010':
                    CASESheterozygous+=1
                elif individual=='001':
                    CASEShomozygous_alternative+=1
                
                     #ypologismos suxnotitwn:
                pCASES= (2*CASEShomozygous_refrence + CASESheterozygous)/1000  #einai /2N opou N=500   
                qCASES= 1-pCASES
            
            #STA CONTROLS
            CONTROLShomozygous_refrence=0
            CONTROLShomozygous_alternative=0
            CONTROLSheterozygous=0
    
            for i in range(5, len(splittedline_controls), 3): #!tsekare to range
                individual= splittedline_controls[i] + splittedline_controls[i+1] + splittedline_controls[i+2]
                if individual=='100':
                    CONTROLShomozygous_refrence+=1
                elif individual=='010':
                    CONTROLSheterozygous+=1
                elif individual=='001':
                    CONTROLShomozygous_alternative+=1
                
                     #ypologismos suxnotitwn:
                pCONTROLS= (2*CONTROLShomozygous_refrence + CONTROLSheterozygous)/1000  #einai /2N opou N=500   
                qCONTROLS= 1-pCONTROLS
            
            #KOITAKSE TO EDIT # kanto me return #PREPEI NA GRAFEI SE ARXEIO
            print(snp, round(pCONTROLS, 3) ,'\t', round(qCONTROLS, 3) ,'\t', round(pCASES,3) ,'\t', round(qCASES,3), '\t', round(pCASES+pCONTROLS ,3), '\t', round(qCONTROLS+qCONTROLS,3), file=output_name)
            output_name.close()  
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
                p= round((2*homozygous_refrence + heterozygous)/1000 , 3)  #einai /2N opou N=500   
                q= round(1-p , 3)     
                
    return snp+'\t'+ str(p) +'\t'+  str(q)        #na dw ta rounds

    #an to eixa se return snp, p, q to ekane tuple! /an to eixa se print meta i epomeni print mmou ekane nera (none)  
#%%
with open('/home/rantaplan/master/togamatoproject/data/mikrakicases.txt') as cases, open('/home/rantaplan/master/togamatoproject/data/mikrakicontrols.txt') as controls: 
    for line_cases, line_controls in zip(cases, controls):
        line_cases=line_cases.rstrip('\n')
        line_controls=line_controls.rstrip('\n')        #!to suxnotita(line_cases) exei type: NoneType kaii otan to kanw str() mou typwnei ena None san deuteri grammi
        
        print(suxnotita(line_cases), suxnotita(line_controls)[6:])        
#===================Ο ΑΛΛΟΣ ΤΡΟΠΟΣ ΜΕ ΤΙΣ 2 FOR===========================================================         
# with open('/home/rantaplan/master/togamatoproject/data/mikrakicases.txt') as cases, open('/home/rantaplan/master/togamatoproject/data/mikrakicontrols.txt') as controls: 
#     for line_controls in controls:
#         line_controls=line_controls.rstrip('\n') 
#         for line_cases in cases:
#             line_cases=line_cases.rstrip('\n')
#             if line_cases[:5]==line_controls[:5]:
#                 print(suxnotita(line_cases))
#                 break
#             else:
#                 continue
#             
#         continue                 
#==============================================================================
        

