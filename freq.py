f=open('/home/rantaplan/master/togamatoproject/data/controls10K.txt')

for line in f: 

    line=line.rstrip('\n') 
#         # if line: (gia an exw thema me tis kenes seires)
    splittedline= line.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
    snp= splittedline[0]
    
    homozygous_refrence=0
    homozygous_alternative=0
    heterozygous=0
    suxnotites=[]
    for i in range(5, len(splittedline), 3): #!tsekare to range
        individual= splittedline[i] + splittedline[i+1] + splittedline[i+2] #individual = splittedline[i:i+3]
        if individual=='100':
            homozygous_refrence+=1
        elif individual=='010':
            heterozygous+=1
        elif individual=='001':
            homozygous_alternative+=1
#==============================================================================
#             
#         #ypologismos suxnotitwn:
        p= (2*homozygous_refrence + heterozygous)/1000  #einai /2N opou N=500   
        q= (2*homozygous_alternative + heterozygous)/1000 #OR einai q=1-p    
#     
#     #HW
#     #upologizw anamenomeno ARITHMO gonotupwn
#     expect_homozygous_refrence= (p**2)*500         
#     expect_homozygous_alternative= (q**2)*500 
#     expect_heterozygous= 2*p*q*500   
#     
#     xtest=((homozygous_refrence-expect_homozygous_refrence)**2)/expect_homozygous_refrence+((heterozygous-expect_heterozygous)**2)/expect_heterozygous+((homozygous_alternative-expect_heterozygous)**2)/expect_heterozygous
    print(snp +' ' + str(p) +' '  + str(q))
#==============================================================================
f.close()    









    #############http://stackoverflow.com/questions/11295171/read-two-textfile-line-by-line-simultaneously-python
#==============================================================================
# f=open('/home/rantaplan/master/togamatoproject/data/cases10K.txt')
# 
# for line in f: 
# 
#     line=line.rstrip('\n') 
# #         # if line: (gia an exw thema me tis kenes seires)
#     splittedline= line.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
#     snp= splittedline[0]
#     
#     homozygous_refrence=0
#     homozygous_alternative=0
#     eterozygous=0
#     
#     for i in range(5, len(splittedline), 3): #!tsekare to range
#         individual= splittedline[i] + splittedline[i+1] + splittedline[i+2]
#         if individual=='100':
#             homozygous_refrence+=1
#         elif individual=='010':
#             eterozygous+=1
#         elif individual=='001':
#             homozygous_alternative+=1
#             
#         #ypologismos suxnotitwn:
#         refrence_allele_frequency= (2*homozygous_refrence + eterozygous)/1000  #einai /2N opou N=500   
#         alternative_allele_frequency= (2*homozygous_alternative + eterozygous)/1000     
#     
#  #allakse to   print(snp +' ' + str(refrence_allele_frequency) +' '  + str(alternative_allele_frequency))
# f.close()    
#     
# 
#==============================================================================
