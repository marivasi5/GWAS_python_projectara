#Exoume treksei to tool sto arxiko dataset gwas.cases.gen kai gwas.controls.gen
#Apothikeusame ta snps me maf>0.05 sta arxeia: STUDYremovedMAFS_cases STUDYremovedMAFS_controls
#Trexoume to tool me auta ta arxeia ws input kai flag -HWE wste na upologistei to pvalue gia kathe thesi
#Output: STUDY.hwe

with open("path to -->'STUDY.hwe'")as file:
    #gia to plot
    X=[]
    Y=[]
    
    #gia to removal: apothikeusi se lista twn snp me pvalue < 0,001
    remove=[]
    
    for line in file:
        line=line.rstrip('\n')        
        splittedline= line.split(' ')
        
        if splittedline[1] == 'cannot':
            remove.append(splittedline[0])
        else:
            #apothikeuse stin lista twn X to snp
            snp_number= int(splittedline[0].split('_')[1])      #krataw mono ton arithmo tou snp (kai oxi to string px snp_142)    
            pvalue= round(float(splittedline[1]), 3)            #stroggulopoiw gia na min pethanw to plot        
            X.append(snp_number)
            Y.append(pvalue)
            
            #afairesi twn snp me pvalue<0.001
            if pvalue< 0.001:
                remove.append(splittedline[0])

#%%
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

fig_size = plt.rcParams["figure.figsize"]  # default diastaseis askonwn: 4, 6
fig_size[0] = 150 # dieurunsh aksona x
fig_size[1] = 10 # aksonas y
plt.rcParams["figure.figsize"] = fig_size 
ax.plot(X,Y, '.', c="purple") 
ax.set_xlabel("SNP ID", fontsize=15)
ax.set_ylabel("Pvalues ", fontsize=15)
ax.set_title("Hardy-Weinberg p-values Distribution", fontsize= 25)       


plt.savefig('hwe_plot.jpg')
plt.show()     
 #%%       