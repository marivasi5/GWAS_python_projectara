#Exoume treksei to tool sto arxiko dataset gwas.cases.gen kai gwas.controls.gen
#Apothikeusame ta snps me maf>0.05 sta arxeia: STUDYremovedMAFS_cases STUDYremovedMAFS_controls
#Trexoume to tool me auta ta arxeia ws input kai flag -HWE wste na upologistei to pvalue gia kathe thesi
#Output: STUDY.hwe

with open('STUDY.hwe')as file:
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
            snp_number= int(splittedline[0].split('_')[1]) #krataw mono ton arithmo tou snp (kai oxi to string px snp_142)    
            pvalue= round(float(splittedline[1]), 3)       #stroggulopoiw gia na min pethanw to plot        
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
ax.set_xlabel("SNP_ID", fontsize=15)
ax.set_ylabel("p-value", fontsize=15)
ax.set_title("Hardy-Weinberg p-values Distribution", fontsize= 25) 
plt.savefig('hwe_plot.jpg')
plt.show()       




      
#%%               REMOVAL
with open('Path to --> STUDYremovedMAFS_cases') as cases, open('Path to --> STUDYremovedHWE_cases', 'w') as output_cases:
        for line_cases in cases:
            if line_cases.split(' ')[0] not in remove:
                output_cases.write(line_cases)
with open('Path to --> STUDYremovedMAFS_controls') as controls, open('Path to --> STUDYremovedHWE_controls', 'w') as output_controls:
        for line_controls in controls:
            if line_controls.split(' ')[0] not in remove:
                output_controls.write(line_controls)