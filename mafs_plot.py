#Exoume treksei to tool sto arxiko dataset gwas.cases.gen kai gwas.controls.gen
#COMMAND: python PROJECTARA.py -controls_file data/gwas.controls.gen -cases_file data/gwas.cases.gen -output STUDY -allele_frequency
start = time.time() 
with open('STUDY.frequency') as output:
    #gia to plot
    X=[]
    controls=[]
    cases=[]
    
    #gia to removal twn snp me maf<0,5
    remove=[]
       
    for line in output:
        line=line.rstrip('\n')        
        splittedline= line.split(' ')
        
        #apothikeuei stin lista twn X to snp
        snp_number= int(splittedline[0].split('_')[1]) #krataw mono ton arithmo tou snp (kai oxi to string snp_142)    
        X.append(snp_number)
        
        controls_ref=float(splittedline[1])
        controls_alt=float(splittedline[2])
        cases_ref=float(splittedline[3])
        cases_alt=float(splittedline[4])
        merged_ref=float(splittedline[5])
        merged_alt=float(splittedline[6])
        #vriskei poioi einai to MAF(anamesa se ref kai alt) 
        if controls_alt < controls_ref:
            controls_maf = controls_alt
        else:
            controls_maf = controls_ref
        if cases_alt < cases_ref:
            cases_maf = cases_alt
        else:
            cases_maf = cases_ref
        
        if merged_alt < merged_ref:
            merged_maf = merged_alt
        else:
            merged_maf = merged_ref
        if merged_maf < 0.05:
            remove.append(splittedline[0])
            
        cases.append(cases_maf)
        controls.append(controls_maf)
        
        

#%%                     PLOT
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

fig_size = plt.rcParams["figure.figsize"]  # default diastaseis askonwn: 4, 6
fig_size[0] = 150 # dieurunsh aksona x
fig_size[1] = 10 # aksonas y
plt.rcParams["figure.figsize"] = fig_size 

Cases, = ax.plot(X, cases , '.')
Controls, = ax.plot(X, controls , '.', c="magenta")  

plt.legend([Cases,Controls], ["Cases", "Controls"], loc=1) # loc=2 shmainei panw aristera 
ax.set_xlabel("SNP ID", fontsize=15)
ax.set_ylabel("Minor Allele Frequencies", fontsize=15)
ax.set_title("Minor Allele Frequency Distribution", fontsize= 25)       


plt.savefig('mafs.jpg')
plt.show()     
        



#%%                 SNP REMOVAL

with open("path to --> 'gwas.cases.gen'") as cases, open("path to --> 'STUDYremovedMAFS_cases'", 'w') as output_cases:
        for line_cases in cases:
            if line_cases.split(' ')[0] not in remove:
                output_cases.write(line_cases)
with open("path to --> 'gwas.controls.gen'") as controls, open("path to -->'STUDYremovedMAFS_controls'", 'w') as output_controls:
        for line_controls in controls:
            if line_controls.split(' ')[0] not in remove:
                output_controls.write(line_controls)


















