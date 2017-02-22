#treksimo
#python allele_freqNEW.py -controls_file data/gwas.controls.gen -cases_file data/gwas.cases.gen -output outputallelefreq -allele_frequency
###########allakse onoma script

#The NCBI "Reference allele" for a given SNP refers to the nucleotide base on the NCBI reference assembly at the SNP’s position
#Minor allele frequency (MAF) refers to the frequency at which the second most common allele occurs in a given population.

import time
start = time.time() 
with open('/home/rantaplan/master/projectara/output/outputallelefreq.frequency') as output:
    #gia to plot
    X=[]
    controls=[]
    cases=[]
    
    #gia to removal twn snp me maf<0,5
    cases_remove=[]
    controls_remove=[]
    
    for line in output:
        line=line.rstrip('\n')        
        splittedline= line.split(' ')
        
        #apothikeuse stin lista twn X to snp
        snp_number= int(splittedline[0].split('_')[1]) #krataw mono ton arithmo tou snp (kai oxi to string snp_142)    
        X.append(snp_number)
        
        controls_ref=float(splittedline[1])
        controls_alt=float(splittedline[2])
        cases_ref=float(splittedline[3])
        cases_alt=float(splittedline[4])
        
        #vres poioi einai to MAF(anamesa se ref kai alt) kai prosthese to stin lista twn Y
        if controls_alt < controls_ref:
            controls_maf = controls_alt
        else:
            controls_maf = controls_ref
        if cases_alt < cases_ref:
            cases_maf = cases_alt
        else:
            cases_maf = cases_ref
        
        cases.append(cases_maf)
        controls.append(controls_maf)
        
        #vres ta snp me maf<0,5
        if cases_maf < 0.05:
            cases_remove.append(splittedline[0])    
        if controls_maf < 0.05:
            controls_remove.append(splittedline[0])            
        
print(time.time()-start)

#%%                     PLOT
import matplotlib.pyplot as plt
fig, ax = plt.subplots()         

#megalwma diastasewn  odigies: https://codeyarns.com/2014/10/27/how-to-change-size-of-matplotlib-plot/  
fig_size = plt.rcParams["figure.figsize"]  #->fig size; 6, 4
fig_size[0] = 10
fig_size[1] = 10
plt.rcParams["figure.figsize"] = fig_size        
        
#name labels            
ax.set_xlabel("SNP ID")
ax.set_ylabel("MAF")
ax.set_title("Minor Allele Frequency Distribution")

ax.plot(X, cases , '.')
ax.plot(X, controls , '.', c="peru")

#=====================APOTIXIMENI PROSPA8EIA NA VALW UPOMNIMA=========================================================
# legends=ax.plot(X, cases , '.', 'b', X , controls,'.', 'r' ) 
# plt.legend(legends, ["cases", "controls"], loc=2)
#==============================================================================
plt.savefig('mafs.jpg')
plt.show()
#%%                 SNP REMOVAL
#   Epeidi thelw na afairesw diaforetika snp apo cases kai diaforetika apo controls (den mporw na to treksw dld apo command line)
import time
start_removal = time.time() 
with open('/home/rantaplan/master/projectara/data/gwas.cases.gen') as cases, open('maf_removal_cases', 'w') as output_cases:
        for line_cases in cases:
            if line_cases.split(' ')[0] not in  cases_remove:
                output_cases.write(line_cases)
with open('/home/rantaplan/master/projectara/data/gwas.controls.gen') as controls, open('maf_removal_controls', 'w') as output_controls:
        for line_controls in controls:
            if line_controls.split(' ')[0] in controls_remove:
                output_controls.write(line_controls)
                 
print(time.time()-start_removal)


#xronos 613.176543712616 (10lepta)
#remove apo controls: 47689
# apo cases 47450
















