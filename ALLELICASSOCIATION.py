import os

os.chdir('gwas')

#%%
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
    return snp ,homozygous_refrence, homozygous_alternative, heterozygous, N, locus
            
#%%                
def allele_freq(x):   
    '''Ypologismos suxnotitas allilomorfwn (p=suxnotita refrence allilomorfou, q=suxnotita alternative allilomorfou)
    *Prepei prwta na treksi i genotype_counts*
    Input: tuple tis morfis (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn, thesi SNP sumfwna me to NCBI build 36)
    Output: tuple (snp_ID, p, q) '''
    
    snp, R, A, het, N, loci = x
    p= round((R*2 + het)/(2*N), 3)  
    q= round((A*2 + het)/(2*N), 3)
    return snp, p, q
    

 #%%            
def major_minor(x):  # controls tuple genotype_counts/allele_freq twn controls, cases genotype_counts/allele_freq twn cases
    '''Ypologismos tou plh8ous twn major kai minor alleles se ka8e snp
    *Prepei prwta na treksei i genotype_counts*
    Input: tuple tis morfis(snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn, thesi SNP sumfwna me to NCBI build 36)
    Output: tuple (snp_ID, ari8mos atomwn gia to major allilomorfo, ari8mos atomwn gia to minor allilomorf, ari8mos eterozugwn atomwn,sunolo atomwn, thesi SNP sumfwna me to NCBI build 36) '''
    snp, R, A, het, N, locus = x
         
    if R < A:  
        minor = R
        major = A
    else:
        minor = A   
        major = R
    
                
    return snp, major, minor, het, N, locus

#%%

def allelic_association_test(a,b): 
    '''Ypologismos tou allelic association test gia kathe SNP se controls kai cases.
    *Prepei prwta na treksoun i genotype_counts kai i major_minor*
    Input: 2 tuple tis morfis (snp_ID, arithmos  atomwn gia to major allilomorfo sta controls/cases, arithmos atomwn gia to minor allilomorfo sta controls/cases, arithmos eterozugwn atomwn sta controls/cases, sunolo atomwn, thesi SNP sumfwna me to NCBI build 36)
    a: twn CONTROLS   b: twn CASES
    Output: tuple (snp_ID, pvalue, odds ratio) '''
    
    from scipy import stats 
    snp, major_controls_c, minor_controls_c, het_controls, N, locus = a
    snp, major_cases_c, minor_cases_c, het_cases, N, locus  = b
    
    counts_major_controls = 2*major_controls_c+het_controls
    counts_minor_controls = 2*minor_controls_c+het_controls
    counts_major_cases = 2*major_cases_c+het_cases
    counts_minor_cases = 2*minor_cases_c+het_cases
    
    sum_major = counts_major_controls + counts_major_cases
    sum_minor = counts_minor_controls + counts_minor_cases
    
    exp_major_cases= (sum_major * N)/(2*N)
    exp_minor_cases=  (sum_minor * N)/(2*N)
    exp_major_controls= (sum_major * N)/(2*N)
    exp_minor_controls = (sum_minor * N)/(2*N)
    #x2test_ga = stats.chisquare([x[1]*z[4], x[2]*z[4], x[3]*z[4], x[4]*z[4]], []
    
                                 
    x2test_aa = stats.chisquare([counts_major_controls, counts_minor_controls ,counts_major_cases , counts_minor_cases], [exp_major_controls, exp_minor_controls, exp_major_cases, exp_minor_cases])
    #2*major_controls_c+het_controls = o ari8mos twn observed atomwn me to major allele
    #2*N*((major_controls_f**2)+2*major_controls_f*minor_controls_f) = expected count tou major allele
    if minor_cases_c != 0 and major_controls_c != 0:
        OR= (major_cases_c * minor_controls_c) / (minor_cases_c * major_controls_c)
    else:
        OR = "Nan"
        
    return snp, locus, x2test_aa.pvalue,OR
    
    
#==============================================================================
# If OR is < 1:
# OR = 1/OR
# 
#     
#==============================================================================
#%%
#Treksimo Allelic assosciation

j=0
with open('gwas.cases.gen') as cases, open('gwas.controls.gen') as controls:            
    for line_cases, line_controls in zip(cases, controls):
             line_cases=line_cases.rstrip('\n')
             line_controls=line_controls.rstrip('\n')
             if j % 2000 == 0:
                 print(j)   
             j += 1
             
             counts_cases = genotype_counts(line_cases)                
             counts_controls = genotype_counts(line_controls)
             
             a = major_minor(counts_cases)#major-minor counts gia cases
             b = major_minor(counts_controls)#major-minor counts gia controls
            
             snp, locus, pvalue, OR = allelic_association_test(a,b)
             print(snp, locus, pvalue, OR)
