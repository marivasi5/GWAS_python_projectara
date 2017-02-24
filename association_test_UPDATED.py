#GIA TREKSIMO
#python allele_freqNEW.py -controls_file data/controlsmikraki.txt -cases_file data/casesmikraki.txt -output testaki -allele_frequency
import os
import sys
import argparse

#einai apo to arxeio argtest. einai mono oi grammes pou xreiazontai gia na treksw to allelefreq
parser = argparse.ArgumentParser(description='This is our projectara')
parser.add_argument('-controls_file', help='File containing control genotypes',required=True)
parser.add_argument('-cases_file', help='File containing cases genotypes',required=True)
parser.add_argument('-output', help='Output file name',required=True)
parser.add_argument('-allele_frequency', action='store_true')           #!ftiakse ta help
#NEW FLAG!!!
parser.add_argument('-HWE', action='store_true')
parser.add_argument('-association_test', action='store_true')
args = parser.parse_args()

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
    return snp , homozygous_refrence, homozygous_alternative, heterozygous, N , locus
            
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
def HWE(x,y):                               # x=counts_cases, y=counts_controls
    '''Ypologismos tou Hardy Weinberg Equilibrium statistic gia to enwmeno cases+controls dataset
    *Prepei prwta na treksi i genotype_counts gia ta dataset ksexwrista*
    Input: 2 tuple tis morfis (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn, thesi SNP sumfwna me to NCBI build 36)
    X: twn CASES    Y: twn contols
    Output: tuple(snp_ID, pvalue)'''       
    from scipy import stats #This test is invalid when the observed or expected frequencies in each category are too small. A typical rule is that all of the observed and expected frequencies should be at least 5.       
    
    counts_merged= (y[0], x[1]+y[1], x[2]+y[2], x[3]+y[3], x[4]+y[4], x[5]) 
    
    if counts_merged[1]!=0 and counts_merged[2]!=0  and counts_merged[3]!=0: 
        snp, p_merged, q_merged = allele_freq(counts_merged)
        N = counts_merged[4]
        expect_homozygous_refrence = (p_merged**2)*N       #pairnoume to plh8os twn anamenomenwn atomwn onozugwn ws pros reference
        expect_homozygous_alternative = (q_merged**2)*N
        expect_heterozygous = p_merged*q_merged*2*N
        x2test_hw = stats.chisquare([counts_merged[1], counts_merged[3],counts_merged[2]], [expect_homozygous_refrence, expect_heterozygous, expect_homozygous_alternative])
    
        return snp, x2test_hw.pvalue
    else:
        snp = counts_merged[0]
        return snp, "cannot compute pvalue"

#%%               #Genotypic Association --> http://www.gwaspi.org/?page_id=332
def association_test(x,y):                   #x=counts_cases, y=counts_controls
    '''Pragmatopoiei Genotypic Association test
    *Prepei prwta na treksi i genotype_counts gia ta dataset ksexwrista*
    Input: 2 tuple tis morfis (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn, thesi SNP sumfwna me to NCBI build 36)
    X: twn CASES    Y: twn contols
    Output: tuple (snp_ID, thesi SNP sumfwna me to NCBI build 36, pvalue) tou test, odds ratio, odds ratio'''
    from scipy import stats            
    snp = x[0]
    loci = x[5]
    
    if x[1]!=0 and y[1]!=0 and x[2]!=0 and y[2]!=0 and x[3]!=0 and y[3]!=0:
        snp, x2test_hw= HWE(x,y)
        snp, p_cases, q_cases= allele_freq(x)
        snp, p_controls, q_controls= allele_freq(y)
        N=x[4]
        x2test_ga = stats.chisquare([x[1],y[1],x[2],y[2],x[3],y[3]], [(p_cases**2)*N,(p_controls**2)*N, (q_cases**2)*N, (q_controls**2)*N ,2*p_cases*q_cases*N, 2*p_controls*q_controls*N]) 
        x2test_ga_final= x2test_ga.pvalue
    else:
        x2test_ga_final = "cannot compute pvalue"
    
    if x[2]!= 0 and y[1]!=0 and y[3]!=0:
        OR_RRAA = round((x[1]*y[2])/(x[2]*y[1]), 5)
        OR_RAAA = round((x[3]*y[2])/(x[2]*y[3]), 5)
    else:
        OR_RRAA = "none"
        OR_RAAA = "none"

    return snp, loci, x2test_ga_final, OR_RRAA, OR_RAAA

#%%                 ASSOCIATION TREKSIMO
#==============================================================================
# #MARIA: import os ; os.chdir('/home/rantaplan/master/projectara/data/')
# with open('gwas.cases.gen') as cases, open('gwas.controls.gen') as controls:
#     j=0
#     #gia to manhattan plot
#     loci_list = []
#     pvalue_list = []
#     
#     for line_cases, line_controls in zip(cases, controls):
#         line_cases=line_cases.rstrip('\n')
#         line_controls=line_controls.rstrip('\n')
#         if j % 2000 == 0:
#             print(j)   
#         j += 1
#                 
#         counts_cases=genotype_counts(line_cases)                
#         counts_controls=genotype_counts(line_controls)
#                 
#         snp, loci, pvalue, OR_RRAA, OR_RAAA  = association_test(counts_cases,counts_controls)
#         loci_list.append(loci)
#         pvalue_list.append(pvalue)
#         print(snp, loci, pvalue, OR_RRAA, OR_RAAA)
# 
#==============================================================================
#%%
if args.association_test:
    import time
    start = time.time()     
    j = 0
    #gia manhattan plot
    loci_list = []
    pvalue_list = []
    
    output= open('{}.association'.format(args.output), 'w')
    with open(args.cases_file) as cases, open(args.controls_file) as controls:
        for line_cases, line_controls in zip(cases, controls):
            line_cases=line_cases.rstrip('\n')
            line_controls=line_controls.rstrip('\n')
            if j % 2000 == 0:
                print(j)   
            j += 1
                
            counts_cases=genotype_counts(line_cases)                
            counts_controls=genotype_counts(line_controls)
                
            snp, loci, pvalue, OR_RRAA, OR_RAAA  = association_test(counts_cases,counts_controls)
            loci_list.append(loci)
            pvalue_list.append(pvalue)
            print(snp, loci, pvalue, OR_RRAA, OR_RAAA, file=output)
                    
    print(time.time()-start)     
#%%
                
















