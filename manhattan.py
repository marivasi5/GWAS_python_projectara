#%%
import os

os.chdir('gwas')
#%%
#%%
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
parser.add_argument('-association_test', action='store_value')
args = parser.parse_args()




#%%
def genotype_counts(datasetLINE):        
    '''Ypologismos arithmou gonotupwn
    Input: grammi arxeiou se Genotype File Format 
    Output: tuple tis morfis (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn, thesi SNP)'''
    
    splittedline= datasetLINE.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
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
def allele_freq(x):  #isws prepei na to treksoume etsi gia na vroume kai gia to merged xwris na exoume kanei merging 
    '''Ypologismos suxnotitas allilomorfwn (p=suxnotita refrence allilomorfou, q=suxnotita alternative allilomorfou)
    *Prepei prwta na treksi i genotype_counts*
    Input: tuple  (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn)
    Output: tuple (snp_ID, p, q) '''
    
    snp, R, A, het, N, loci = x
    p= round((R*2 + het)/(2*N), 3)  
    q= round((A*2 + het)/(2*N), 3)
    return snp, p, q
    
#%%


def HWE(x,y):# x=counts_cases, y=counts_controls
    from scipy import stats
    
    counts_merged= (y[0], x[1]+y[1], x[2]+y[2], x[3]+y[3], x[4]+y[4], x[5]) 
    
    if counts_merged[1] > 0 and counts_merged[2] > 0  and counts_merged[3] > 0: #This test is invalid when the observed or expected frequencies in each category are too small. A typical rule is that all of the observed and expected frequencies should be at least 5.
        snp, p_merged, q_merged = allele_freq(counts_merged)
        N = counts_merged[4]
        expect_homozygous_refrence = (p_merged**2)*N #pairnoume to plh8os twn anamenomenw reference
        expect_homozygous_alternative = (q_merged**2)*N
        expect_heterozygous = p_merged*q_merged*2*N
        x2test_hw = stats.chisquare([counts_merged[1], counts_merged[3],counts_merged[2]], [expect_homozygous_refrence, expect_heterozygous, expect_homozygous_alternative])
    
        return snp, x2test_hw.pvalue, round(expect_homozygous_refrence,3), round(expect_homozygous_alternative,3), round(expect_heterozygous,3)
    else:
        snp = counts_merged[0]
        return snp, "cannot compute pvalue", "ehr", "eha", "ehet"
#%%
#%%
#Genotypic Association --> http://www.gwaspi.org/?page_id=332, https://en.wikipedia.org/wiki/Genome-wide_association_study


'''After odds ratios and P-values have been calculated for all SNPs, a common approach is to create a 
Manhattan plot. In the context of GWA studies, this plot shows the negative logarithm of the P-value
as a function of genomic location. Thus the SNPs with the most significant association stand out on
the plot, usually as stacks of points because of haploblock structure. Importantly, the P-value 
threshold for significance is corrected for multiple testing issues. The exact threshold varies by 
study,(Wittkowski et al. 2014), but the conventional threshold is 5×10−8 to be significant in the face of hundreds of 
thousands to millions of tested SNPs (Bush & Moore 2012, Clarke et al. 2011,Barsh et al. 2012). 
GWA studies typically perform the first analysis in a discovery cohort, followed by validation of 
the most significant SNPs in an independent validation cohort.'''

def association_test(x,y): # taizw tuples tou genotype_counts # x=counts_cases, y=counts_controls
            from scipy import stats            
            snp = x[0]
            loci = x[5]
            if x[1] > 0 and y[1] > 0 and x[2] >0 and y[2] > 0 and x[3] > 0 and y[3] > 0:
                snp,x2test_hw,ehr,eha,ehet = HWE(x,y)
                snp, p_cases, q_cases= allele_freq(x)
                snp, p_controls, q_controls= allele_freq(y)
                N=x[4]
                x2test_ga = stats.chisquare([x[1],y[1],x[2],y[2],x[3],y[3]], [(p_cases**2)*N,(p_controls**2)*N, (q_cases**2)*N, (q_controls**2)*N ,2*p_cases*q_cases*N, 2*p_controls*q_controls*N]) 
                x2test_ga_final= x2test_ga.pvalue
            else:
                x2test_ga_final = "cannot compute pvalue"
            if x[2]!= 0 and y[1]!=0 and y[3]!=0:
    
                OR_RRAA = (x[1]*y[2])/(x[2]*y[1])
                OR_RAAA = (x[3]*y[2])/(x[2]*y[3])
            else:
                OR_RRAA = "none"
                OR_RAAA = "none"

            return snp, loci, x2test_ga_final, OR_RRAA, OR_RAAA
            
#%%             

    

#%%                  HWE treksimo

with open('gwas.cases.gen') as cases, open('gwas.controls.gen') as controls:
    j=0
    
    for line_cases, line_controls in zip(cases, controls):
        
                line_cases=line_cases.rstrip('\n')
                line_controls=line_controls.rstrip('\n')
                #gia elegxo taxutitas, mou deixnei poses fores exei treksei(ara poses grammes)
                
                if j % 2000 == 0:
                     print(j)   
                j += 1
                
                counts_cases=genotype_counts(line_cases)                
                counts_controls=genotype_counts(line_controls)
                
                snp, x2test_hw,ehr,eha,ehet = HWE(counts_cases,counts_controls)
                print(snp, x2test_hw)

#%%                 ASSOCIATION TREKSIMO
with open('gwas.cases.gen') as cases, open('gwas.controls.gen') as controls:
    j=0
    loci_list = []
    pvalue_list = []
    
    for line_cases, line_controls in zip(cases, controls):
        
                line_cases=line_cases.rstrip('\n')
                line_controls=line_controls.rstrip('\n')
                #gia elegxo taxutitas, mou deixnei poses fores exei treksei(ara poses grammes)
                
                if j % 2000 == 0:
                     print(j)   
                j += 1
                
                counts_cases=genotype_counts(line_cases)                
                counts_controls=genotype_counts(line_controls)
                
                snp, locus, pvalue, OR  = allelic_association_test(counts_cases,counts_controls)
                loci_list.append(locus)
                pvalue_list.append(pvalue)
                #print(snp, loci, pvalue, OR_RRAA, OR_RAAA)






#%% 
#==============================================================================
#Manhattan plot
def manhattan(x,y): #x einai to pvalue_list kai y to loci_list p exoume parei ap to treksimo ths allele_association_test
    import matplotlib.pyplot as plt
    import math
    
    pvalues = list(map(lambda x:(-math.log(x)),pvalue_list))#efarmozw ton arnhtiko logari8mo sola ta pvalues
    plt.plot(loci_list, pvalues, ls='', marker='.', color='green')
    plt.xlim([0,int(loci_list[-1])])#oria aksona x,mia 8esh akrivws meta to telutaio snp
    plt.title("Association Manhattan Plot", fontsize=15)
    plt.xlabel("Position on chr 20", fontsize=10)
    plt.ylabel("-log(Pvalue)", fontsize=10)
    plt.savefig("manhattan")
    plt.show()

#%%
#qq plot --> https://en.wikipedia.org/wiki/Q%E2%80%93Q_plot
'''In statistics, a Q–Q plot[1] ("Q" stands for quantile) is a probability plot, which is a 
graphical method for comparing two probability distributions by plotting their quantiles against 
each other.A Q–Q plot is used to compare the shapes of distributions, providing a graphical view of how
properties such as location, scale, and skewness are similar or different in the two distributions.
Q–Q plots can be used to compare collections of data, or theoretical distributions.The use of Q–Q 
plots to compare two samples of data can be viewed as a non-parametric approach to comparing their
 underlying distributions. A Q–Q plot is generally a more powerful approach to do this than the
common technique of comparing histograms of the two samples, but requires more skill to interpret. 
Q–Q plots are commonly used to compare a data set to a theoretical model( Gnanadesikan 1977,
Thode 2002). This can provide an assessment of "goodness of fit" that is graphical, rather than 
reducing to a numerical summary. Q–Q plots are also used to compare two theoretical distributions 
to each other((Gibbons & Chakraborti 2003).Since Q–Q plots compare distributions, there is no need
 for the values to be observed as pairs, as in a scatter plot, or even for the numbers of values 
 in the two groups being compared to be equal.'''


def qqplot(x): #x einai to p value_list opws kai sth manhattan
    from scipy import  stats
    import math
    import pylab 
    pvalues = list(map(lambda x:(-math.log(x)),pvalue_list))#efarmozw ton arnhtiko logari8mo sola ta pvalues
    stats.probplot(pvalues, dist = stats.exponnorm,sparams=(2.5,), plot=pylab)
    pylab.set_title("Probplot for exponential distr with shape parameter 2.5")
    pylab.savefig("qqplot")
    pylab.show()
#%%
#OR  #POLY PIO GRHGOROS!!
def qqplot(x): #x einai to p value_list opws kai sth manhattan
    import statsmodels.api as sm
    import pylab
    import numpy as np
    import math
    pvalues = list(map(lambda x:(-math.log(x)),pvalue_list))#efarmozw ton arnhtiko logari8mo sola ta pvalues
    pv = np.asarray(pvalues)
    sm.qqplot(pv, line='s')# etoimo module gia qqplot
    pylab.set_title("Probplot for exponential distr with shape parameter 2.5")
    pylab.savefig("qqlot2")
    pylab.show()


        



