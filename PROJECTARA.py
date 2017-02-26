
import argparse

parser = argparse.ArgumentParser(description='This is our projectara')
parser.add_argument('-controls_file', help='File containing control genotypes in Genotype File Format.',required=True)
parser.add_argument('-cases_file', help='File containing cases genotypes in Genotype File Format.',required=True)
parser.add_argument('-output', help='Output file name.',required=True)
parser.add_argument('-keep_snps', help= 'File containing snp_IDs per line. Those will be preserved in the dataset. The IDs must be in snp_NUMBER format.')
parser.add_argument('-remove_snps', help= 'File containing snp_IDs per line. Those will be excluded from the dataset. The IDs must be in snp_NUMBER format.')
parser.add_argument('-keep_samples', help= 'File containing sample IDs per line. Those will be preserved in the dataset. The IDs must be in either cases_NUMBER or controls_NUMBER format.')
parser.add_argument('-remove_samples' , help= 'File containing sample IDs per line. Those will be excluded from the dataset. The IDs must be in either cases_NUMBER or controls_NUMBER format')
parser.add_argument('-allele_frequency', action='store_true', help='Returns one-line-per-SNP information: snp_ID, refrence allele frequency in controls, alternative allele frequency in controls, refrence allele frequency in cases, alternative allele frequency in cases, refrence allele frequency in both controls and cases, alternative allele frequency in both controls and cases')           
parser.add_argument('-HWE', action='store_true', help='Returns one-line-per-SNP information: snp_ID, Hardy Weinberg Equilibrium statistic') 
parser.add_argument('-get_info', help='snp_ID to request from Ensembl Variant Predictor. The IDs must be in snp_NUMBER format.')                   

args = parser.parse_args()

#%%
def createSNPlist(file):
    '''Input:File containing snp_IDs per line
    Output: Lista me ta snp_IDs'''
    
    with open(file) as snpfile:
        snplist=[]
        for line in snpfile:      
            if line[0:4]=='snp_':
                line=line.rstrip('\n') 
                snp= line.split(' ')[0] 
                snplist.append(snp)
    snplist.sort()     
    return snplist
#%%
def get_samples(file):
    '''Input: File containing sample IDs per line
    Output: Mia lista me tous kwdikous pou anaferontai se deigmata cases kai mia me kwdikous pou aforoun controls'''
    
    with open(file) as samplelist:
        case_samples=[]    
        control_samples=[]
        for line in samplelist:
            line=line.rstrip('\n')
            if line.split('_')[0] == 'case':
                case_samples.append(int(line.split('_')[1]))   #int giati to diavazei san string 
            elif line.split('_')[0] == 'control':
                control_samples.append(int(line.split('_')[1]))    
        return case_samples, control_samples

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
    
#%%                         KEEP SNPS
if args.keep_snps is not None: 
    snplist=createSNPlist(args.keep_snps) 
    if len(snplist)==0:
        print('Please provide file containing snp_IDs per line in the required format. Try "-help" for more information.')
    
    else:
        with open(args.cases_file) as cases, open('{}.cases.gen'.format(args.output), 'w') as output_cases:
            for line_cases in cases:
                if line_cases.split(' ')[0] in snplist:
                    output_cases.write(line_cases)
                    
        with open(args.controls_file) as controls, open('{}.controls.gen'.format(args.output), 'w') as output_controls:
            for line_controls in controls:
                if line_controls.split(' ')[0] in snplist:
                    output_controls.write(line_controls)
                
#%%                        REMOVE SNPS
if args.remove_snps is not None:
    snplist=createSNPlist(args.remove_snps) 
    if len(snplist)==0:
        print('Please provide file containing snp_IDs per line in the required format. Try "--help" for more information.')
    
    else:
        with open(args.cases_file) as cases, open('{}.cases.gen'.format(args.output), 'w') as output_cases:
            for line_cases in cases:
                if line_cases.split(' ')[0] not in snplist:
                    output_cases.write(line_cases)
        with open(args.controls_file) as controls, open('{}.controls.gen'.format(args.output), 'w') as output_controls:
            for line_controls in controls:
                if line_controls.split(' ')[0] in snplist:
                    output_controls.write(line_controls)
    
#%%                         KEEP SAMPLES
if args.keep_samples is not None:    
    control_samples= sorted(get_samples(args.keep_samples)[1])
    case_samples= sorted(get_samples(args.keep_samples)[0])
    if len(case_samples)==0 and len(control_samples)==0:
         print('Please provide file containing sample IDs per line in the required format. Try "--help" for more information.')
    else: 
        with open(args.controls_file) as controls, open('{}.controls.gen'.format(args.output), 'w') as output_controls:
            for line_controls in controls:
                line_controls=line_controls.rstrip('\n')
                splittedline_controls= line_controls.split(' ')   
                newline_controls= splittedline_controls[0:5]  
                for i in control_samples:
                    index=(i-1)*3+5
                    individual= splittedline_controls[index]+ splittedline_controls[index+1]+splittedline_controls[index+2]
                    newline_controls+=individual
                newline_controls_string=' '.join(newline_controls) + '\n'
                output_controls.write(newline_controls_string) 
        
        with open(args.cases_file) as cases, open('{}.cases.gen'.format(args.output), 'w') as output_cases:
            for line_cases in cases:
                line_cases=line_cases.rstrip('\n')
                splittedline_cases= line_cases.split(' ')   
                newline_cases= splittedline_cases[0:5]  
                for i in case_samples:
                    index=(i-1)*3+5
                    individual= splittedline_cases[index]+ splittedline_cases[index+1]+splittedline_cases[index+2]
                    newline_cases+=individual
                newline_cases_string=' '.join(newline_cases) + '\n'
                output_cases.write(newline_cases_string) 

#%%                       REMOVE SAMPLES
if args.remove_samples is not None:    
    control_samples= sorted(get_samples(args.remove_samples)[1], reverse=True)
    case_samples= sorted(get_samples(args.remove_samples)[0], reverse=True)
    if len(case_samples)==0 and len(control_samples)==0:
         print('Please provide file containing sample IDs per line in the required format. Try "--help" for more information.')
    else: 
        with open(args.controls_file) as controls, open('{}.controls.gen'.format(args.output), 'w') as output_controls:
            for line_controls in controls:
                line_controls=line_controls.rstrip('\n')
                splittedline_controls= line_controls.split(' ')   
                for i in control_samples:
                    index=(i-1)*3+5
                    del splittedline_controls[index:index+3]
                newline_controls_string=' '.join(splittedline_controls) + '\n'
                output_controls.write(newline_controls_string) 
        
        with open(args.cases_file) as cases, open('{}.cases.gen'.format(args.output), 'w') as output_cases:
            for line_cases in cases:
                line_cases=line_cases.rstrip('\n')
                splittedline_cases= line_cases.split(' ')  
                for i in case_samples:
                    index=(i-1)*3+5
                    del splittedline_cases[index:index+3] 
                newline_cases_string=' '.join(splittedline_cases) + '\n'
                output_cases.write(newline_cases_string) 

#%%                     ALLELE FREQUENCY 
if args.allele_frequency:
    with open(args.cases_file) as cases, open(args.controls_file) as controls,open('{}.frequency'.format(args.output), 'w') as output :
        for line_cases,line_controls in zip(cases, controls):
            line_cases=line_cases.rstrip('\n')
            line_controls=line_controls.rstrip('\n')        
            
            counts_cases=genotype_counts(line_cases)                
            counts_controls=genotype_counts(line_controls)
            counts_merged= (counts_controls[0], counts_cases[1]+counts_controls[1], counts_cases[2]+counts_controls[2], counts_cases[3]+counts_controls[3], counts_cases[4]+counts_controls[4], counts_controls[5]) 
            
            snp, p_controls, q_controls = allele_freq(counts_controls)
            snp, p_cases, q_cases = allele_freq(counts_cases)
            snp, p_merged, q_merged = allele_freq(counts_merged)
            
            print(snp, p_controls, q_controls, p_cases, q_cases, p_merged, q_merged, file=output)     

#%%                             HWE
if args.HWE:
    
    def HWE(x,y):                               # x=counts_cases, y=counts_controls
    '''Ypologismos tou Hardy Weinberg Equilibrium statistic gia to enwmeno cases+controls dataset
    *Prepei prwta na treksi i genotype_counts gia ta dataset ksexwrista*
    Input: 2 tuple tis morfis (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn, thesi SNP sumfwna me to NCBI build 36)
    X: twn Cases    Y: twn contols
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
        
    with open(args.cases_file) as cases, open(args.controls_file) as controls, open('{}.hwe'.format(args.output), 'w') as output:
        for line_cases,line_controls in zip(cases, controls):
            line_cases=line_cases.rstrip('\n')
            line_controls=line_controls.rstrip('\n')
                
            counts_cases=genotype_counts(line_cases)                
            counts_controls=genotype_counts(line_controls)
            
            snp, x2test_hw = HWE(counts_cases,counts_controls)
            print(snp, x2test_hw, file=output)

#%%















#%%                             GET INFO
if args.get_info is not None:
    snp= args.get_info
    if snp[:4]=='snp_':
        import requests
        from pyliftover import LiftOver     #gia to liftover
        lo=LiftOver('hg18', 'hg38')         # downloads the hg18-to-hg38 coordinate conversion chain file from UCSC
        
        #get position
        with open(args.cases_file) as cases:
            for line in cases:
                splittedline= line.split(' ')
                if splittedline[0] == snp:
                    position = int(splittedline[2])
                    refrence= splittedline[3]
                    alternative= splittedline[4]
        #LIFTOVER    
        conversion= lo.convert_coordinate('chr20', position)
        #The result is either None (if the source chromosome name is unrecognized) or a list of target positions in the new assembly. The list may be empty (locus is deleted in the new assembly), have a single element (locus matched uniquely)
        try:
            if len(conversion)==1:
                position_new=conversion[0][1]
                #GET REQUEST
                url= 'http://rest.ensembl.org/vep/human/hgvs/20:g.' + str(position_new) + refrence + '>' + alternative + '?' #to xromoswma twra einai to 20 alla tha eprepe na to vazw
                headers={ "Content-Type" : "application/json"}
                r = requests.get(url, headers=headers)
                            
                if r.ok:
                    data = r.json()
                    with open('{}.info'.format(args.output), 'w') as output:
                        output.write('The information acquired from Ensembl\'s Variant Effect Predictor for {} is: \n {}'.format(snp, data))    
        
                else:    
                    print('The request for {} was not successful.'.format(snp))    
            else:
                print('The request for {} was not successful. Could not convert position to USCS hg38.'.format(snp))
        
        except TypeError:       #object of type 'NoneType' has no len()
            print('The request for {} was not successful. Could not convert position to USCS hg38.'.format(snp))
    else:
        print('Please provide one snp_ID in the required format. Try "--help" for more information.')
   






