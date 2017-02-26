import sys
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
parser.add_argument('-HWE', action='store_true', help='Returns one-line-per-SNP information: snp_ID, Hardy Weinberg Equilibrium statistic.')
parser.add_argument('-get_info', help='snp_ID to request from Ensembl Variant Predictor. The IDs must be in snp_NUMBER format. REQUIREMENTS: Please install requests (http://docs.python-requests.org/en/master/user/install/) AND pyLiftover from "https://github.com/konstantint/pyliftover".')                   
parser.add_argument('-association_test', action='store_true', help= 'Performs allelic association test between cases and controls. Returns one-line-per-SNP information: snp_ID, locus, p value calculated from the test, odds ratio between major and minor allele in cases and controls.')
parser.add_argument('-manhattan', action='store_true', help= 'Requires the "-association_test" beforehand. Creates an manhattan plot of the pvalues calculated from the assosiation test.')
parser.add_argument('-qqplot', action='store_true', help= 'Requires the "-association_test" beforehand. Creates an qqplot of the pvalues calculated from the assosiation test.')
parser.add_argument('-LD', nargs=2, help='Requires two snp_IDs. Checks the linkage disequilibrium in the dataset of cases. Returns one line information: snp_IDs, rsquare, Dprime. The IDs must be in snp_NUMBER format.')
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

#%%                             GET INFO

if args.get_info is not None:
    snp= args.get_info
    if snp[:4]=='snp_':
        try:
            import requests
        except ModuleNotFoundError:
            print('Please install requests Try "--help" for more information.')            
            sys.exit()
        try:    
            from pyliftover import LiftOver     #gia to liftover
            lo=LiftOver('hg18', 'hg38')         # downloads the hg18-to-hg38 coordinate conversion chain file from UCSC
        except ModuleNotFoundError:
            print('Please install pyLiftover. Try "--help" for more information.')            
            sys.exit()
            
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

#%%                          LD
if args.LD is not None:
    snpA, snpB= args.LD
    if snpA[:4]=='snp_' and snpB[:4]=='snp_':
        
        def diplotupoi(snplist):
            '''Ypologismos twn diplotupwn pou xreiazontai gia na pragmatopoih8ei to haplotype estimation wste na upologistei to LD 2 thesewn
            Input: lista me 2 antikeimena: grammes twn 2 thesewn se Genotype File Format  
            Output: tuple tis morfis (arithmos atomwn omozugwn gia to refrence kai stis 2 theseis, arithmos atomwn omozugwn gia to refrence stin prwti thesi kai eterozugwn stin deuteri, arithmos atomwn eterozugwn stin prwti thesi kai omozugwn gia to refrence stin deuteri, arithmos atomwn eterozugwn kai stis duo theseis, arithmos atomwn'''    
            snpA, snpB=snplist    
            snpA_splitted=snpA.split(' ')        
            snpB_splitted=snpB.split(' ')
            N=(len(snpA_splitted)-5)/3
            RR=0 ; Rh=0 ; 
            hR=0 ; hh=0 ; 
             
            for i in range(5, len(snpA_splitted) -len(snpA_splitted)%3 ,3):
                #oi gonotupoi tou kathe atomou stis 2 theseis (theseis a kai p) 
                a_individual= snpA_splitted[i] + snpA_splitted[i+1] + snpA_splitted[i+2]
                b_individual= snpB_splitted[i] + snpB_splitted[i+1] + snpB_splitted[i+2]
                
                if a_individual == '100':
                    if b_individual == '100':
                        RR+=1
                    elif b_individual == '010':
                        Rh+=1
                elif a_individual == '010':
                    if b_individual == '100':
                        hR+=1
                    elif b_individual == '010':
                        hh+=1
            return RR, Rh, hR, hh, N
        
        def EM(pAB,pA, pB, diplotype): #diplotupe: einai tuple apo tin sunartisi 
            '''Ypologismos tis ektimomenis aplotupikis suxnotitas ektelwntas ton EM algorithmo gia dyo snp
            *Prepei prwta na exei treksie i sunartisi diplotupoi gia ta snp pou theloume na meletisoume*
            Input: arxiki timi gia tin suxnotita tou aplotupou gia ton EM, suxnotita refrence allilomorfou stin mia thesi, suxnotita refrence allilomorfou stin alli thesi, tuple tis morfis (arithmos atomwn omozugwn gia to refrence kai stis 2 theseis, arithmos atomwn omozugwn gia to refrence stin prwti thesi kai eterozugwn stin deuteri, arithmos atomwn eterozugwn stin prwti thesi kai omozugwn gia to refrence stin deuteri, arithmos atomwn eterozugwn kai stis duo theseis    
            Output: ektimomeni aplotupiki suxnotita'''
            EMrun=True
            RR, Rh, hR, hh , N= diplotype
            
            while EMrun:
                E= (2*RR) + Rh + hR + (pAB * (1 + pAB - pA - pB) * hh)/ (((pA - pAB) * (pB - pAB)) + pAB * (1 + pAB - pA - pB))
                pAB_new= E/(2*N)    
                
                if abs(pAB_new - pAB)< 10**(-6):
                    EMrun = False
                else:
                    pAB=pAB_new
            return pAB_new


        with open(args.cases_file) as cases:
        #Apothikeusi twn lines twn thesewn pros meleti     
            snps=[]
            for line in cases:
                splittedline=line.split(' ')
                
                if splittedline[0] ==snpA or splittedline[0]==snpB:
                    line=line.rstrip('\n')
                    snps.append(line)
         
        #Ypologismos pA pB
        counts=list(map(genotype_counts, snps))
        alleles=list(map(allele_freq, counts))
        pA=alleles[0][1] ;pa=alleles[0][2]
        pb=alleles[1][2] ;pB=alleles[1][1] 
    
        import random
        if pA==1 or pB==1 or pa==1 or pb==1:        
            D=0
            Dtonos=0
            r2= 'cannot compute r2'
        else:    
            diplotype= diplotupoi(snps)
            pABs=[random.uniform(0,1) for x in range(100)]
            pAB_newS=[]
            for i in pABs:
                pAB_new= EM(i, pA, pB, diplotype)
                pAB_newS.append(pAB_new)
            pAB_mean= sum(pAB_newS)/len(pAB_newS)    
            
            D= pAB_mean - (pA * pB)
            r2= D**2 /( pA * pa * pB * pb)
            #Ypologismos Dtonos
            if D <0:
                lista= [pA*pB, (1-pA) * (1-pB)]
                Dmax= min(lista)                        
            else:
                lista=[pA*(1-pB), (1-pA)* pB]
                Dmax= min(lista)
            Dtonos= D/Dmax      
        
        with open('{}.ld'.format(args.output), 'w') as output:
            print(snpA, snpB, r2, Dtonos, file=output)
    else:
        print('Please provide two snp_IDs in the required format. Try "--help" for more information.')
   
#%%                         ASSOCIATION TEST
if args.association_test:
    
    def genotypic_association_test(x,y):                    #Giati AYTO elege stis odigies....
        '''Pragmatopoiei Genotypic Association test
        *Prepei prwta na treksi i genotype_counts gia ta dataset ksexwrista*
        Input: 2 tuple tis morfis (snp_ID, arithmos omozugwn atomwn gia to refrence allilomorfo, arithmos omozugwn atomwn gia to alternative allilomorfo, arithmos eterozugwn atomwn, sunolo atomwn, thesi SNP sumfwna me to NCBI build 36)
        X: twn CASES    Y: twn contols
        Output: tuple (snp_ID, thesi SNP sumfwna me to NCBI build 36, pvalue) tou test, odds ratio, odds ratio'''
        from scipy import stats            
        snp = x[0]
        loci = x[5]
        
        if x[1]!=0 and y[1]!=0 and x[2]!=0 and y[2]!=0 and x[3]!=0 and y[3]!=0:
            #snp, x2test_hw= HWE(x,y)
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
    
    def allelic_association_test(a,b): 
        '''Ypologismos tou allelic association test gia kathe SNP se controls kai cases.
        *Prepei prwta na treksoun i genotype_counts kai i major_minor*
        Input: 2 tuple tis morfis (snp_ID, arithmos  atomwn gia to major allilomorfo sta controls/cases, arithmos atomwn gia to minor allilomorfo sta controls/cases, arithmos eterozugwn atomwn sta controls/cases, sunolo atomwn, thesi SNP sumfwna me to NCBI build 36)
        a: major-minor countstwn CONTROLS   b: major-minor counts twn CASES
        Output: tuple (snp_ID, locus, pvalue tou association test, odds ratio metaksu major kai minor allele gia cases kai controls ) '''
        from scipy import stats 
        snp, major_controls_c, minor_controls_c, het_controls, N, locus = a
        snp, major_cases_c, minor_cases_c, het_cases, N, locus  = b
        
        #Ypologismos plithous atomwn pou feroun ta allilomorfa (major kai minor)
        counts_major_controls = 2*major_controls_c+het_controls
        counts_minor_controls = 2*minor_controls_c+het_controls
        counts_major_cases = 2*major_cases_c+het_cases
        counts_minor_cases = 2*minor_cases_c+het_cases
        
        #Apaitountai gia ton upologismo twn expected counts twn allilomorfwn
        sum_major = counts_major_controls + counts_major_cases
        sum_minor = counts_minor_controls + counts_minor_cases
        #Expected counts allilomorfwn   
        exp_major_cases= (sum_major * N)/(2*N)      #N=Ncases=Ncontrols
        exp_minor_cases=  (sum_minor * N)/(2*N)
        exp_major_controls= (sum_major * N)/(2*N)
        exp_minor_controls = (sum_minor * N)/(2*N)
                                     
        x2test_aa = stats.chisquare([counts_major_controls, counts_minor_controls ,counts_major_cases , counts_minor_cases], [exp_major_controls, exp_minor_controls, exp_major_cases, exp_minor_cases])
        if minor_cases_c != 0 and major_controls_c != 0:
            OR= (major_cases_c * minor_controls_c) / (minor_cases_c * major_controls_c)
            if OR < 1 and OR!=0:
                OR = 1/OR
        else:
            OR = "Nan"
            
        return snp, locus, x2test_aa.pvalue,OR
        
    with open(args.cases_file) as cases, open(args.controls_file) as controls, open('{}.association'.format(args.output), 'w') as output:    
        #gia to manhattan plot
        loci_list = []
        pvalue_list = []
    
        for line_cases, line_controls in zip(cases, controls):
            line_cases=line_cases.rstrip('\n')
            line_controls=line_controls.rstrip('\n')
               
            counts_cases = genotype_counts(line_cases)                
            counts_controls = genotype_counts(line_controls)
             
            a = major_minor(counts_cases)       #major-minor counts gia cases
            b = major_minor(counts_controls)    #major-minor counts gia controls
            
            snp, locus, pvalue, OR = allelic_association_test(a,b)
            loci_list.append(locus)
            pvalue_list.append(pvalue)
            print(snp, locus, pvalue, OR, file=output)

#%%                         MANHATTAN PLOT
if args.manhattan and args.association_test:
    
    def manhattan(loci_list, pvalue_list): 
        '''Dimiourgia manhattan plot apo ta pvalues ths HWE h' ths association test
        *Prepei na xei treksei h HWE h' h allelic_Association_test antistoixa*
        Input: 2 listes (1 me ta pvalues olwn twn snp (pvalue_list) kai 1 me tis thesis SNP sumfwna me to NCBI build 36)
        x: loci_list   y: pvalue_list
        Outupt: manhattan plot'''
        import matplotlib.pyplot as plt
        import math
        
        pvalues = list(map(lambda x:(-math.log(x)),pvalue_list))        #efarmozw ton arnhtiko logari8mo sola ta pvalues
        plt.plot(loci_list, pvalues, ls='', marker='.', color='green')
        #plt.xlim([0,int(loci_list[-1])])#oria aksona x,mia 8esh akrivws meta to telutaio snp
        plt.title("Association Manhattan Plot", fontsize=15)
        plt.xlabel("Chromosome 20 position", fontsize=10)
        plt.ylabel("-log(Pvalue)", fontsize=10)
        plt.savefig("manhattan.jpg")
        plt.show()
    
    manhattan(loci_list, pvalue_list)

elif args.manhattan and args.association_test==False:
    print('Can not creat manhattan plot if no association test is performed.')

#%%                             QQPLOT
if args.qqplot and args.association_test:
    
    def qqplot(pvalue_list): 
        '''Dimiourgia qq plot apo ta pvalues ths association test
        *Prepei na xei treksei h allelic_association_test *
        Input: 1 lista me ta pvalues olwn twn snp (pvalue_list) 
        Outupt: qq plot'''
        import statsmodels.api as sm
        import pylab
        import numpy as np
        import math
        pvalues = list(map(lambda x:(-math.log(x)),pvalue_list))#efarmozw ton arnhtiko logari8mo sola ta pvalues
        pv = np.asarray(pvalues)
        sm.qqplot(pv, line='s')# etoimo module gia qqplot
        pylab.title("QQplot")
        pylab.savefig("qqlot.jpg")
        pylab.show()

    qqplot(pvalue_list)
    
elif args.manhattan and args.association_test==False:
    print('Can not creat qqplot if no association test is performed.')









