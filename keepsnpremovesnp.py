#COMMAND GIA TREKSIMO:
#python keepsnpremovesnp.py -controls_file data/controlsmikraki.txt -cases_file data/casesmikraki.txt -output testaki -keep_snps data/snplist.txt
#to snplist.txt exei kwdikous snp ana grammmi. snp_2 '\n' snp_5

#einai apo to arxeio argtest 
#einai mono oi grammes pou xreiazontai gia na treksw to arxeio
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='This is our projectara')
parser.add_argument('-controls_file', help='File containing control genotypes',required=True)
parser.add_argument('-cases_file', help='File containing cases genotypes',required=True)
parser.add_argument('-output', help='Output file name',required=True)
#diko tou flag:
parser.add_argument('-keep_snps')           #!ftiakse ta help
parser.add_argument('-remove_snps')
parser.add_argument('-keep_samples')
parser.add_argument('-remove_samples')
args = parser.parse_args()

#%%
def createSNPlist(file):
    '''input: arxeio me snps ana seira /output: lista twn snp'''
    
    with open(file) as snpfile:
        snplist=[]
        for line in snpfile:        #ftiaxnw mia lista me ta snp pou thelw na paiksw mpala
            line=line.rstrip('\n') 
            snp= line.split(' ')[0] #check GIATI to ekanes auto
            snplist.append(snp)
    snplist.sort()  #giati mporei o xristie na exei dwsei ta snp me anakati seira      
    return snplist

#%%                         KEEP SNPS

if args.keep_snps is not None:  #elegxw an o xristis exei dwsei to flag -keep_snps kai ena arxeio 
    snplist=createSNPlist(args.keep_snps) 
    with open(args.cases_file) as cases, open('{}.cases.gen'.format(args.output), 'w') as output_cases:
        for line_cases in cases:
            if line_cases.split(' ')[0] in snplist:
                output_cases.write(line_cases)
                
    with open(args.controls_file) as controls, open('{}.controls.gen'.format(args.output), 'w') as output_controls:
        for line_controls in controls:
            if line_controls.split(' ')[0] in snplist:
                output_controls.write(line_controls)
            
#==================ALLOS TROPOS============================================================
#             for snp in createSNPlist(args.keep_snps):
#                 if  snp in line_cases:  #mporw kai:  if  line_cases.split(' ')[0]==snp:
#                     output_cases.write(line_cases)
#                     break   #min koitakseis alles times gia snp apo tin lista
#==============================================================================
                   
#%%                        REMOVE SNPS

if args.remove_snps is not None:
    snplist=createSNPlist(args.remove_snps) 
    with open(args.cases_file) as cases, open('{}.cases.gen'.format(args.output), 'w') as output_cases:
        for line_cases in cases:
            if line_cases.split(' ')[0] not in snplist:
                output_cases.write(line_cases)
    with open(args.controls_file) as controls, open('{}.controls.gen'.format(args.output), 'w') as output_controls:
        for line_controls in controls:
            if line_controls.split(' ')[0] in snplist:
                output_controls.write(line_controls)
                 
#%%                 SYNARTISI 
#xreiazetai arxeio pou exei samples ana grammi (px) case_1 '\n' control_13 etc.     
def get_samples(file):
    '''Input:  Ένα αρχείο με κωδικούς sample ανά γραμμή. Ένα sample έχει κωδικό: (case ή control)αριθμός. Π.χ. το sample control_5 είναι το 5ο control. Το case_100 είναι το 100ό case.
    Output: mia lista me tous kwdikous pou anaferontai se deigmata cases kai mia me kwdikous controls'''
    
    with open(file) as samplelist:
        case_samples=[]    #tha ftiaksw duo listes kai gia tin kathe mia tha apothikeuw ton aritho tou sample(stin antistoixi lista) 
        control_samples=[]
        for line in samplelist:
            line=line.rstrip('\n')
            if line.split('_')[0] == 'case':
                case_samples.append(int(line.split('_')[1]))   #int giati to diavazei san string 
            elif line.split('_')[0] == 'control':
                control_samples.append(int(line.split('_')[1]))    
        return case_samples, control_samples
    
    
#%%             KEEPSAMPLES -SPYDER
#to idio prepei na kanw kai gia ta cases me tin alli lista
#trexw to DATAmpourdaki: einai upergamato giati mou dexnei akrivws ti pairnw (kai oxi midenika kai asous)
#==============================================================================
# control_samples= sorted(get_samples('/home/rantaplan/master/projectara/data/sampleslist.txt')[1])
# 
# with open('/home/rantaplan/master/projectara/DATAmpourdaki') as controls:
#     for line_controls in controls:
#         line_controls=line_controls.rstrip('\n')
#         splittedline= line_controls.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
#         
#         newline= splittedline[0:5]  #einai lista!!
#         for i in control_samples:
#             index=(i-1)*3+5
#             individual= splittedline[index]+ splittedline[index+1]+splittedline[index+2]
#             newline+=individual
#         newline_string=' '.join(newline) + '\n'
#         print(newline_string)       #printare to se file
# 
#==============================================================================
#%%                   keep  SAMPLES -ARGS
if args.keep_samples is not None:    
    control_samples= sorted(get_samples(args.keep_samples)[1])
    case_samples= sorted(get_samples(args.keep_samples)[0])
    
    with open(args.controls_file) as controls, open('{}.controls.gen'.format(args.output), 'w') as output_controls:
        for line_controls in controls:
            line_controls=line_controls.rstrip('\n')
            splittedline_controls= line_controls.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
            newline_controls= splittedline_controls[0:5]  #einai lista!!
            for i in control_samples:
                index=(i-1)*3+5
                individual= splittedline_controls[index]+ splittedline_controls[index+1]+splittedline_controls[index+2]
                newline_controls+=individual
            newline_controls_string=' '.join(newline_controls) + '\n'
            output_controls.write(newline_controls_string) 
    
    with open(args.cases_file) as cases, open('{}.cases.gen'.format(args.output), 'w') as output_cases:
        for line_cases in cases:
            line_cases=line_cases.rstrip('\n')
            splittedline_cases= line_cases.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
            newline_cases= splittedline_cases[0:5]  #einai lista!!
            for i in case_samples:
                index=(i-1)*3+5
                individual= splittedline_cases[index]+ splittedline_cases[index+1]+splittedline_cases[index+2]
                newline_cases+=individual
            newline_cases_string=' '.join(newline_cases) + '\n'
            output_cases.write(newline_cases_string) 

#%%                 REMOVE SAMPLES-SPYDER
#============================================================================== 
# control_samples= sorted(get_samples('/home/rantaplan/master/projectara/data/sampleslist.txt')[1], reverse=True)
# with open('/home/rantaplan/master/projectara/DATAmpourdaki') as controls:
#     for line_controls in controls:
#         line_controls=line_controls.rstrip('\n')
#         splittedline= line_controls.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings        
#         for i in control_samples:
#             index=(i-1)*3+5
#             del splittedline[index:index+3] #einai san to range, paei mexri to proigoumeno stoixeio tou anw oriou
#             
#         newline=' '.join(splittedline) + '\n'
#         print(newline)           
#==============================================================================
#%%                 REMOVE SAMPLES-ARGS

if args.remove_samples is not None:    
    control_samples= sorted(get_samples(args.remove_samples)[1], reverse=True)
    case_samples= sorted(get_samples(args.remove_samples)[0], reverse=True)
    
    with open(args.controls_file) as controls, open('{}.controls.gen'.format(args.output), 'w') as output_controls:
        for line_controls in controls:
            line_controls=line_controls.rstrip('\n')
            splittedline_controls= line_controls.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
            for i in control_samples:
                index=(i-1)*3+5
                del splittedline_controls[index:index+3] #einai san to range, paei mexri to proigoumeno stoixeio tou anw oriou
            newline_controls_string=' '.join(splittedline_controls) + '\n'
            output_controls.write(newline_controls_string) 
    
    with open(args.cases_file) as cases, open('{}.cases.gen'.format(args.output), 'w') as output_cases:
        for line_cases in cases:
            line_cases=line_cases.rstrip('\n')
            splittedline_cases= line_cases.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
            for i in case_samples:
                index=(i-1)*3+5
                del splittedline_cases[index:index+3] #einai san to range, paei mexri to proigoumeno stoixeio tou anw oriou
            newline_cases_string=' '.join(splittedline_cases) + '\n'
            output_cases.write(newline_cases_string) 


