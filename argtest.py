#!/usr/bin/python

#GIA TREKSIMO
#python argtest.py -controls_file data/mikrakicontrols.txt -cases_file data/mikrakicases.txt -output testaki -allele_frequency

import sys
import argparse
#dimiourgw ena ArgumentParser antikeimeno
parser = argparse.ArgumentParser(description='This is our projectara')

#prosthetw ta args. Vazw flags kai parametrous
#gia ta args: default type=string /default action=store
parser.add_argument('-controls_file', help='File containing control genotypes',required=True)
parser.add_argument('-cases_file', help='File containing cases genotypes',required=True)
parser.add_argument('-output', help='Output file name',required=True)

#ta flags gia ta options tou programmatos. An o xristis dinei ena, tote pairnei tin timi True
parser.add_argument('-allele_frequency', action='store_true')           #!ftiakse ta help
parser.add_argument('-HWE', action='store_true')
#TO LD TELEI 2 VALUES SNPL 
parser.add_argument('-association_test', action='store_true')
# ginetai dekto mono an exw kanei assosiation: parser.add_argument('-manhattan', action='store_true')
# paromoiws parser.add_argument('-qqplot', action='store_true')
# to -get_info SNP thelei timi SNP

args = parser.parse_args()      #exei type Namespace
#etsi i metavliti args.controls_file exei parei tin antistoixi timi pou dwsame apo CL 
#ara oi metavlites mou einai oi args.controls_file, args.cases_file, args.output
        
#gia na dw an exw toulaxiston ena apo ta options allelefr, HWE etc.        
if len(sys.argv)>=8:
    pass
else:
    print('option missing')

    
        #--------ORIZW SUNARTISEIS EDW?????----------#

if args.allele_frequency:    
    print('kane ta kopla sou')
    
    #........
    
#==============for testing=====================================================
print(args.allele_frequency)
print(len(sys.argv))
print(sys.argv)
controlsfile=open(args.controls_file)
for line in controlsfile:
        print(line[:10])
    
    
        





 