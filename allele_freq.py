#GIA TREKSIMO
#python argtest.py -controls_file data/mikrakicontrols.txt -cases_file data/mikrakicases.txt -output testaki -allele_frequency

import sys
import argparse

#einai apo to arxeio argtest. einai mono oi grammes pou xreiazontai gia na treksw to allelefreq
parser = argparse.ArgumentParser(description='This is our projectara')
parser.add_argument('-controls_file', help='File containing control genotypes',required=True)
parser.add_argument('-cases_file', help='File containing cases genotypes',required=True)
parser.add_argument('-output', help='Output file name',required=True)
parser.add_argument('-allele_frequency', action='store_true')           #!ftiakse ta help
args = parser.parse_args()

#%%     Synartisi
def suxnotita(datasetLINE):    #prepei na to taiso lines    
    splittedline= datasetLINE.split(' ')   #einai lista, ta items tis anagnwrizontai ws strings
    snp= splittedline[0]
            
    homozygous_refrence=0
    homozygous_alternative=0
    heterozygous=0
    for i in range(5, len(splittedline), 3): #!tsekare to range
                individual= splittedline[i] + splittedline[i+1] + splittedline[i+2]
                if individual=='100':
                    homozygous_refrence+=1
                elif individual=='010':
                    heterozygous+=1
                elif individual=='001':
                    homozygous_alternative+=1
        #ypologismos suxnotitwn:
                p= round((2*homozygous_refrence + heterozygous)/1000 , 3)  #einai /2N opou N=500   
                q= round(1-p , 3)     
                    
                
    return snp+'\t'+ str(p) +'\t'+  str(q)        #na dw ta rounds

    #an to eixa se return snp, p, q to ekane tuple! /an to eixa se print meta i epomeni print mmou ekane nera (none)  
#%%                 START
if args.allele_frequency:
    
    # YPARXEI LATHOS ME TO GRAPSIMO SE ARXEIO outputfile=open(args.output, 'w')
            
    with open(args.controls_file) as cases, open(args.controls_file) as controls: 
        for line_cases, line_controls in zip(cases, controls):
            line_cases=line_cases.rstrip('\n')
            line_controls=line_controls.rstrip('\n')        #!to suxnotita(line_cases) exei type: NoneType kaii otan to kanw str() mou typwnei ena None san deuteri grammi
            
            print(suxnotita(line_cases), suxnotita(line_controls)[6:]) 
            
            #file=outputfile
            #outputfile.close()  
    


























    
    
    