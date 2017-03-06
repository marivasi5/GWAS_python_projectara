#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 23:58:56 2017

@author: rantaplan
"""
import numpy as np
import time
start=time.time()

with open('/home/rantaplan/master/projectara/convert_test.gen') as file, open('convert.map' , 'w') as output:
    all=[]      #lista me ta cases_list olwn twn lines
    for line in file:
        line.rstrip('\n')
        
        splittedline= line.split(' ') 
        snp=splittedline[0]
        locus=splittedline[2]
        ref=splittedline[3]
        alt=splittedline[4]
        
        #gia to map
        print( '20 {} 0 {}'.format(snp, locus), file=output)
        
        cases_list=[]   #lista me gonotupous twn cases ana index  (ksexwristi gia to line pou eimaste) 
        
        for i in range(5, len(splittedline) -len(splittedline)%3, 3): #!tsekare to range
            individual= splittedline[i] + splittedline[i+1] + splittedline[i+2]
            if individual=='100':
                cases_list.append(ref+ ' '+ref)
            elif individual=='010':
                cases_list.append(ref+ ' '+alt)
            elif individual=='001':
                cases_list.append(alt+ ' '+alt)
                
        all.append(cases_list)

pinakas_snps_cases=np.array(all)
pinakas_cases_snps= pinakas_snps_cases.transpose()
ped_list= pinakas_cases_snps.tolist()

print(time.time() -start)

with open('converted.ped', 'w') as output:
    
    for i in range(0, len(ped_list)):
        start= '0 {} 0 0 0 0 '.format(i)
        genotypes=' '.join(ped_list[i])
        new_line=start + genotypes + '\n'
        print(new_line, file=output)
        
    