with open('/home/rantaplan/master/togamatoproject/data/mikrakicases.txt') as cases, open('/home/rantaplan/master/togamatoproject/data/mikrakicontrols.txt') as controls, open('/home/rantaplan/master/togamatoproject/data/snplist.txt') as snpfile:

    snplist=[]
    for line in snpfile:        #ftiaxnw mia lista me ta snp pou thelw na paiksw mpala
        line=line.rstrip('\n') 
        snp= line.split(' ')[0]
        snplist.append(snp)
    
#%%             
def createSNPlist(file):
    '''input: arxeio me snps ana seira /output: lista twn snp'''
    
    with open(file) as snpfile:
        snplist=[]
        for line in snpfile:        #ftiaxnw mia lista me ta snp pou thelw na paiksw mpala
            line=line.rstrip('\n') 
            snp= line.split(' ')[0]
            snplist.append(snp)
    return snplist

#%% CALLING    
createSNPlist('/home/rantaplan/master/togamatoproject/data/snplist.txt')
#%%          KEEP SNP
with open('/home/rantaplan/master/togamatoproject/data/mikrakicases.txt') as cases, open('/home/rantaplan/master/togamatoproject/data/mikrakicontrols.txt') as controls:
    f = open('output', 'w')
    for line_cases in cases:
        for snp in createSNPlist('/home/rantaplan/master/togamatoproject/data/snplist.txt'):
            if  line_cases.split(' ')[0]==snp:
                f.write(line_cases)
                break   #min koitakseis alles times gia snp apo tin lista
            
#==============================================================================
#     for line_controls in controls:    
#         line_controls=line_controls.rstrip('\n')    
#==============================================================================




























