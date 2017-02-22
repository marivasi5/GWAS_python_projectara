
#%%                                                         TEST LIFTOVER
rs= 'rs4814683'
thesi= '9795'

from pyliftover import LiftOver
lo=LiftOver('hg18', 'hg38')
lo.convert_coordinate('chr20', 9795)
#[('chr20', 81154, '+', 5643036713)]

#%%                     GET kantale
import requests
url = 'http://rest.ensembl.org/vep/human/hgvs/9:g.22125504G>C?'
headers={ "Content-Type" : "application/json"}

r = requests.get(url, headers=headers)
print (r.ok)
data = r.json()

#==============================================================================
# to data einai lista me ena antikeimeno
# pou einai dictionary(11 elements) 
# uparxei key: transcript_consequences 
# me value mia lista 
# me 11 antikeimena-kathe ena einai dictionary!!! (diaferon metaksi tous sto transcript id)
#==============================================================================

#%%
#input snp
snp='snp_9'

#get position
with open('/home/rantaplan/master/projectara/data/casesmikraki.txt') as cases:
    for line in cases:
        splittedline= line.split(' ')
        if splittedline[0] == snp:
            position = int(splittedline[2])
            refrence= splittedline[3]
            alternative= splittedline[4]
#liftover    
from pyliftover import LiftOver
lo=LiftOver('hg18', 'hg38')
conversion= lo.convert_coordinate('chr20', position)
position_new=conversion[0][1]

#get request
import requests
url= 'http://rest.ensembl.org/vep/human/hgvs/20:g.' + str(position_new) + refrence + '>' + alternative + '?' #to xromoswma twra einai to 20 alla tha eprepe na to vazw

#url = 'http://rest.ensembl.org/vep/human/hgvs/9:g.15507T>C?'    #snp9
headers={ "Content-Type" : "application/json"}
r = requests.get(url, headers=headers)

if not r.ok:
    print('The request was not successful.')    
    r.raise_for_status()    #apo tis odigies tou API
    sys.exit()
 
else:
    data = r.json()
    #!!!!!!!!!!!!!!!!to dictionary exei diaforetika keys apo to paradeigma
    #!!!!!!!!!!!!!!!!!!!!!!!!!!MPALES!!!!!!!!!!!!!!!!!!!























