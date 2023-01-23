#!/usr/bin/env python 

############################################################
### MODULES ################################################

import sys
from optparse import OptionParser
import re

############################################### /MODULES ###
############################################################




############################################################
### FUNCTIONS ##############################################

#def allowed_cigar(p,cigar):
#	return not p.search(cigar) 

def allowed_strata(strata):
	return strata in ['1','0:1','0:0:1','0:0:0:1','0:0:0:0:1']

############################################# /FUNCTIONS ###
############################################################




############################################################
### ARGUMENTS,OPTIONS ######################################

parser = OptionParser(usage="\n%prog [options]", version="%prog 0.1")
parser.add_option("-i", metavar="FILE", type="string", dest="input_file", default = 'stdin', help="input gem alignment filename (default = 'stdin')")
parser.add_option("-o", metavar="FILE", type="string", dest="output_file", default = 'reads', help="output filename prefix (default = 'reads')")

(opt, args) = parser.parse_args()
if opt.input_file==None:
	parser.print_help()
	sys.exit(-1)
        
##################################### /ARGUMENTS,OPTIONS ###
############################################################




############################################################
### CONSTANTS ##############################################

 
############################################# /CONSTANTS ###
############################################################




############################################################
### MAIN ###################################################

of1=open(opt.output_file+'.1.fastq','w')
of2=open(opt.output_file+'.2.fastq','w')

if opt.input_file!='stdin':
	f=open(opt.input_file,'r')
else:
	f=sys.stdin


while True:
    line1 = f.readline()
    if not line1:
    	break  # EOF
    line2 = f.readline()
    line3 = f.readline()
    line4 = f.readline()
    print>>of1,line1.strip()
    print>>of1,line2.strip()
    print>>of1,line3.strip()
    print>>of1,line4.strip()
    
    line5 = f.readline()
    line6 = f.readline()
    line7 = f.readline()
    line8 = f.readline()
    print>>of2,line5.strip()
    print>>of2,line6.strip()
    print>>of2,line7.strip()
    print>>of2,line8.strip()


if opt.input_file!='stdin':
	f.close()

of1.close()
of2.close()
	
################################################ /MAIN ###
############################################################
