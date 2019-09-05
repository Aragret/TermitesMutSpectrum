#!/home/tools/bin/python

import sys
import re
import subprocess

file_name = sys.argv[1]
gene_name = sys.argv[2]

with open(file_name) as infile: #, open('temp_file', 'w') as outfile:
	#while True:
	for i in range(10):
        	name = infile.readline().strip()[1:-1]
        	bases = infile.readline().strip()

        	sp_name = name.split('I") ')[1].split(' mitochondrion')[0]
        	species_name = re.sub(' ', '_', sp_name)

        	with open('nuc_seq.fa', 'w+') as nucl:
        		nucl.write('>{0}\n{1}\n'.format(species_name, bases))
        	
        	subprocess.call("~/bin/emboss/bin/transeq -sequence nuc_seq.fa -outseq prot_seq.pep -table 5", shell=True)

        	with open("prot_seq.pep", 'r') as prot:
        		protein = ''
        		for line in prot:
        			if line[0] != '>':
        				protein += line.strip()
		print(species_name, protein, gene_name)

        	subprocess.call("SCRIPT2.2.1.sh species_name protein gene_name", shell=True)

        	#if not bases:
            		#break


