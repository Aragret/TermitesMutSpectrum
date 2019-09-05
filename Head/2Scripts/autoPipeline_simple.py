#!/home/tools/bin/python

import sys
import re
import subprocess

file_name = sys.argv[1] #this is RefSeq CDS (genes) file containing one gene type, for example filtered (by gene names considering synonyms) GenesRefSeqs.fasta file
gene_name = sys.argv[2] #this is short gapless gene name, for example CYTB

with open(file_name) as File:
    text = File.read()
    Sequences = text.split("\n>")

for Seq in Sequences:
	lines = Seq.splitlines()
	sp_name = lines[0].split('") ')[1].split(' mitochondrion')[0]
	species_name = re.sub(' ', '_', sp_name)

	with open('nuc_seq.fa', 'w+') as nucl:
		nucl.write('>{0}\n'.format(Seq))

        subprocess.call("~/bin/emboss/bin/transeq -sequence nuc_seq.fa -outseq prot_seq.pep -table 5", shell=True)

        prot = open ("prot_seq.pep", "r")
        lines = prot.readlines()
        prot.close()
        protein = "".join(map(str.strip, lines[1:]))
        protein_clear = re.sub('\*', '', protein)

        subprocess.call('./SCRIPT2.2.1.sh %s %s %s' % (species_name, protein_clear, gene_name), shell=True)
