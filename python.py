#!/usr/bin/env python3

import sys
import gzip

def read_filename():		
	# read the filename
	if len(sys.argv) == 2:
		filename = sys.argv[-1]
	else:
		filename = input('Enter <filename>: ')
	return filename 

def get_genes(filename):
	with open(filename, 'r') as filein:
		genes = []
		dna = ''
		for line in filein:
			if line[0] == '>':
				if dna != '':
					genes.append(dna)
				dna = ''
			else:
				dna += line[:-1]
	return genes

def get_gene_kmer(genes, k):
	### returns a set of all k-mers taken from all resistance genes
	gene_kmer = set()
	for gene in genes:
		for i in range(len(gene)-k+1):		# len-k+1 takes the last symbol too
			gene_kmer.add(gene[i:i+k])		# [i:i+19] -- 19-mer
	return gene_kmer

def get_ngsread_kmer_list(dna, k):
	### returns a LIST of all k-mers extracted from an NGS read in the 
	### respective order

	ngsread_kmer = []
	for i in range(len(dna)-k):				# len-k ignores \n at the end
		ngsread_kmer.append(dna[i:i+k])
	return ngsread_kmer


def get_ngsread_kmer_set(dna, k):
	### returns a SET of all k-mers extracted from an NGS read

	ngsread_kmer = set()
	for i in range(len(dna)-k):
		ngsread_kmer.add(dna[i:i+k])
	return ngsread_kmer



### Main Program ###

k = 19					# k-mer size

gene_list = get_genes('resistance_genes.fsa')
gene_kmer = get_gene_kmer(gene_list, k)
print(gene_list[0], len(gene_kmer))
count = 0					# counts dna reads in the NGS file
countin = 0					# counts dna reads with at least one k-mer match

filename = read_filename()
filein = gzip.open(filename, 'r')
for b_line in filein:
	line = b_line.decode('ASCII')

	if line[0] == '@':
		flag = True
	elif flag:
		if line[0] == '+':
			flag = False

		else:
			# get the dna seq
			dna_read = line
			count += 1

			# check if at least 1 k-mer from read is present in genes
			ngsread_kmer = get_ngsread_kmer_list(dna_read, k)
			for kmer in ngsread_kmer:
				if kmer in gene_kmer:
					countin += 1
					break					# takes 2m9s

#			ngsread_kmer = get_ngsread_kmer_set(dna_read, k)		
#			if len( ngsread_kmer.intersection(gene_kmer) ) > 0:		# or > 1?
#				countin += 1				# takes 1m59.9s
#				ngsread_kmer = get_ngsread_kmer_list(dna_read, k)	# takes 2m9s

print(countin, '/', count)					# 6,717 / 3,469,171			
