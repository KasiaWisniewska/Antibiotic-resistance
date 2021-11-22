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
		gene_names = []
		dna = ''
		poop = []
		for line in filein:
			if line[0] == '>':
				gene_names.append(line[:-1])
				if dna != '':
					genes.append(dna)
				dna = ''
			else:
				dna += line[:-1]
		genes.append(dna)
	return [genes, gene_names]

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

def match_ngs_to_gene(ngsread_kmer, gene):
	### returns the list of depths of nucleotides that matched for the gene
	match = False

	for i in range(len(ngsread_kmer)):
		k = len(ngsread_kmer[0])			# k as in k-mer, usually 19

		# see if the k-mer is found in the gene and find its pos
		start_gene = gene.find(ngsread_kmer[i])	
		if start_gene != -1:			# -1 if the k-mer is not found
			start_read = i

			if start_gene > 0:
				# start_read must be 0
				if start_read == 0:
					
					# list with zeroes, length = len(gene)
					depth = [0 for j in range(len(gene))]
					print('before',len(depth))				####### print for visualization, delete later
					# depth = 1 for the first matched k-mer
					depth[start_gene:start_gene+k] = [1 for j in range(k)]

					pos_read = start_read + 1
					pos_gene = start_gene + 1
					while pos_read < len(ngsread_kmer) and pos_gene < len(gene)-k:
						if ngsread_kmer[pos_read] == gene[pos_gene:pos_gene+k]:
							depth[pos_gene:pos_gene+k] = [1 for j in range(k)]		# optimize later?
							#depth[pos_gene+k-1] = 1
							pos_read += k
							pos_gene += k
						
						else:
							print('after',len(depth))		####### print for visualization, delete later
							return depth




				else:
					return None

			elif start_gene == 0:
				pass



	return None


#	while kmer[startpos] not in gene_list:
#		startpos += 1
#		if startpos == 



### Main Program ###

k = 19					# k-mer size, usually 19

[gene_list, gene_names] = get_genes('resistance_genes.fsa')
gene_kmer = get_gene_kmer(gene_list, k)
#print(gene_list[0], len(gene_kmer))
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
					for gene_num in range(len(gene_list)):
						depth_update = match_ngs_to_gene(ngsread_kmer, gene_list[gene_num])
						if depth_update is not None:
							print(gene_names[gene_num], '|', countin, '| Read #', count)		####### print for visualization, delete later
							print(depth_update, '\n')						####### print for visualization, delete later
					break					# takes 2m9s

#			ngsread_kmer = get_ngsread_kmer_set(dna_read, k)		
#			if len( ngsread_kmer.intersection(gene_kmer) ) > 0:		# or > 1?
#				countin += 1				# takes 1m59.9s
#				ngsread_kmer = get_ngsread_kmer_list(dna_read, k)	# takes 2m9s

print(countin, '/', count)					# 6,717 / 3,469,171			
