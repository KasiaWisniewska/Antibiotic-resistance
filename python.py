#!/usr/bin/env python3
'''
Res genes: Dict, key = kmer, value = positions CHECK
Then compare kmers against it (first and last of the read)
Two files need to work at a time!
Create reverse complementary strands of reads!

To possibly improve:
in read_filename() make it possible for as many files to be added as we want
'''
import sys
import gzip

def read_filename():		
	# read the filename
	if len(sys.argv) == 3:
		f1 = sys.argv[-2]
		f2 = sys.argv[-1]
	else:
		f1 = input('Enter <filename>: ')
		f2 = input('Enter another <filename>: ')
	filenames = [f1, f2]
	return filenames
	
	
def get_ngsread_kmer_list(dna, k):
	### returns a LIST of all k-mers extracted from an NGS read in the 
	### respective order

	ngsread_kmer = []
	for i in range(len(dna)-k):				# len-k ignores \n at the end
		ngsread_kmer.append(dna[i:i+k])
	return ngsread_kmer
	
def get_genes(filename):
	with open(filename, 'r') as infile:
		genes = []
		gene_names = []
		dna = ''
		for line in infile:
			if line[0] == '>':
				gene_names.append(line[:-1])
				if dna != '':
					genes.append(dna)
				dna = ''
			else:
				dna += line[:-1]
		genes.append(dna)
	return genes, gene_names
'''	
def complementary_strand(read):
	### returns a list of original genes and reverse complementary genes together
	transTable = str.maketrans('ATCG', 'TAGC')
	for gene in genes:
		complementDNA = read.translate(transTable)
		rev_dna = complementDNA[::-1]
		read.append(rev_dna)  
	
	return read
'''

def complementary_strand(read):
	### returns a reverse complementary version of the read
	transTable = str.maketrans('ATCG', 'TAGC')
	complementDNA = read.translate(transTable)
	rev_dna = complementDNA[::-1]  
	
	return rev_dna

#a = read_filename()	#delete
#print(a[1])			#delete

def get_gene_kmer(genes, k):
	### returns a set of all k-mers taken from all resistance genes
	gene_kmer_dict = dict()
	
	for gene_num in range(len(genes)):
		for i in range(len(genes[gene_num])-k+1):		# len-k+1 takes the last symbol too
			check = gene_kmer_dict.get(genes[gene_num][i:i+k])
			if check is not None: 
				gene_kmer_dict[genes[gene_num][i:i+k]].append([gene_num, i])
			else:
				gene_kmer_dict[genes[gene_num][i:i+k]] = [[gene_num, i]]	# value = kmer position [which gene][where in the gene]

	return gene_kmer_dict




'''MAIN CODE'''
	
k = 19 #setting the kmer length

###processing resistance genes

#Obtaining lists of resistance gene sequences and names
[gene_list, gene_names] = get_genes('resistance_genes.fsa')

#Obtaining a dict of res_gene kmers with their respective locations as values
gene_dict = get_gene_kmer(gene_list, k)

#print(gene_dict)	#check - delete


###creating a template for output

#Creating depth list
depth = []
for i in range(len(gene_list)):
	depth.append([])
	for j in range(len(gene_list[i])):
		depth[i].append(0)

###processing the reads

#opening both read files

filenames = read_filename()

for file in filenames:	#to iterate through both files
	with gzip.open(file, "r") as infile:
		for b_line in infile:
			line = b_line.decode('ASCII')
			read_count = 0
			if line[0] == '@':
				flag = True
			elif flag:
				
				if line[0] == '+':
					flag = False

				else:
					read_count += 1
					# get the dna seq
					dna_read = line
					#get the reverse complement for the read
					rev_read = complementary_strand(dna_read)
					both_reads = [dna_read, rev_read]
					for read in both_reads:
						#create kmer lists for both
						read_kmer = get_ngsread_kmer_list(read, k)
						#checking if extremities of the read fit the dict (initial read elimination)
						start_check = gene_dict.get(read_kmer[0])
						end_check = gene_dict.get(read_kmer[-1])
						if start_check is not None or end_check is not None: 
							for kmer in read_kmer: 
								pos = gene_dict.get(kmer)
								if pos is not None: 
									for hit in pos:
										for i in range(hit[1],hit[1]+k):
											depth[hit[0]][i] += 1
'''				
				#reinitializing hit numbers (one hit for one read only)					
				for single_depth in depth:
					for position in single_depth:
						if position >= read_count:
							position = read_count
'''

print(depth)			
					
			
					
					
					
					
					
				
