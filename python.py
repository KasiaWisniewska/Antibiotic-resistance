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
#			if check is not None: 
#				gene_kmer_dict[genes[gene_num][i:i+k]].add({gene_num: i})
#			else:
			gene_kmer_dict[genes[gene_num][i:i+k]] = {gene_num: i}	# value = kmer position [which gene][where in the gene]

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
read_count = 0
for file in filenames:	#to iterate through both files
	with gzip.open(file, "r") as infile:
		for b_line in infile:
			line = b_line.decode('ASCII')
			if line[0] == '@':
				flag = True
			elif flag:
				
				if line[0] == '+':
					flag = False

				else:
					read_count += 1
					#if read_count > 500000:
					#	break
					# get the dna seq
					dna_read = line
					#get the reverse complement for the read
					rev_read = complementary_strand(dna_read)
					both_reads = [dna_read, rev_read]
					for read_n in both_reads:
						read = read_n[:-1]
						#create kmer lists for both
						read_kmer = get_ngsread_kmer_list(read, k)
						#checking if extremities of the read fit the dict (initial read elimination)
						start_check = gene_dict.get(read_kmer[0])
						end_check = gene_dict.get(read_kmer[-1])
						# list of gene numbers to which the read will be aligned
						genes_to_analyze = set()
						if start_check is not None:
							print(read_count)
							for hit in start_check.keys():
								genes_to_analyze.add(hit)
						if end_check is not None:
							for hit in end_check.keys():
								genes_to_analyze.add(hit)

							for gene_num in genes_to_analyze:

								if start_check is not None and start_check.get(gene_num) is not None:
									if end_check is not None and end_check.get(gene_num) is not None:
										#print(read)
										#print(gene_list[gene_num][ start_check[gene_num]:end_check[gene_num]+k+1 ])
										#print(read == gene_list[gene_num][ start_check[gene_num]:end_check[gene_num]+k+1 ])
										#print()
										if read == gene_list[gene_num][ start_check[gene_num]:end_check[gene_num]+k+1 ]:
											update_reg = [start_check[gene_num], end_check[gene_num]+k+1]
										else:
											update_reg = None
									else:
										end_pos = 0
										check_next = gene_dict.get(read_kmer[1])
										if check_next is None or check_next.get(gene_num) is None:
											check_indicator = None
										else:
											check_indicator = check_next.get(gene_num)
										while check_indicator is not None:
											end_pos += 1
											check_next = gene_dict.get(read_kmer[end_pos+1])
											if check_next is None or check_next.get(gene_num) is None:
												check_indicator = None
											else:
												check_indicator = check_next.get(gene_num)
										check = gene_dict.get(read_kmer[end_pos])
										print(end_pos, end_pos+k)
										print(read[:end_pos+k], len(read[:end_pos+k]))
										#print(gene_list[gene_num][ start_check[gene_num]:check[gene_num]+k], len(gene_list[gene_num][ start_check[gene_num]:check[gene_num]+k]))
										print(gene_list[gene_num][ start_check[gene_num]: ], len(gene_list[gene_num][ start_check[gene_num]: ]))
										if read[:end_pos+k] == gene_list[gene_num][ start_check[gene_num]: ]:
											print(True)
											input()
										else:
											print(False)
										#if read[:end_pos+k] == gene_list[gene_num][ start_check[gene_num]: ]:
										#	update_reg = [start_check[gene_num], check[gene_num]+k]
										#else:
										#	update_reg = None



								else:
									if end_check is not None and end_check.get(gene_num) is not None:
										start_pos = 0
										check_next = gene_dict.get(read_kmer[1])
										if check_next is None or check_next.get(gene_num) is None:
											check_indicator = None
										else:
											check_indicator = check_next.get(gene_num)
										while check_indicator is None:
											start_pos += 1
											check_next = gene_dict.get(read_kmer[start_pos+1])
											if check_next is None or check_next.get(gene_num) is None:
												check_indicator = None
											else:
												check_indicator = check_next.get(gene_num)
										if read[start_pos:] == gene_list[gene_num][ check_next[gene_num]:end_check[gene_num]]:
											update_reg = [start_check[gene_num], check[gene_num]+k]
										else:
											update_reg = None
									else:
										update_reg = None


								if update_reg is not None:
									for i in range(update_reg[0], update_reg[1]):
										depth[gene_num][i] += 1
								

print(depth)

coverageCount = 0
coverage = dict()
for gene_num in range(len(depth)):
	hitCount = 0
	for nt in depth[gene_num]:
		if nt >= 10:
			hitCount += 1
	a = hitCount/len(depth[gene_num])
	if a > 0.95:
		coverageCount += 1
	coverage[gene_num] = a

print(coverage)

sorted_gene_num = sorted(coverage.keys(), key=coverage.get, reverse=True)

print(coverageCount, " genes have achieved coverage above 95%, meaning they are very likely present in the sample.", sep='')

res_list =	[["Beta-lactam", 0], ["Phenicol", 0], ["Aminoglycoside", 0], 
		["Tetracycline", 0], ["Fluoroquinolone and aminoglycoside", 0], 
		["Sulphonamide resistance", 0], ["Fosfomycin resistance", 0]]
for gene_num in sorted_gene_num[:coverageCount]:
	print(gene_names[gene_num], "; gene coverage: ", coverage[gene_num]*100, "%", sep='')	#To see the gene names that are covered
	print(depth[gene_num])																	#To see corresponding depth distribution 
	for res in res_list:
		if res[0] in gene_names[gene_num]:
			res[1] += 1

print("Many of these genes encode resistance for the same antibiotics, notably for:")			
for res in res_list:
	print(res[0], "resistance :", res[1], "genes present.")									#Shows only how many genes for which resistance are present

#more_info = Input("Do you wish to see the gene names (type g) or the gene depth distribution (type d)? If not: type 'n': ")

#for gene_num in sorted_gene_num[:coverageCount]:
#	if more_info == "g":
#		print(gene_names[gene_num], "; gene coverage: ", coverage[gene_num]*100, "%", sep='')
#	elif more_info == "d":
#		print(depth[gene_num])
#	else:
#		break
	
			
					
					
					
					
					
