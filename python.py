#!/usr/bin/env python3

#open the resistance gene file
	#extract all genes into a dict variable
	

#open gzip file
#for line in file
	#find the sequence between @ and + lines
	#extract all 19-mers into a list

	#for 


""" ### Just playing with 'gzip' library

import gzip

with gzip.open('file.txt.gz', 'w') as fileout:
	for i in range(10):
		# has to be a byte string -- b'something something'
		fileout.write(b"Oh that's line number ")
		# convert number into string so I can encode it into bytes
		fileout.write(str(i).encode('ASCII'))
		fileout.write(b'\n')

with gzip.open('file.txt.gz', 'r') as filein:
	for line in filein:
		print(line[:-1].decode('ASCII'))

s_in = b"Lots of content here"
s_out = gzip.compress(s_in)

print(s_out)

"""

import gzip

def get_dnaread(filename):
	### Extracts a DNA line from the next NGS read in the file

	flag = False

	# reads 4 lines i.e. the whole NGS read
	for j in range(4):
		# read line and decode from bytes
		line = filein.readline().decode('ASCII')

		# I thought we could use stateful parsing to make suuuuuure 
		# that we get the DNA line and not something else, or is it stupid?
			# other ideas for input control: regex
			# and of course add try-except
				# for example if flag is still True in the end 
				# of function -- something wrong

		if line[0] == '@':
			flag = True
		elif flag:
			if line[0] == '+':
				flag = False
			else:
				dna = line

	return dna
	# return a signal if reached end of file


### Main Program ###

filein = gzip.open('Unknown3_raw_reads_1.txt.gz', 'r')

dnaread = get_dnaread(filein)
# in the end switch the loop to 'while not end-of-file'
for i in range(5):					# get top 5 DNA entries
	print(dnaread)
	dnaread = get_dnaread(filein)
print(i)
filein.close()

# I tried 1 million DNA reads without printing them -- it takes 7 seconds heheh

