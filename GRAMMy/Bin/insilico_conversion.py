# Does insilico conversion ...
import sys
import os
import itertools

infile=sys.argv[1]
outfile=sys.argv[2]

conversion=sys.argv[3]

# If ever non-directional libraries are used ...
#if conversion == 'C2T':
#	incharacter='C'
#	outcharacter='T'
#if conversion == 'G2A':
#	incharacter='G'
#	outcharacter='A'

incharacter='C'
outcharacter='T'

with open(infile) as f, open(outfile, 'w') as w:
	for line in f:
		if line[0] == '>':
			w.write(line)
		else:
			w.write(line.replace(incharacter, outcharacter))
