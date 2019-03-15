import sys
sample=sys.argv[1]

with open(sample+'.overlapping.tblat1.C2T') as c2t, open(sample+'.overlapping.tblat1.G2A') as g2a:
	for c2tline, g2aline in zip(c2t, g2a):
		c2t_evalue=c2tline.split('\t')[10]
		g2a_evalue=g2aline.split('\t')[10]
		with open('V1/C2T/exe/'+sample+'/'+sample+'.tblat.2', 'a') as c2t2, open('V1/G2A/exe/'+sample+'/'+sample+'.tblat.2', 'a') as g2a2:
			if c2t_evalue>=g2a_evalue:
				c2t2.write(c2tline)
			else:
				g2a2.write(g2aline)
