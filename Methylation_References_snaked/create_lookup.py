#!/usr/bin/env python3

"""
Create lookup table to associate file with new name and tissue group
Alexandre Pellan Cheng 2018
"""
import os
import sys
from fnmatch import fnmatch
# List tissue groups for this project
tissue_groups = ["macrophage", "BCell", "bladder","TCell","monocyte","NKCell", \
					"dendritic","eosonophil","erythroblast", "liver/hepatocyte", \
                    "pancreas/islet","kidney/podocyte/mesangial","colon/intestine",\
                    "neutrophil","skin","spleen"]

# List all methylation references
def find_references(desired_path, ext1, ext2, ext3, ext4):
	references=[]
	for path, subdirs, files in os.walk(desired_path):
			for name in files:
				if fnmatch(name, ext1) or fnmatch(name, ext2) or fnmatch(name, ext3) or fnmatch(name, ext4):
					references.append(os.path.join(path, name)[2:])
	return(references)

def create_tissue_groups(tissue_groups):
	groups={}
	i=1
	for tissue in tissue_groups:
		groups[tissue]="G"+str(i)
		i=i+1
	return(groups)

def create_lookup(groups, references, tissue_groups):
	reference_list=[]
	group=[]
	group_id=[]
	new_name=[]
	for file in references:
		reference_list.append(file)
		for tissue in tissue_groups:
			acceptable_tissues = tissue.split('/')
			for t in acceptable_tissues:
				if t in file:
					tissue_group = groups[tissue]
					group.append(tissue_group)
					group_id.append(groups[tissue]+'_'+str(group.count(tissue_group)))
					new_name.append(acceptable_tissues[0]+str(group.count(tissue_group)))
					break
	return([reference_list, group, group_id, new_name])

def main():
	outfile=sys.argv[1]
	references=find_references('./', '*.bed.gz', '*.bw', '*.bigWig', '*.bam')
	groups = create_tissue_groups(tissue_groups)
	final=create_lookup(groups, references, tissue_groups)
	ref=final[0]
	group=final[1]
	group_id=final[2]
	new_name=final[3]
	lines=[]
	for r, g, gid, nn in zip(ref, group, group_id, new_name):
		line=[g,gid,nn, r]
		lines.append('\t'.join(line))
	with open(outfile, 'w') as w:
		w.write('\n'.join(lines))
	os.system('sort -V ' + outfile + ' > tmp && mv tmp ' + outfile)

if __name__ == '__main__':
	main()
