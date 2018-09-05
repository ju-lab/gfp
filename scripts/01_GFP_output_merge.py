import sys
in_file=open('fusion_summary.txt')
ref_file=open('fusion_splits_homologyChecked.txt')
out_file=open('fusion_summary_splits.txt','w')


in_line=in_file.readline().strip()
pair_N_dic={}
pair_split_dic={}
while in_line:
	in_indi=in_line.split('\t')
	pairread_N=int(in_indi[2])
	if pairread_N > 3:  # set cutoff
		pair_N_dic[in_indi[0]+'\t'+in_indi[1]]=pairread_N
		pair_split_dic[in_indi[0]+'\t'+in_indi[1]]=''
	else:
		'blank'
	in_line=in_file.readline().strip()


ref_line=ref_file.readline().strip()
ref_list=[]
while ref_line:
	ref_indi=ref_line.split('\t')
	chr1=int(((ref_indi[2].split(':')[0][4:]).replace('X','23')).replace('Y','24'))
	pos1=int((ref_indi[2].split(':')[1]).split('-')[0])
	chr2=int(((ref_indi[4].split(':')[0][4:]).replace('X','23')).replace('Y','24'))
	pos2=int((ref_indi[4].split(':')[1]).split('-')[0])
	r1=ref_indi[2][1:];g1=ref_indi[1]
	r2=ref_indi[4][1:];g2=ref_indi[3]
	if g1+'\t'+g2 in pair_split_dic:
		if pair_split_dic[g1+'\t'+g2]=='':
			pair_split_dic[g1+'\t'+g2]=r1+'\t'+r2
		else:
			pair_split_dic[g1+'\t'+g2]=pair_split_dic[g1+'\t'+g2]+';'+r1+'\t'+r2
	elif g2+'\t'+g1 in pair_split_dic:
		if pair_split_dic[g2+'\t'+g1]=='':
			pair_split_dic[g2+'\t'+g1]=r2+'\t'+r1
		else:
			pair_split_dic[g2+'\t'+g1]=pair_split_dic[g2+'\t'+g1]+';'+r2+'\t'+r1
	else:  #split_read which is not exist in fusion_summary will be removed
		'blank'
	
	ref_line=ref_file.readline().strip()


for gene_fusion in pair_split_dic: 
	if pair_split_dic[gene_fusion]=='':   #gene_fusion_list w/o split read will be removed
		'blank'
	else:
		split_list=pair_split_dic[gene_fusion].split(';')
		for split in split_list:
			out_file.write(gene_fusion+'\t'+str(pair_N_dic[gene_fusion])+'\t'+split+'\n')
	


