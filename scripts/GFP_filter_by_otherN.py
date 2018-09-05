import sys

in_file=open(sys.argv[1])  #fusion_summary_split
ref_file=open(sys.argv[2]) 
out_file=open(sys.argv[1]+'.fi','w')
in_line=in_file.readline().strip()
ref_line=ref_file.readline().strip()

ref_list=[]

while ref_line:
	ref_indi=ref_line.split('\t')
	r1=ref_indi[0].split('-')[0][:-2]+'-'+ref_indi[0].split('-')[1][:-2]
	r2=ref_indi[1].split('-')[0][:-2]+'-'+ref_indi[1].split('-')[1][:-2]
	ref_list.append(r1+'\t'+r2)
	ref_line=ref_file.readline().strip()

print(len(ref_list))

n=0
while in_line:
	n=n+1
	in_indi=in_line.split('\t')
	r1=in_indi[3].split('-')[0][:-2]+'-'+in_indi[3].split('-')[1][:-2]
	r2=in_indi[4].split('-')[0][:-2]+'-'+in_indi[4].split('-')[1][:-2]
	if r1+'\t'+r2 in ref_list or r2+'\t'+r1 in ref_list:
		'blank'
	else:
		out_file.write(in_line+'\n')

	if n%1000==0:
		print(in_line)
	in_line=in_file.readline().strip()
