import sys, os

fn_file=os.popen("ls | grep "+sys.argv[1]) # input: grep item
fn_line=fn_file.readline().rstrip()
os.system('mkdir annovar_intermediate')
#print dbsnp_dic
while fn_line:
	in_file=open(fn_line)
	out_file=open('./annovar_intermediate/'+fn_line+'.anv_input','w')
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0]=='#':
			'blank'
		else:
			in_indi=in_line.split('\t')
			chr1=in_indi[3].split(':')[0][3:]
			pos1=(int((in_indi[3].split(':')[1]).split('-')[0])+int(in_indi[3].split('-')[1]))//2
			chr2=in_indi[4].split(':')[0][3:]
			pos2=(int((in_indi[4].split(':')[1]).split('-')[0])+int(in_indi[4].split('-')[1]))//2
			out_file.write(chr1+'\t'+str(pos1)+'\t'+str(pos1)+'\tA\tC\t0\n')
			out_file.write(chr2+'\t'+str(pos2)+'\t'+str(pos2)+'\tA\tC\t0\n')
		in_line=in_file.readline().strip()
	out_file.close()
	vcf_fn='./annovar_intermediate/'+fn_line+'.anv_input'
	out_format='./annovar_intermediate/'+fn_line
	
	annovar_dir="/home/users/tools/annovar/annovar/"
	db_dir="/home/users/data/02_annotation/02_annovar/humandb_hg19"
	

	os.system("%sannotate_variation.pl --geneanno -out %s -buildver hg19 %s %s" % (annovar_dir, out_format, vcf_fn, db_dir))

	var_dic={};exonic_dic={}
	in_file.close()
	in_file=open('./annovar_intermediate/'+fn_line+'.variant_function')
	in_line=in_file.readline().strip()
	while in_line:	
		in_indi=in_line.split('\t')
		var_dic[in_indi[2]+'\t'+in_indi[3]]=in_indi[0]+'\t'+in_indi[1]
		in_line=in_file.readline().strip()
	in_file.close()

	in_file=open(fn_line)
	out_file=open(fn_line+'.anv','w')
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0:5]=='#CHR1':
			out_file.write(in_line+'\tANV_loca1\tANV_gene1\tANV_loca2\tANV_gene2\n')
		else:
			in_indi=in_line.split('\t')
			chr1=in_indi[3].split(':')[0][3:]
			pos1=(int((in_indi[3].split(':')[1]).split('-')[0])+int(in_indi[3].split('-')[1]))//2
			chr2=in_indi[4].split(':')[0][3:]
			pos2=(int((in_indi[4].split(':')[1]).split('-')[0])+int(in_indi[4].split('-')[1]))//2
			idx1=chr1+'\t'+str(pos1)
			idx2=chr2+'\t'+str(pos2)
			if chr1==chr2:
				dist=str(abs(pos1-pos2)/float(1000))+'k'
			else:
				dist='TRA'
			out_file.write(in_line+'\t'+var_dic[idx1]+'\t'+var_dic[idx2]+'\t'+str(dist)+'\n')
		in_line=in_file.readline().strip()

	fn_line=fn_file.readline().strip()
