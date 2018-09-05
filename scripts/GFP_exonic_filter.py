import sys
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.fi','w')
in_line=in_file.readline().strip()
while in_line:
	in_indi=in_line.split('\t')
	if in_indi[5][0:4]=='exon' and in_indi[7][0:4]=='exon':
		out_file.write(in_line+'\n')
	else:
		'blank'
	in_line=in_file.readline().strip()

	
