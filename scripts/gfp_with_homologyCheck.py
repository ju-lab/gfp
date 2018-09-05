#!/usr/bin/python

### script by Won-Chul Lee

import sys, getopt, os, string, gzip


def main():

    picard= "java -jar /home/users/tools/picard/dist/picard.jar"
    bl2seqPATH= "/home/users/tools/blast-2.2.24/bin/bl2seq"

    sampleID= "ND"
    read1_fastq, read2_fastq= '', ''
    trim_len= 25
    read_mm= 1  # mismatches allowed for each of discordant read-pairs
    ref_fasta= ''
    gene_bed= ''
    tx_fasta= ''
    output_root_dir= ''
    num_threads= 1

    if len(sys.argv) == 1:
        print "\nUsage:"
        print "\t--sample                sample ID.\n"
        print "\t--fastq1                path to read1 fastq (gzipped allowed) file.\n"
        print "\t--fastq2                path to read2 fastq (gzipped allowed) file.\n"
        print "\t--ref-fa                path to reference fasta file.\n"
        print "\t--gene-bed              exon/intron information bed file.\n"
        print "\t--tx-fa                 path to transcript fasta file.\n"
        print "\t--out-dir               path to analysis output root directory.\n"
        print "\t--trim-len              read trimming length. (default: %d)\n" %trim_len
        print "\t--read-mm               mismathces allowed for each of discordant read-pairs. (default: %d)\n" %read_mm
        print "\t--num-threads           number of cpu threads. (default: %d)\n" %num_threads
        sys.exit(1)
    opts, args= getopt.getopt(sys.argv[1:], "", ["sample=", "fastq1=", "fastq2=", "ref-fa=", "gene-bed=", "tx-fa=", "out-dir=", "trim-len=", "read-mm=", "num-threads="])
    for opt, arg in opts:
        if opt == "--sample": sampleID= arg
        elif opt == "--fastq1": read1_fastq= os.path.abspath(arg)
        elif opt == "--fastq2": read2_fastq= os.path.abspath(arg)
        elif opt == "--ref-fa": ref_fasta= os.path.abspath(arg)
        elif opt == "--gene-bed": gene_bed= os.path.abspath(arg)
        elif opt == "--tx-fa": tx_fasta= os.path.abspath(arg)
        elif opt == "--out-dir": output_root_dir= os.path.abspath(arg)
        elif opt == "--trim-len": trim_len= int(arg)
        elif opt == "--read-mm": read_mm= int(arg)
        elif opt == "--num-threads": num_threads= int(arg)

    #########################################
    # Make output root directory if not exist
    #########################################
    try:
        os.system("mkdir "+output_root_dir)
        os.chdir(output_root_dir)
    except:
        os.chdir(output_root_dir)

    ##
    # Adapter trimming using cutadapt
    ##
    os.system("cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o read1_adapterTrimmed.fastq.gz -p read2_adapterTrimmed.fastq.gz --minimum-length 20 %s %s" %(read1_fastq, read2_fastq))

    ##
    # Align original reads
    ##
    os.system("bwa mem -t %d %s read1_adapterTrimmed.fastq.gz read2_adapterTrimmed.fastq.gz > origin.sam" %(num_threads, ref_fasta))
    os.system("%s SortSam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=picard_tmp_dir I=origin.sam O=origin.sorted.sam" %(picard))   # sort reads
    os.system("%s MarkDuplicates M=picard_MarkDuplicates_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT TMP_DIR=picard_tmp_dir I=origin.sorted.sam O=origin.markdup.sam" %(picard))   # mark duplicates
    os.system("samtools view -bS origin.markdup.sam > origin.markdup.bam")
    os.system("samtools index origin.markdup.bam")
    os.system("rm origin.sam origin.sorted.sam picard_MarkDuplicates_metrics.txt")

    #####################
    # Align trimmed reads
    #####################
    trim_fastq("read1_adapterTrimmed.fastq.gz", output_root_dir, trim_len, "read1_adapterTrimmed_3pTrimmed")    # trim read1 fastq
    trim_fastq("read2_adapterTrimmed.fastq.gz", output_root_dir, trim_len, "read2_adapterTrimmed_3pTrimmed")    # trim read2 fastq
    os.system("bwa mem -T %d -t %d %s read1_adapterTrimmed_3pTrimmed.fastq read2_adapterTrimmed_3pTrimmed.fastq > trimmed.sam" %(trim_len-read_mm, num_threads, ref_fasta))
    os.system("%s SortSam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=picard_tmp_dir I=trimmed.sam O=trimmed.sorted.sam" %(picard)) # sort reads
    os.system("%s MarkDuplicates M=picard_MarkDuplicates_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=picard_tmp_dir I=trimmed.sorted.sam O=trimmed.dedup.sam" %(picard))      # remove duplicates
    os.system("rm read?_adapterTrimmed_3pTrimmed.fastq trimmed.sam trimmed.sorted.sam picard_MarkDuplicates_metrics.txt")
    os.system("%s CollectInsertSizeMetrics H=insert_size.pdf I=trimmed.dedup.sam O=insert_size.stats" %(picard))
    ins_mean, ins_std= 0, 0
    ins_size_fi= open("insert_size.stats", 'r')
    while 1:
        line= ins_size_fi.readline()
        if not line: break
        if line.startswith("MEDIAN_INSERT_SIZE"):
            f= ins_size_fi.readline().rstrip().split("\t")
            ins_mean= float(f[4])
            ins_std= float(f[5])
            break

    #################################
    # Queryname sort & convert to bed
    #################################
    os.system("%s SortSam SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT TMP_DIR=picard_tmp_dir I=trimmed.dedup.sam O=trimmed.sortByName.sam" %(picard))
    os.system("rm -rf picard_tmp_dir")
    sam_fi= open("trimmed.sortByName.sam", 'r')
    bed_fo= open("discordant_reads.bed", 'w')
    while 1:
        line= sam_fi.readline()
        if not line: break
        if line.startswith('@'): continue
        fields= line.rstrip().split("\t")
        if fields[2] == '*': continue   # unmapped reads

        MAPQ= int(fields[4])
        if MAPQ == 0: continue  # only consider uniquely mapped reads

        TLEN= abs(int(fields[8]))
        TLEN_z= float(TLEN-ins_mean)/ins_std
        if -1.96 < TLEN_z < 1.96: continue  # skip reads in a concordant pair

        FLAG= int(fields[1])
        ref_cur_pos= int(fields[3]) # reference current position
        CIGAR= fields[5]
        CIGAR_cur_pos= 0
        ref_coords= []
        for i in range(len(CIGAR)):
            char= CIGAR[i:i+1]
            if char in ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']:
                num= int(CIGAR[CIGAR_cur_pos:i])
                CIGAR_cur_pos= i+1
                if char == 'M':
                    ref_coords.append(ref_cur_pos)
                    ref_coords.append(ref_cur_pos+num-1)
                    ref_cur_pos= ref_cur_pos+num
                elif char in ['S', 'H', 'I']: continue
                elif char == 'D': ref_cur_pos= ref_cur_pos+num
        ref_coords.sort()

        mismatch= "NA"
        for field in fields:
            if field.startswith("NM:"): mismatch= field.split(':')[-1]

        read_ID= "read1"
        if FLAG & 64 != 64: read_ID= "read2"
        if len(ref_coords) != 0:
            bed_fo.write("%s\t%d\t%d\t%s\t%s\t%s\t%d\n" %(fields[2], ref_coords[0]-1, ref_coords[-1], fields[0], read_ID, mismatch, TLEN))
    bed_fo.close()
    sam_fi.close()

    ###########################
    # Annotate gene information
    ###########################
    os.system("bedtools intersect -a discordant_reads.bed -b %s -wa -wb > discordant_reads.anno.bed" %(gene_bed))
    os.system("rm discordant_reads.bed")

    ####################################
    # Find read-pairs supporting fusions
    ####################################
    bed_fi= open("discordant_reads.anno.bed", 'r')
    fo= open("fusion_read_pairs.txt", 'w')
    cur_name= ''
    read_pair_lines= []
    while 1:
        line= bed_fi.readline()
        if not line: break
        fields= line.rstrip().split("\t")
        if fields[3] != cur_name:
            if cur_name != '':
                find_fusion_pair(read_pair_lines, fo)
            cur_name= fields[3]
            read_pair_lines= [line]
        else: read_pair_lines.append(line)
    find_fusion_pair(read_pair_lines, fo)
    fo.close()
    bed_fi.close()

    #####################################
    # Find split reads supporting fusions
    #####################################
    extract_clipped_reads("origin.markdup.sam", "soft_clipped_regions.fa")
    os.system("/home/users/tools/blat -out=blast8 -minIdentity=95 %s soft_clipped_regions.fa soft_clipped_regions.out" %(ref_fasta))
    blat2bed("soft_clipped_regions.out", "soft_clipped_regions.bed")
    os.system("bedtools intersect -a soft_clipped_regions.bed -b %s -wa -wb > soft_clipped_regions.anno.bed" %(gene_bed))
    os.system("rm soft_clipped_regions.bed")

    fi= open("soft_clipped_regions.anno.bed")
    fo= open("fusion_splits.txt", 'w')
    cur_name= ''
    sub_lines= []
    while 1:
        line= fi.readline()
        if not line: break
        fields= line.rstrip().split("\t")
        if fields[4] != cur_name:
            if cur_name != '':
                find_fusion_splits(sub_lines, fo)
            cur_name= fields[4]
            sub_lines= [line]
        else: sub_lines.append(line)
    find_fusion_splits(sub_lines, fo)
    fo.close()
    fi.close()

    ##
    # Check sequence homology
    ##
    tx_fasta_dict= {}
    build_tx_fasta_dict(tx_fasta, tx_fasta_dict)
    homologous_tx_pairs= {}

    fi= open("fusion_read_pairs.txt", 'r')
    fo= open("fusion_read_pairs_homologyChecked.txt", 'w')
    while 1:
        line= fi.readline()
        if not line: break
        f= line.rstrip().split("\t")
        tx_pair= [f[-2], f[-1]]
        tx_pair.sort()
        if tx_pair[0] + "\t" + tx_pair[1] in homologous_tx_pairs: continue
        is_homologous= check_homology(tx_pair[0], tx_pair[1], bl2seqPATH, tx_fasta_dict)
        if is_homologous:
            homologous_tx_pairs[tx_pair[0] + "\t" + tx_pair[1]]= ''
            continue
        fo.write(line)
    fo.close()
    fi.close()

    fi= open("fusion_splits.txt", 'r')
    fo= open("fusion_splits_homologyChecked.txt", 'w')
    while 1:
        line= fi.readline()
        if not line: break
        f= line.rstrip().split("\t")
        tx_pair= [f[-2], f[-1]]
        tx_pair.sort()
        if tx_pair[0] + "\t" + tx_pair[1] in homologous_tx_pairs: continue
        is_homologous= check_homology(tx_pair[0], tx_pair[1], bl2seqPATH, tx_fasta_dict)
        if is_homologous:
            homologous_tx_pairs[tx_pair[0] + "\t" + tx_pair[1]]= ''
            continue
        fo.write(line)
    fo.close()
    fi.close()

    ##
    # Finalize fusion discovery
    ##
    pairs= open("fusion_read_pairs_homologyChecked.txt").readlines()
    splits= open("fusion_splits_homologyChecked.txt").readlines()

    summary= {}

    for pair in pairs:
        f= pair.rstrip().split("\t")
        pair_id= f[0]
        genes= [f[1], f[2]]
        genes.sort()
        fusion= genes[0]+"\t"+genes[1]
        if fusion in summary:
            if not pair_id in summary[fusion]: summary[fusion].append(pair_id)
        else: summary[fusion]= [pair_id]

    for split in splits:
        f= split.rstrip().split("\t")
        pair_id= f[0]
        genes= [f[1], f[3]]
        genes.sort()
        fusion= genes[0]+"\t"+genes[1]
        if fusion in summary:
            if not pair_id in summary[fusion]: summary[fusion].append(pair_id)
        else: summary[fusion]= [pair_id]

    fo_summary= open("fusion_summary.txt", 'w')
    for fusion in summary:
        fo_summary.write(fusion+"\t"+str(len(summary[fusion]))+"\n")
    fo_summary.close()

    # cjy added 2018.08.10
    os.system('rm origin.markdup.sam read1_adapterTrimmed.fastq.gz read2_adapterTrimmed.fastq.gz trimmed.dedup.sam trimmed.sortByName.sam')

def check_homology(tx1, tx2, bl2seqPATH, tx_fasta_dict):
    is_homologous= False

    if not tx1 in tx_fasta_dict: return is_homologous
    if not tx2 in tx_fasta_dict: return is_homologous

    open("tx1.fasta", 'w').write(">%s\n%s" %(tx1, tx_fasta_dict[tx1]))
    open("tx2.fasta", 'w').write(">%s\n%s" %(tx2, tx_fasta_dict[tx2]))
    os.system(bl2seqPATH+" -p blastn -e 0.01 -D 1 -i tx1.fasta -j tx2.fasta -o temp.bl2seqout")
#   if len(open("temp.bl2seqout").readlines()) > 3:
#       is_homologous= True
    try: #try, except added, ysj 20160504
        if len(open("temp.bl2seqout").readlines()) > 3:
            is_homologous= True
    except:
        is_homologous = is_homologous

    os.system("rm tx1.fasta tx2.fasta temp.bl2seqout")

    return is_homologous


def build_tx_fasta_dict(tx_fasta, tx_fasta_dict):
    name, seq= '', ''
    fi= open(tx_fasta, 'r')
    while 1:
        line= fi.readline()
        if not line:
            tx_fasta_dict[name]= seq
            break
        if line.startswith('>'):
            if name != '': tx_fasta_dict[name]= seq
            name= line.rstrip()[1:]
            seq= ''
        else: seq= seq + line.rstrip()


def find_fusion_splits(lines, fo):

    left_genes, right_genes= [], []
    left_transcripts, right_transcripts= [], []
    left_interval, right_interval= '', ''
    read_name= ''

    for line in lines:
        f= line.rstrip().split("\t")
        interval= "%s%s:%s-%s" %(f[3], f[0], f[1], f[2])
        read_name= f[4]

        if f[5] == "left":
            left_interval= interval
            for anno in f[-1].split('|'):
                gene= string.join(anno.split('.')[1:-1], '.')
                transcript= anno.split('.')[0]
                if not gene in left_genes:
                    left_genes.append(gene)
                    left_transcripts.append(transcript)
        else:   # right
            right_interval= interval
            for anno in f[-1].split('|'):
                gene= string.join(anno.split('.')[1:-1], '.')
                transcript= anno.split('.')[0]
                if not gene in right_genes:
                    right_genes.append(gene)
                    right_transcripts.append(transcript)

    is_fusion= True
    for left_gene in left_genes:
        if left_gene in right_genes: is_fusion= False
    if len(left_genes) == 0 or len(right_genes) == 0: is_fusion= False

    if is_fusion:
        for i in range(len(left_genes)):
            for j in range(len(right_genes)):
                fo.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(read_name, left_genes[i], left_interval, right_genes[j], right_interval, left_transcripts[i], right_transcripts[j]))


def blat2bed(infile, outfile):
    name_dict= {}
    fi= open(infile, 'r')
    while 1:
        line= fi.readline()
        if not line: break
        fields= line.rstrip().split("\t")
        if not fields[0] in name_dict: name_dict[fields[0]]= 1
        else: name_dict[fields[0]]+= 1
    fi.close()
    fi= open(infile, 'r')
    fo= open(outfile, 'w')
    while 1:
        line= fi.readline()
        if not line: break
        fields= line.rstrip().split("\t")
        if name_dict[fields[0]] != 1: continue  # skip queries mapped to multiple regions

        pident= float(fields[2])
        if pident < 95: continue    # % identity cutoff

        read_name= fields[0][:fields[0].rfind('_')]     
        if fields[0].endswith('|'): # right clipped
            interval= fields[0][fields[0].rfind('_')+1:][:-1]
            fo.write("%s\t%s\t%s\t%s\t%s\tleft\n" %(interval[1:].split(':')[0], interval.split(':')[1].split('-')[0], interval.split(':')[1].split('-')[1], interval[:1], read_name))
            temp= [int(fields[8]), int(fields[9])]
            strand= '+'
            if temp[0] > temp[1]: strand= '-'
            temp.sort()
            fo.write("%s\t%d\t%d\t%s\t%s\tright\n" %(fields[1], temp[0], temp[1], strand, read_name))
        else:
            interval= fields[0][fields[0].rfind('|')+1:]
            temp= [int(fields[8]), int(fields[9])]
            strand= '+'
            if temp[0] > temp[1]: strand= '-'
            temp.sort()
            fo.write("%s\t%d\t%d\t%s\t%s\tleft\n" %(fields[1], temp[0], temp[1], strand, read_name))
            fo.write("%s\t%s\t%s\t%s\t%s\tright\n" %(interval[1:].split(':')[0], interval.split(':')[1].split('-')[0], interval.split(':')[1].split('-')[1], interval[:1], read_name))
    fo.close()
    fi.close()


def extract_clipped_reads(infile, outfile):
    fi= open(infile, 'r')
    fo= open(outfile, 'w')
    while 1:
        line= fi.readline()
        if not line: break
        if line.startswith('@'): continue
        fields= line.rstrip().split("\t")

        FLAG= int(fields[1])

        if FLAG & 1024 == 1024: continue    # skip duplicated reads

        QNAME= fields[0]+":read2"
        if FLAG & 64 == 64: QNAME= QNAME.replace("read2", "read1")

        RNAME= fields[2]
        ref_cur_pos= int(fields[3])
        CIGAR= fields[5]
        SEQ= fields[9]
        strand= '+'
        if FLAG & 16 == 16: strand= '-'

        if CIGAR.find('I') != -1 or CIGAR.find('D') != -1: continue # skip reads with insertion or deletion
        if CIGAR.count('M') != 1: continue  # skip reads mapped to multiple regions of reference
        if CIGAR.count('S') != 1: continue  # skip reads without only one soft clipped region
        MAPQ= int(fields[4])
        if MAPQ == 0: continue  # only consider uniquely mapped reads

        CIGAR_cur_pos= 0
        read_cur_pos= 0
        ref_mapped_interval= ''
        soft_clipped_seq= ''
        for i in range(len(CIGAR)):
            char= CIGAR[i:i+1]
            if char in ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']:
                num= int(CIGAR[CIGAR_cur_pos:i])
                CIGAR_cur_pos= i+1
                if char == 'M':
                    if soft_clipped_seq == '': ref_mapped_interval= strand+RNAME+':'+str(ref_cur_pos)+'-'+str(ref_cur_pos+num-1)+'|'
                    else: ref_mapped_interval= '|'+strand+RNAME+':'+str(ref_cur_pos)+'-'+str(ref_cur_pos+num-1)
                    ref_cur_pos= ref_cur_pos+num
                    read_cur_pos= read_cur_pos+num
                elif char == 'S':
                    soft_clipped_seq= SEQ[read_cur_pos:read_cur_pos+num]
                    read_cur_pos= read_cur_pos+num
        if len(soft_clipped_seq) < 20: continue # the recommended minimum length by blat
        fo.write(">%s_%s\n%s\n" %(QNAME, ref_mapped_interval, soft_clipped_seq))
    fo.close()
    fi.close()


def find_fusion_pair(lines, fo):
    num_read1, num_read2= 0, 0
    read1_genes, read2_genes= [], []
    read1_transcripts, read2_transcripts= [], []
    for line in lines:
        f= line.rstrip().split("\t")
        if f[4] == "read1":
            num_read1+= 1
            for anno in f[-1].split('|'):
                gene= string.join(anno.split('.')[1:-1], '.')
                transcript= anno.split('.')[0]
                if not gene in read1_genes:
                    read1_genes.append(gene)
                    read1_transcripts.append(transcript)
        elif f[4] == "read2":
            num_read2+= 1
            for anno in f[-1].split('|'):
                gene= string.join(anno.split('.')[1:-1], '.')
                transcript= anno.split('.')[0]
                if not gene in read2_genes:
                    read2_genes.append(gene)
                    read2_transcripts.append(transcript)

    if num_read1 != 1 or num_read2 != 1: return # both reads are uniquely mapped

    is_fusion= True
    for read1_gene in read1_genes:
        for read2_gene in read2_genes:
            if read1_gene == read2_gene: is_fusion= False

    if is_fusion:
        for i in range(len(read1_genes)):
            for j in range(len(read2_genes)):
                fo.write("%s\t%s\t%s\t%s\t%s\n" %(f[3], read1_genes[i], read2_genes[j], read1_transcripts[i], read2_transcripts[j]))


def trim_fastq(fastq, out_dir, trim_len, out_prefix):
    fi= ''
    if fastq.endswith(".gz"): fi= gzip.open(fastq, "rb")
    else: fi= open(fastq, 'r')
    fo= open(os.path.join(out_dir, out_prefix+".fastq"), 'w')
    while 1:
        line= fi.readline()
        if not line: break
        fo.write(line)
        line= fi.readline()
        if len(line.rstrip()) < trim_len or trim_len == 0:
            fo.write(line)
            fo.write(fi.readline())
            fo.write(fi.readline())
        else:
            fo.write(line.rstrip()[:trim_len]+"\n")
            fo.write(fi.readline())
            fo.write(fi.readline().rstrip()[:trim_len]+"\n")
    fo.close()
    fi.close()
    

if __name__ == "__main__":
    main()
