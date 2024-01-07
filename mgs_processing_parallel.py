import argparse as ap
import os
from pathlib import Path
import pandas as pd
import time
import multiprocessing as mp

# 1. Adapter trimming
def exec_trim(fwd_file, rev_file, output_dir, adapter, qual_scr, mismatches, minlength):
    fwd_filename = os.path.split(fwd_file)[-1]
    rev_filename = os.path.split(rev_file)[-1]
    try:
        os.makedirs(output_dir)
    except FileExistsError:
        pass
    fwd_paired = os.path.abspath(os.path.join(output_dir, fwd_filename))
    rev_paired = os.path.abspath(os.path.join(output_dir, rev_filename))
    fwd_unpaired = os.path.normpath(os.path.join(output_dir, 'unpaired_'+fwd_filename))
    rev_unpaired = os.path.normpath(os.path.join(output_dir, 'unpaired_'+rev_filename))
    
    if os.path.isfile(fwd_paired) and os.path.isfile(rev_paired):
        print('Trimmed files are already exists at the location. Skipping rerun trimming operation!')
    else:
        qual_scr = '-'+qual_scr
        cmd_trim = 'java -jar Tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE ' + qual_scr + ' ' + fwd_file + ' ' + rev_file + ' ' + \
            fwd_paired + ' ' + fwd_unpaired + ' ' + rev_paired + ' ' + rev_unpaired + \
            ' ILLUMINACLIP:' + adapter + ':' + str(mismatches) + ':30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:' + str(minlength) 
        print(cmd_trim)
        os.system(cmd_trim)
        os.remove(fwd_unpaired)
        os.remove(rev_unpaired)
    return fwd_paired, rev_paired

# 2. Decontamination
def remove_hostgenome(sampleid, fwd_file, rev_file, output_dir):
    fwdname = os.path.split(fwd_file)[-1]
    revname = os.path.split(rev_file)[-1]
    try:
        os.makedirs(output_dir)
    except FileExistsError:
        pass    
    
    if os.path.isfile('../database/hostDB/idx_hostgenomes.1.bt2'):
        print('Index of the host genome(s) are already exists at the location. To generate new index, delete the index file and rerun!')
    else:  
        cmd_merge = 'cat ../database/hostDB/*.fa > ../database/hostDB/hostgenomes.fna'
        os.system(cmd_merge)
        bowtie_idx = 'bowtie2-build ../database/hostDB/hostgenomes.fna ../database/hostDB/idx_hostgenomes'
        os.system(bowtie_idx)

    fwd_output = os.path.join(output_dir,fwdname)
    rev_output = os.path.join(output_dir,revname)
    if os.path.isfile(fwd_output) and os.path.isfile(rev_output):
        print('Files are already exists at the location. Skipping rerun removing host genome!')
    else:
        # bowtie2 mapping against host sequence database, keep both aligned and unaligned reads (paired-end reads)
        samfile = os.path.join(output_dir,sampleid+'_mapped_and_unmapped.sam')
        bowtie_cmd = 'bowtie2 -x ../database/hostDB/idx_hostgenomes -1 ' + fwd_file + ' -2 ' + rev_file + ' -S ' + samfile + ' --quiet'
        print(bowtie_cmd)
        os.system(bowtie_cmd)
        # convert file .sam to .bam
        bamfile = os.path.join(output_dir,sampleid+'_mapped_and_unmapped.bam')
        samcmd = 'Tools/samtools/bin/samtools view -bS ' + samfile  + ' > ' + bamfile  
        os.system(samcmd)
        # get unmapped pairs (both reads R1 and R2 unmapped)
        # -f  12    # Extract only ( -f ) alignments with both reads unmapped: <read unmapped><mate unmapped>
        # -F 256    # Do not (  -F  ) extract alignments which are: <not primary alignment>
        unmappair_bam = output_dir + '/' + sampleid + '_bothReadsUnmapped.bam'
        samcmd = 'Tools/samtools/bin/samtools view -b -f 12 -F 256 ' + bamfile + ' > ' + unmappair_bam
        os.system(samcmd)
        # sort bam file by read name ( -n ) to have paired reads next to each other
        sorted_bam = output_dir+ '/' + sampleid + '_bothReadsUnmapped_sorted.bam'
        samcmd = 'samtools sort -n  -f -o ' + unmappair_bam + ' ' + sorted_bam
        os.system(samcmd)
        # samcmd = 'Tools/samtools/bin/samtools fastq -@ 8 '  + ' -1 ' + fwd_output + ' -2 ' + rev_output + ' -0 /dev/null -s /dev/null -n ' + sorted_bam
        # os.system(samcmd)
        samcmd = 'bedtools bamtofastq -i ' + sorted_bam  + ' -fq ' + fwd_output + ' -fq2 ' + rev_output 
        print(samcmd)
        os.system(samcmd)
        try:
            os.remove(samfile)
            os.remove(bamfile)
            os.remove(unmappair_bam)
            os.remove(sorted_bam)
        except:
            print('FAILED')
    
    return fwd_output, rev_output

# 1. Quality control
def qc(sampleid,fwdfile_path,revfile_path,outdir):
    fwd = fwdfile_path
    rev = revfile_path
    # 1: trimming adapters 
    outdir_trim = os.path.join(outdir,'1.2_trimmed_adapter')
    adapter = 'Tools/Trimmomatic-0.39/adapters/TruSeq3-PE-Novogene.fa' 
    qual_scr = 'phred33'    # phred33 or phred64
    mismatches = 0  # Number between 0 to 100
    minlength = 1   # Number between 1 to 2000
    fwd, rev = exec_trim(fwd, rev, outdir_trim, adapter, qual_scr, mismatches, minlength)

    # 2: host genome filtering
    outdir_hostfilered = os.path.join(outdir,'1.3_filtered_hostgenome')
    start = time.time()
    fwd, rev = remove_hostgenome(sampleid,fwd, rev, outdir_hostfilered)

    return fwd, rev

def parallel_pipeline(sample):
    sampleid,fwdfile,revfile,outdir = sample[0],sample[1],sample[2],sample[3]
    print('Input: forward file {}, Reverse file {}'.format(fwdfile, revfile))

    basedir = '1_quality_control'    
    outdir_qc = os.path.join(outdir,basedir)
    fwdfile, revfile =  qc(sampleid,fwdfile,revfile,outdir_qc)
    print('Processed: forward file {}, Reverse file {}'.format(fwdfile, revfile))

    
parser = ap.ArgumentParser()
parser.add_argument('-o','--outdir',dest='output_dir', type=str, required=True, help='Output directory')
parser.add_argument('-m','--metafile',dest='metafile', type=str, required=True, 
                    help='An Excel file contains a list of samplename and the location of the corresponding fwd and rev file')

# input and output directories
args = parser.parse_args()
OUTDIR = args.output_dir
METAFILE = args.metafile

metadata = pd.read_csv(METAFILE, sep=' ')

# Using multiprocessing to execute parallel
samples = []
basedir = '1_quality_control'
outdir_qa = os.path.join(OUTDIR,basedir)
for sample_idx in metadata.index:
    sampleID = metadata.SampleID[sample_idx].strip()
    fwdfile = metadata.Forward_read[sample_idx].strip()
    revfile = metadata.Reverse_read[sample_idx].strip()
    samples.append([sampleID,fwdfile,revfile,OUTDIR])

pool = mp.Pool(mp.cpu_count())
result = pool.map(parallel_pipeline, samples)