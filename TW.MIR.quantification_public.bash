#!/bin/bash

###
### Pipeline for alignment and quantification of single-end MIR-sequencing data to human mature MIR database 
### Author: Terrence Wu
### Department: Genomic Medicine
### Organization: University of Texas, MD Anderson Cancer Center, Houston TX, USA
### 
### Date of first version:   August, 2018
### Date of current release: March,  2021
###
### Usage: TW.MIR.quantification_public.bash [-p 4]  input.fastq.gz  output_prefix  reference_bowtie2_index
###
### Input:
###      input.fastq.gz:            gzipped fastq from single-end MIR sequencing for alignment
###      output_prefix:             prefix of output files
###      reference_bowtie2_index:   reference MIR database from mirbase (fasta and bowtie2 indexes) 
###
### Options:
###      -p 4:  Use 4 threads for bowtie2 aligning (default: 4)
###
### Output:
###      {output_prefix}.MIR.bam            bam file from bowtie2
###      {output_prefix}.MIR.mapping.gz     list of aligned reads and matching 
###      {output_prefix}.MIR.counts.tsv     raw counts of reference MIRs
###      {output_prefix}.MIR.stats.tsv      basic statistics
###      
### Requirements:
###      samtools
###      bowtie2
###      gzip
###
### Test run with example data (expected run time < 1 min in a standard Desktop computer):
###      ./TW.MIR.quantification_public.bash  example.fastq.gz  example_output  hsa.mature.uniq.fa
###
### Licensed under the GNU General Public License v3.0
###
###

# module load bowtie2
# module load samtools

if [ "$1" == "-h" -o "$1" == "--help" ]; then
    cat "$0" | grep '^###'
    exit 0
fi

THREADS=4

if [ "$1" == "-p" -o "$1" == "--threads" ]; then
    THREADS=$2
    shift
    shift
fi

if [ -z "$(which bowtie2)" ]; then
    >&2 echo -e "\n## ERROR: bowtie2 is required \n"  
    exit 1
fi 

if [ -z "$(which samtools)" ]; then
    >&2 echo -e "\n## ERROR: samtools is required \n"  
    exit 1
fi 

if [ -z "$(which gzip)" ]; then
    >&2 echo -e "\n## ERROR: gzip is required \n"  
    exit 1
fi 

FQFN=$1

if [ -z "$FQFN" ]; then
    if [ -f "./example.fastq.gz" ]; then
        FQFN=example.fastq.gz
        >&2 echo -e "\n## Usage: TW.MIR.quantification_public.bash example.fastq.gz"
        >&2 echo -e "## Use example fastq as input \n"
    fi
fi

if [ ! -f "$FQFN" ]; then
    >&2 echo -e "\n## ERROR: Missing input FASTQ file $FQFN \n"  
    >&2 echo -e "## Usage: TW.MIR.quantification_public.bash example.fastq.gz \n"
    exit 1
fi

OFN=$2
if [ -z "$OFN" ]; then
    OFN=$(echo "$FQFN" | sed 's/.gz$//'  | sed 's/.bz2$//' | sed 's/.fq$//'  | sed 's/.fastq$//' | sed 's/$/_output/' )
fi

reference_fa=$3
if [ -z "$reference_fa" ]; then
    if [ -f "./hsa.mature.uniq.fa" ]; then
        reference_fa=$PWD/hsa.mature.uniq.fa
    fi
fi

if [ ! -f "$reference_fa" -a ! -f "$reference_fa.1.bt2" ]; then
    >&2 echo -e "\n## ERROR: Missing reference_bowtie2_index \n" 
    exit 1
fi

if [ ! -f "$reference_fa.1.bt2" ]; then
    
    bowtie2-build $reference_fa $reference_fa
    
    if [ ! -f "$reference_fa.1.bt2" ]; then
        >&2 echo -e "\n## ERROR: Could not build bowtie2 index for $reference_fa \n" 
        exit 1
    fi
fi

echo "Input FASTQ file = $FQFN"
echo "Output files = $OFN"
echo "Reference = $reference_fa"
echo "Number of Threads = $THREADS"

# reference_fa=/scratch/references/MIR/mirbase.v22_20180311/homo_sapiens/hsa.mature.uniq.fa

## Download from ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
## Current project used mirbase v22 (2018-03-11)
## Retrieve reference mature MIRs of Homo sapiens
## Convert U to T
## Then combine MIR/MIR families with identical sequences
## e.g. hsa-miR-548aa_MIMAT001844 and hsa-miR-548t-3p_MIMAT0022730 (AAAAACCACAATTACTTTTGCACCA)
## Combined fa file like this: (__n25 is the length of the MIR sequence, this number must be in the name)
## ...
## >hsa-miR-548aa_MIMAT0018447,hsa-miR-548t-3p_MIMAT0022730__n25 Homo sapiens 
## AAAAACCACAATTACTTTTGCACCA
## >hsa-miR-548ab_MIMAT0018928__n22 Homo sapiens 
## AAAAGTAATTGTGGATTTTGCT
## ...

MIR_bam="$OFN.MIR.bam"
MIR_map="$OFN.MIR.mapping.gz"
MIR_raw_counts="$OFN.MIR.counts.tsv"
MIR_stats="$OFN.MIR.stats.tsv"

bowtie2 -p $THREADS -N 1 -L 16 -k 5 --local --norc   -x $reference_fa -U "$FQFN" | samtools view -bh - > $MIR_bam
##  -p 4        number of alignment threads to launch = 4
##  -N 1        max 1 mismatches in seed alignment
##  -L 16       length of seed substrings = 16
##  -k 5        report up to 5 alns per read; MAPQ not meaningful
##  --local     local alignment; ends might be soft clipped 
##  --norc      do not align reverse-complement version of read 

samtools view $MIR_bam | 
    gawk -v FS="\t" '$6~/^[1-9][0-9]*M/' |  ### No clipping at 5' end
    gawk '''
            ### For each MIR read, filter out alignments (in case of multiple alignment) 
            ###  whose AS score were lower than the best one 
            ###  (highest AS, best alignment, always occurred first)
            
            BEGIN {
                FS="\t"
                nowAS=""
                nowID=""
            }
            /^@/ {
                print $0
            }
            !/^@/ {
                newID=$1
                newAS=$0
                if (newAS~/\tAS:i:[0-9]/) {
                    sub("^.*\tAS:i:", "", newAS)
                    sub("[^0-9].*$", "", newAS)
                    newAS=strtonum(newAS)
                } else {
                    newAS=0
                }
                if (newID!=nowID || NR==1) {
                    nowID=newID
                    nowAS=newAS
                    print $0
                } else {
                    if (newAS>=nowAS) print $0
                }
            }
         '''  | 
    gawk '''
            BEGIN {
                FS="\t"
                OFS="\t"
            }
            {
                nMatch_5prime=$6
                sub("M.*$", "", nMatch_5prime); 
                
                ref_MIR_size=$3;
                sub(".*_n", "", ref_MIR_size); 
                
                if (ref_MIR_size==nMatch_5prime) {
                    print $1, $3, $6
                }
            }
         ''' |  ### filter out alignments that map only to part of a reference MIR (number of matches < length of reference MIR)
    gawk '''
            ## concatenate field2 with semicomma if current field1 is equal to previous field1
            BEGIN {
                FS="\t"
                OFS="\t"
                f1="";
                f2="";
            }
            END {
                print f1, f2;
            }
            {
                if ($1!=f1) { 
                    if (NR>1) print f1, f2; 
                    f1=$1;
                    f2=$2;
                } else { 
                    f2=f2";"$2;
                }
            }
         ''' | 
    gzip > $MIR_map

gzip -d -c $MIR_map | 
    gawk -v FS="\t" '$2!~/;/' |  ## get rid of reads that mapped to more than one reference MIRs (best AS scores tied)
    cut -f 2 | 
    sort | 
    gawk '''
            ## count occurrence of each reference MIR id
            BEGIN {
                nowseq="";
                nowoccur=0;
                FS="\t"
                OFS="\t"
            }
            END {
                print nowseq, nowoccur;
            }
            {
                if ($1!=nowseq) { 
                    if (NR>1) print nowseq, nowoccur; 
                    nowoccur=1;
                    nowseq=$1;
                } else { 
                    nowoccur=nowoccur+1;
                }
            }
         ''' | 
    sed 's/__/\t/' > $MIR_raw_counts

echo -e "Input fastq:\t$FQFN" > $MIR_stats
echo -e "Output prefix:\t$OFN" >> $MIR_stats
samtools view $MIR_bam | cut -f 1 | uniq | wc -l | sed 's/ //g' | sed 's/^/Total reads:\t/' >> $MIR_stats
gzip -d -c $MIR_map | grep ';' | wc -l | sed 's/ //g' | sed 's/^/Ambiguously aligned reads(not counted):\t/' >> $MIR_stats
gzip -d -c $MIR_map | grep -v ';' | wc -l | sed 's/ //g' | sed 's/^/Uniquely aligned reads:\t/' >> $MIR_stats

echo -e "\n"
cat $MIR_stats

exit 0

