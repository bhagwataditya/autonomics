# Download - install - configure SRA toolkit
#===========================================
# https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
# https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
# https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration
# https://github.com/ncbi/sra-tools/wiki/HowTo:-Access-SRA-Data

# Download Stem Cell Comparison Data
#===================================
prefetch SRR1909613   # BM +
prefetch SRR1909637   #    +
prefetch SRR1909638   #    + 
prefetch SRR1909639   #    +
prefetch SRR1926136   # EM +
prefetch SRR1926132   #    +
prefetch SRR1926133   #    +
prefetch SRR1926134   # E  +
prefetch SRR1926135   #    +
prefetch SRR1867792   #    +

# Convert SRA -> FASTQ
#=====================
# https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
fasterq-dump SRR1909613   # BM  17 M
fasterq-dump SRR1909637   #     16 M
fasterq-dump SRR1909638   #     28 M
fasterq-dump SRR1909639   #     25 M
fasterq-dump SRR1926136   # EM  31 M
fatserq-dump SRR1926132   #     32 M
fasterq-dump SRR1926133   #     16 M
fasterq-dump SRR1926134   # E    4 M
fasterq-dump SRR1926135   #     15 M
fasterq-dump SRR1867792   #     32 M

rm -r SRR1909613   # BM
rm -r SRR1909637   #   
rm -r SRR1909638   #   
rm -r SRR1909639   #   
rm -r SRR1926136   # EM
rm -r SRR1926132   #   
rm -r SRR1926133   #   
rm -r SRR1926134   # E 
rm -r SRR1926135   #   
rm -r SRR1867792   #   
