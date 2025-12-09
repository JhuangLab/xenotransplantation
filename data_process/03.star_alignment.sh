# star index monkey ref
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/mmul_index \
     --genomeFastaFiles /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/GCF_003339765.1_Mmul_10_genomic.fna \
     --sjdbGTFfile /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/genomic.gtf \
     --sjdbOverhang 99

# star index pig ref
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/sscrofa_index \
     --genomeFastaFiles /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.fna \
     --sjdbGTFfile /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/genomic.gtf \
     --sjdbOverhang 99

samples=(LC LT RC RT)

# Reference STAR index paths
monkey_index="/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/mmul_index"
pig_index="/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/sscrofa_index"

# Output directory base
out_base="./results/star"

# Loop through each sample and each species
for sample in "${samples[@]}"; do
    for species in monkey pig; do
        # Choose correct STAR index
        if [ "$species" == "monkey" ]; then
            ref="$monkey_index"
        else
            ref="$pig_index"
        fi

        # Make output directory
        out_dir="${out_base}/${species}"
        mkdir -p "$out_dir"

        # STAR alignment
        echo "STAR --runThreadN 8 --genomeDir $ref --readFilesIn /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/${sample}_1_R1.fq.gz /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/${sample}_1_R2.fq.gz --readFilesCommand zcat --outFileNamePrefix ${out_dir}/${sample}_ --outSAMtype BAM SortedByCoordinate"
    done
done > star_alignment_pipe.sh
job -s star_alignment_pipe.sh -c 8


# ======================== not to run =======================
# align to monkey ref
STAR --runThreadN 8 \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/mmul_index \
     --readFilesIn /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LC_1_R1.fq.gz /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LC_1_R2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/star_monkey/LC_ \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 8 \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/mmul_index \
     --readFilesIn /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LT_1_R1.fq.gz /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LT_1_R2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/star_monkey/LT_ \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 8 \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/mmul_index \
     --readFilesIn /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/RC_1_R1.fq.gz /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/RC_1_R2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/star_monkey/RC_ \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 8 \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/mmul_index \
     --readFilesIn /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/RT_1_R1.fq.gz /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/RT_1_R2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/star_monkey/RT_ \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate

# check raw RC fastq
STAR --runThreadN 8 \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/mmul_index \
     --readFilesIn /cluster/home/jhuang/projects/liver/data/zhangwei/pig/rnaseq/Data/RC_1/LC25B146nova_WK-2025021202959_S494_L002_R1.fastq.gz /cluster/home/jhuang/projects/liver/data/zhangwei/pig/rnaseq/Data/RC_1/LC25B146nova_WK-2025021202959_S494_L002_R2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/star_monkey/RC_raw_ \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 8 \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/sscrofa_index \
     --readFilesIn /cluster/home/jhuang/projects/liver/data/zhangwei/pig/rnaseq/Data/RC_1/LC25B146nova_WK-2025021202959_S494_L002_R1.fastq.gz /cluster/home/jhuang/projects/liver/data/zhangwei/pig/rnaseq/Data/RC_1/LC25B146nova_WK-2025021202959_S494_L002_R2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/star_pig/RC_raw_ \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate


# star index pig ref
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/sscrofa_index \
     --genomeFastaFiles /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.fna \
     --sjdbGTFfile /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/genomic.gtf \
     --sjdbOverhang 99

# align to pig ref
STAR --runThreadN 8 \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/sscrofa_index \
     --readFilesIn /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LC_1_R1.fq.gz /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LC_1_R2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/star_pig/LC_ \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 8 \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/sscrofa_index \
     --readFilesIn /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LT_1_R1.fq.gz /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LT_1_R2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/star_pig/LT_ \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 8 \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/sscrofa_index \
     --readFilesIn /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/RC_1_R1.fq.gz /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/RC_1_R2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/star_pig/RC_ \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 8 \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/sscrofa_index \
     --readFilesIn /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/RT_1_R1.fq.gz /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/RT_1_R2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/star_pig/RT_ \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate

