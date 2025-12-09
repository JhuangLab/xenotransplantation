# ================== Sscrofa11.1 =====================
# create rename map
cd /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.6_Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_assembly_report.txt
awk '!/^#/ {print $7, $11}' GCF_000003025.6_Sscrofa11.1_assembly_report.txt > rename_map.txt # manually add chrM at the last line

# rename fasta
awk 'BEGIN {
    while ((getline < "rename_map.txt") > 0) map[$1]=$2
}
{
    if ($0 ~ /^>/) {
        id = substr($0, 2); split(id, id_parts, " ");
        print ">" map[id_parts[1]] map[id_parts[0]]
    } else {
        print
    }
}' GCF_000003025.6_Sscrofa11.1_genomic.fna > GCF_000003025.6_Sscrofa11.1_genomic_renamed.fna
sed -i 's/\r$//' rename_map.txt # remove the return character of windows

# rename annotation
awk 'BEGIN {
    while ((getline < "rename_map.txt") > 0) map[$1]=$2
}
{
    if ($1 in map) $1 = map[$1];
    print
}' genomic.gtf > GCF_000003025.6_Sscrofa11.1_genomic_renamed.gtf

awk 'BEGIN {
    while ((getline < "rename_map.txt") > 0) map[$1]=$2
}
{
    if ($1 in map) $1 = map[$1];
    print
}' genomic.gff > GCF_000003025.6_Sscrofa11.1_genomic_renamed.gff

# star index
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/sscrofa_index \
     --genomeFastaFiles /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.fna \
     --sjdbGTFfile /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Sscrofa11.1/genomic.gtf \
     --sjdbOverhang 99

# ================== Mmul10 =====================
# create rename map
cd /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_assembly_report.txt
awk '!/^#/ {print $7, $11}' GCF_003339765.1_Mmul_10_assembly_report.txt > rename_map.txt # manually add chrM at the last line
sed -i 's/\r$//' rename_map.txt # remove the return character of windows

# rename fasta
awk 'BEGIN {
    while ((getline < "rename_map.txt") > 0) map[$1]=$2
}
{
    if ($0 ~ /^>/) {
        id = substr($0, 2); split(id, id_parts, " ");
        print ">" map[id_parts[1]] map[id_parts[0]]
    } else {
        print
    }
}' GCF_003339765.1_Mmul_10_genomic.fna > GCF_003339765.1_Mmul_10_genomic_renamed.fna

# rename annotation
awk 'BEGIN {
    while ((getline < "rename_map.txt") > 0) map[$1]=$2
}
{
    if ($1 in map) $1 = map[$1];
    print
}' genomic.gtf > GCF_003339765.1_Mmul_10_genomic_renamed.gtf

awk 'BEGIN {
    while ((getline < "rename_map.txt") > 0) map[$1]=$2
}
{
    if ($1 in map) $1 = map[$1];
    print
}' genomic.gff > GCF_003339765.1_Mmul_10_genomic_renamed.gff

# star index monkey ref
mkdir mmul_index
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/mmul_index \
     --genomeFastaFiles /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/GCF_003339765.1_Mmul_10_genomic_renamed.fna \
     --sjdbGTFfile /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/GCF_003339765.1_Mmul_10_genomic_renamed.gtf \
     --sjdbOverhang 99
mv Log.out /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/Mmul10/mmul_index/Log.out

# ===================== build a combined genome =======================
# combine genome fasta
cd /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/combined
sed -e 's/^>/>Mmul_/g' ../Mmul10/GCF_003339765.1_Mmul_10_genomic_renamed.fna > Mmul.fa
sed -e 's/^>/>Sscrofa_/g' ../Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic_renamed.fna > Sscrofa.fa
cat Sscrofa.fa Mmul.fa > Sscrofa_Mmul_combined_genome.fa

# combine gtf
grep -v "^#" ../Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic_renamed.gtf > pig_no_header.gtf
sed 's/^/Sscrofa_/g' pig_no_header.gtf > Sscrofa.gtf
grep -v "^#" ../Mmul10/GCF_003339765.1_Mmul_10_genomic_renamed.gtf > monkey_no_header.gtf
sed 's/^/Mmul_/g' monkey_no_header.gtf > Mmul.gtf
cat Sscrofa.gtf Mmul.gtf > Sscrofa_Mmul_combined.gtf # manully add header and end symbol
awk '
{
  split($1, a, "_");
  prefix = a[1];

  for (i = 1; i < NF; i++) {
    if ($i == "gene_id" && $(i+1) != "\"\";") {
      id = $(i+1);
      gsub(/^"|"[;]?$/, "", id); 
      $(i+1) = "\"" prefix "_" id "\";";
    }
    if ($i == "transcript_id" && $(i+1) != "\"\";") {
      id = $(i+1);
      gsub(/^"|"[;]?$/, "", id);
      $(i+1) = "\"" prefix "_" id "\";";
    }
  }

  out = $1;
  for (i = 2; i <= NF; i++) {
    out = out " " $i;
  }
  print out;
}' Sscrofa_Mmul_combined.gtf > Sscrofa_Mmul_combined_annotation.gtf


# build transcriptome fasta
# gffread -F -w Sscrofa_Mmul_combined_transcripts.fa -g Sscrofa_Mmul_combined_genome.fa Sscrofa_Mmul_combined.gtf
awk 'BEGIN{OFS="\t"} 
$3=="transcript" {
    # Extract transcript_id from attributes
    split($0, attr, "transcript_id \"");
    split(attr[2], tid, "\"");
    ttid = tid[1];
    if (++count[ttid] > 1) ttid = ttid "_dup" count[ttid]-1;

    # Extract gene_id
    split($0, gattr, "gene_id \"");
    split(gattr[2], gid, "\"");

    # transcript_biotype
    split($0, tattr, "transcript_biotype \"");
    split(tattr[2], bt, "\"");

    # Print BED format: chrom, start, end, name, score, strand
    print $1, $4-1, $5, ttid "|" gid[1] "|" $1 "|" $5-$4+1 "|" bt[1] "|", 0, $7
}' Sscrofa_Mmul_combined_annotation.gtf > transcripts.bed
bedtools getfasta -fi Sscrofa_Mmul_combined_genome.fa -bed transcripts.bed -nameOnly -s -fo Sscrofa_Mmul_combined_transcripts.fa


# salmon index
grep "^>" Sscrofa_Mmul_combined_genome.fa | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat Sscrofa_Mmul_combined_transcripts.fa Sscrofa_Mmul_combined_genome.fa > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -p 20 -i salmon_index --gencode

