# combine quantification
salmon_index=/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/combined/salmon_index
gtf=/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/combined/Sscrofa_Mmul_combined.gtf

cd /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei
mkdir -p results/salmon_quant
salmon quant -p 20 -i $salmon_index -l A --validateMappings -o results/salmon_quant/LC -1 /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LC_1_R1.fq.gz -2 /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LC_1_R2.fq.gz -g $gtf

salmon quant -p 20 -i $salmon_index -l A --validateMappings -o results/salmon_quant/LT -1 /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LT_1_R1.fq.gz -2 /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/fastq/LT_1_R2.fq.gz -g $gtf

# tx2gene
cd /cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/ref/combined
awk 'BEGIN{OFS="\t"} 
$3 == "transcript" || $3 == "gene" {
	split($0, attr, "transcript_id \"");
    split(attr[2], tid, "\"");

    split($0, gattr, "gene_id \"");
    split(gattr[2], gid, "\"");

    split($0, tattr, "transcript_biotype \"");
    split(tattr[2], bt, "\"");

    split($1, sp, "_");

    if ($3 == "transcript" && tid[1] != "" && gid[1] != "") {
        print tid[1], gid[1], bt[1], sp[1];
    }
    else if ($3 == "gene" && gid[1] != "") {
        print gid[1], gid[1], bt[1], sp[1]; 
    }
}' Sscrofa_Mmul_combined_annotation.gtf > tx2gene_combined.tsv
