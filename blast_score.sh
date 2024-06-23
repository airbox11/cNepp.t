#!/usr/bin/env bash
# originate from /home/ubt21/lyuya/git/hex/20220322_blast/run_blast.sh
module load ncbi-blast/2.7.1

## public varaibles ==== ==== 


dir8=$1/8_chose_neoepitode
dirBlast=$dir8/blast_score
mkdir -p $dirBlast

sourceDir=/omics/groups/OE0422/internal/yanhong/git/hex/20220322_blast
virus_fasta=$sourceDir/sequences_human_virus_freseq_all.fasta
virus_database_dir=$sourceDir/virus_proteins_database

mhc_input=$2
query_fa=$dirBlast/${mhc_input}_query.fa
tmp=$dirBlast/${mhc_input}_tmp.csv
csv=$dirBlast/${mhc_input}_output.csv
txt=$dirBlast/${mhc_input}_output.txt


## function ==== ====

function prepare_query_fa {
    progf=/omics/groups/OE0422/internal/yanhong/git/awk_tool/t3.awk
    awk -f $progf -v col='Mutant_peptide' $dir8/$mhc_input |sort -u | awk '{print ">" $0; print $0}' > $query_fa
}


function make_database {
    makeblastdb \
    -title "virus-proteins" \
    -dbtype prot \
    -in  $virus_fasta \
    -out $virus_database_dir/vp
    # -parse_seqids \
    # -blastdb_version 5 \
}


function blastp_run {
    rm -f $csv

    blastp -task blastp-short -query $query_fa -db $virus_database_dir/vp -out $tmp -outfmt 6
    blastp -task blastp-short -query $query_fa -db $virus_database_dir/vp -out $txt -outfmt 1

    echo -e 'query_acc_ver\tsubject_acc_ver\tidentity\talignment_length\tmismatches\tgap_opens\tq_start\tq_end\ts_start\ts_end\tevalue\tbit_score' >> $csv
    cat $tmp >> $csv
    rm -f $tmp
}


## main execute ==== ====

prepare_query_fa
# make_database
blastp_run
