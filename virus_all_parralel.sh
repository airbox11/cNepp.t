#!/bin/env bash

workDir=$1
inputID=$2

# workDir=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/MO79
# inputID=MHCI_epitopes_RNAseq_netMHCpan4_1.tab

## 
scriptDir=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/
runDir=$workDir/8_chose_neoepitode/blast_score
qurs_inpput=$runDir/${inputID}_renameCol_loh_query.fa
file_output=$runDir/${inputID}_renameCol_loh_query.fa_blscore

qurs_inpput=/omics/groups/OE0422/internal/yanhong/git/20220406_paula_hex_test/more_test_batch2/batch2.txt
file_output=/omics/groups/OE0422/internal/yanhong/git/20220406_paula_hex_test/more_test_batch2/batch2.txt_output

function bl62_align {

    rm -f $file_output
    readarray array < <(grep -v '>' $qurs_inpput)
    db_types=('virus' 'fungi' 'bacteria')

    script1=${scriptDir}/calculate_score.py
    echo '==== all virus alignment start [[[['
    for qur in ${array[@]};do
        for db_type in ${db_types[@]};do
            python $script1 $qur $file_output $db_type &
        done
    done
    wait
    echo '==== all virus alignment done ]]]]'

    cat $file_output
}

function get_microbialInfo {
    file_input=/omics/groups/OE0422/internal/yanhong/git/20220406_paula_hex_test/more_test_batch2/batch2.txt_output
    file_output=/omics/groups/OE0422/internal/yanhong/git/20220406_paula_hex_test/more_test_batch2/batch2.txt_output_virusName.tsv
    file_reference=/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py_all_virus_protein/all_virus_protein_sequences.fasta_seq_virusName

    rm -f $file_output
    echo $file_output
    echo -e 'qurey_seq\tref_seq\tBLscore\tmicrobial_type\tmicrobial_info' > $file_output

    while IFS= read -r line;do
        echo -ne "$line\t" >> $file_output
        ref=`echo $line | awk '{print $2}'`
        db_type=`echo $line | awk '{print $4}'`

        if [ $db_type == 'virus' ];then
            file_reference="/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py_all_virus_protein/all_virus_protein_sequences.fasta_seq_virusName"
        elif [ $db_type == 'bacteria' ];then
            file_reference="/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py_uniProt/all_protien_bacteria.fasta_nameSeq"
        elif [ $db_type == 'fungi' ];then
            file_reference="/omics/groups/OE0422/internal/yanhong/git/hex/virus_bac_database/splice_py_uniProt/all_protien_fungi.fasta_nameSeq"
        fi

        grep $ref $file_reference | 
        awk '{
            # a = gensub(/>.*\|(.*): .*/, "\\1; ", "g", $0)
            a = gensub(">(.*): .*", "\\1; ", "g", $0)
            printf "%s", a
        }' >> $file_output

        echo >> $file_output

    done < $file_input
    cat $file_output
}

function merge_files {
    Rscript ${scriptDir}/align_BL62_score_allVirus.r $workDir $inputID $file_output
}


bl62_align
# get_microbialInfo
merge_files