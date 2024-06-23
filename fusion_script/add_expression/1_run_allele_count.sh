#! /bin/bash

VCF=$1
bam_input=$2
output_file=$3



pos_list=$VCF

reference_input=$4

score=0

# consider only somatic mutation
#echo $run | xargs -I {} sh -c "{}"


#VCF=snvs_H059-MRAP5C_somatic_functional_snvs_conf_8_to_10.vcf
#bam_input=/icgc/dkfzlsdf/project/hipo/hipo_059/sequencing/whole_genome_sequencing/view-by-pid/H059-MRAP5C/tumor4/paired/merged-alignment/tumor4_H059-MRAP5C_merged.mdup.bam

grep -v "#CHROM" ${pos_list} | samtools mpileup -l - ${bam_input} -f ${reference_input} > ${output_file}_

#cut -f 4,5 ${pos_list} | samtools mpileup -RE -q 30 -ug -l - ${bam_input} -f ${reference_input} > ${output_file}

#echo "samtools mpileup -l - ${bam_input} -f ${reference_input} > ${output_file}" 

# the calculating is based on reads instead of fragment.
perl /icgc/dkfzlsdf/analysis/D120/scripts/add_expression/pileup2base_no_strand.pl ${output_file}_ $score ${output_file}

#perl ~/scripts/pile2base/pileup2base_no_strand.pl ${output_file}_ 1 ${output_file}_BS1
#perl ~/scripts/pile2base/pileup2base_no_strand.pl ${output_file}_ 5 ${output_file}_BS5
#perl ~/scripts/pile2base/pileup2base_no_strand.pl ${output_file}_ 10 ${output_file}_BS10
#perl ~/scripts/pile2base/pileup2base_no_strand.pl ${output_file}_ 15 ${output_file}_BS15
#perl ~/scripts/pile2base/pileup2base_no_strand.pl ${output_file}_ 20 ${output_file}_BS20

rm ${output_file}_

