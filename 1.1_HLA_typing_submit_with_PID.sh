# submit with PID, for example: S2914Nr1

sampleID=$1
workDir=$2

HLA_script=/icgc/dkfzlsdf/analysis/D120/scripts/HLA_typing/run_hla_typing.sh
fq1=/icgc/dkfzlsdf/analysis/D120/yanhong/data_collection_fastq.gz/*${sampleID}*1.fastq.gz
fq2=/icgc/dkfzlsdf/analysis/D120/yanhong/data_collection_fastq.gz/*${sampleID}*2.fastq.gz
zcat $fq2 | head 

# ln1=`zcat $fq2 | head | wc -l `
# echo $ln1


sh $HLA_script -i ${sampleID} -1 $fq1 -2 $fq2 -o ${workDir}/1.1_hla_typing

# cp -r /icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline_collection/HLA_predition/${sampleID}/phlat/phlat_${sampleID}_HLA.sum ${workDir}/1.1_hla_typing/