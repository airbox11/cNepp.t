
# Use previous script to take fusion peptide as input to do fusion epitope prediction. 
## must be improved using scripts in /icgc/dkfzlsdf/analysis/D120/scripts/fusion/

exit

folder=/icgc/dkfzlsdf/analysis/D120/huangz/Fusion_liposarcoma_project/fusionPepDir/
outputDir=/icgc/dkfzlsdf/analysis/D120/huangz/Fusion_liposarcoma_project/outputDir
hlaFilePath=/icgc/dkfzlsdf/analysis/D120/huangz/Fusion_liposarcoma_project/hla_file_path2

for file in `ls $folder`
do
	sampleID=`echo $file | awk '{split($1,a,"_");print a[1]}'`
	hlaFile=`grep $sampleID $hlaFilePath | sort | uniq | head -n 1 `
	pathPepFa=$folder/$file
	mutType=`echo $file | sed 's/.fa//g'`
	echo "sh /home/huangz/scripts/immunoinfo/fusion_neoantigen/main_run_fusion_neoantigen_prediction.sh $pathPepFa $mutType $outputDir $hlaFile"

done
