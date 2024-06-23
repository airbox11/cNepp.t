set -e

pid=$1
result=$2
workDir=$3/3_add_expression
	# workDir=/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/${pid}/3_add_expression

cd $workDir
rm -rf $result



bam=`find . -name '*merged.mdup.bam'| xargs -I {} realpath {}`
ln -sf $bam.bai $workDir



## gene col6a3
# region="2:238232646-238323018"
# len=$((238323018-238232646)) 


## exon 6
region='2:238285415-238285987'
len=19624 # extronic


## reads :
reads=`samtools view -c -F 1540 $bam $region`
echo "reads=$reads" >> $result


## total reads: 
totalReads=`samtools view -c -F 1540 $bam`
echo "totalReads=$totalReads"  >> $result

## FPKM calculation
fpkm=`awk -vr=$reads -vR=$totalReads -vL=$len 'BEGIN{print r*10^9/(R*L)}'`
echo "fpkm=$fpkm" >> $result

## result
cat $result