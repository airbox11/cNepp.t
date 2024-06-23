

## 27mer
#snv
/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/update_netMHCpan/netMHCpan_4.1/main_run_neoepitope_snvs.sh
/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/update_netMHCpan/netMHCpan_4.1/run_netMHC.sh

#indel:
/omics/groups/OE0422/internal/scripts/neoepitope_indels/indels_netmhcpan.sh
i4a_indel_predict

#fusion:
/omics/groups/OE0422/internal/yanhong/git/neoepitope_prediction-master1/immunopipe/run_netMHC_on_consensus.sh

######


all_in_one_pipeline.sh S2975Nr1_COAD S2975Nr1 JD52 No TCGA-COAD &


runID=S2975Nr1_COAD
sampleID=S2975Nr1
patientID=JD52
RNAseq=No
tcga=TCGA-COAD


all_in_one_pipeline.sh S2975Nr3 S2975Nr3 ZK55 No TCGA-PAAD &
runID=S2975Nr3
sampleID=S2975Nr3
patientID=ZK55
RNAseq=No
tcga=TCGA-PAAD
scriptsDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline
workDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline_collection/${runID}


all_in_one_pipeline.sh S3005Nr1 S3005Nr1 GO53 No TCGA-LIHC &
runID=S3005Nr1
sampleID=S3005Nr1
patientID=GO53
RNAseq=No
tcga=TCGA-LIHC
scriptsDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline
workDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline_collection/${runID}



all_in_one_pipeline.sh S2975Nr1_READ S2975Nr1 JD52 No TCGA-READ &


all_in_one_pipeline.sh S2975Nr3_DkmsHLA S2975Nr3 ZK55 No TCGA-PAAD &


runID=S2975Nr3_DkmsHLA
sampleID=S2975Nr3
patientID=ZK55
RNAseq=No
tcga=TCGA-PAAD
scriptsDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline
workDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline_collection/${runID}



all_in_one_pipeline.sh S3005Nr1_4 S3005Nr1 GO53 No TCGA-LIHC &
runID=S3005Nr1_4.1
sampleID=S3005Nr1
patientID=GO53
RNAseq=No
tcga=TCGA-LIHC
scriptsDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline
workDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline_collection/${runID}
mhcIversion=4.1


all_in_one_pipeline.sh S3005Nr1_4 S3005Nr1 GO53 No TCGA-LIHC &
sampleID=i034
workDir=/icgc/dkfzlsdf/analysis/D120/yanhong/all_in_one_pipeline_collection/i034


############

all_in_one_pipeline.sh S2914Nr2 S2914Nr2 AbiH64 No TCGA-COAD 4.1 &
all_in_one_pipeline.sh S2914Nr4 S2914Nr4 MD63 Yes RNAseq 4.1 &
all_in_one_pipeline.sh S2975Nr1-COAD S2975Nr1 JD52 No TCGA-COAD 4.1 &
all_in_one_pipeline.sh S2975Nr1-READ S2975Nr1 JD52 No TCGA-READ 4.1 &
all_in_one_pipeline.sh S2975Nr3 S2975Nr3 ZK55 No TCGA-PAAD 4.1 &
all_in_one_pipeline.sh S2975Nr3-Dkms S2975Nr3 ZK55 No TCGA-PAAD 4.1 &
all_in_one_pipeline.sh S3005Nr1  S2914Nr1  GO53 No TCGA-LIHC 4.1 &

############
all_in_one_pipeline.sh S2914Nr2 S2914Nr2 AbiH64 No TCGA-COAD 4.1
runID=S2914Nr2
sampleID=S2914Nr2
patientID=AbiH64
RNAseq=No
tcga=TCGA-COAD
mhcIversion=4.1


all_in_one_pipeline.sh S2914Nr4 S2914Nr4 MD63 Yes RNAseq 4.1 &
runID=S2914Nr4
sampleID=S2914Nr4
patientID=MD63
RNAseq=Yes
tcga=RNAseq
mhcIversion=4.1


all_in_one_pipeline.sh S2914Nr4 S2914Nr4 MD63 Yes RNAseq 4.0 &
runID=S2914Nr4
sampleID=S2914Nr4
patientID=MD63
RNAseq=Yes
tcga=RNAseq
mhcIversion=4.0

############

all_in_one_pipeline.sh S2914Nr2 S2914Nr2 AbiH64 No TCGA-COAD 4.1 &
all_in_one_pipeline.sh S2914Nr4 S2914Nr4 MD63 Yes RNAseq 4.1 &

all_in_one_pipeline.sh S2975Nr1-COAD S2975Nr1 JD52 No TCGA-COAD 4.1 &

all_in_one_pipeline.sh S2975Nr1-READ S2975Nr1 JD52 No TCGA-READ 4.1 &
all_in_one_pipeline.sh S2975Nr3 S2975Nr3 ZK55 No TCGA-PAAD 4.1 &
all_in_one_pipeline.sh S2975Nr3-Dkms S2975Nr3 ZK55 No TCGA-PAAD 4.1 &

all_in_one_pipeline.sh S3005Nr1 S3005Nr1 GO53 No TCGA-LIHC 4.1 &

###  TRC

all_in_one_pipeline.sh trcpoc004HLA trcpoc004HLA trcpoc004HLA Yes RNAseq 4.1 &

all_in_one_pipeline.sh trcpoc004HLA trcpoc004HLA trcpoc004HLA Yes RNAseq 4.1 s4&


### buffy coat:
all_in_one_pipeline.sh bc01p bc01p bc01p No RNAseq 4.1 s1b&
all_in_one_pipeline.sh bc01p2 bc01p2 bc01p2 No RNAseq 4.1 s2s3s4 &
all_in_one_pipeline.sh bc01pr1_17 bc01pr1_17 bc01pr1_17 No RNAseq 4.1 s1.2&
all_in_one_pipeline.sh bc01pr1_18 bc01pr1_18 bc01pr1_18 No RNAseq 4.1 s1.2&


###
all_in_one_pipeline.sh bc01p2_low bc01p2_low bc01p2_low Yes RNAseq 4.1 s5.1s5.2 &
all_in_one_pipeline.sh bc01p2_upper bc01p2_upper bc01p2_upper Yes RNAseq 4.1 s5.1s5.2 &

all_in_one_pipeline.sh bc01p2_low bc01p2_low bc01p2_low Yes RNAseq 4.1 s3 &

## HLA for IPNCT_RW67, IPNCT_PW60 (s3217)
## indel 

all_in_one_pipeline.sh S3217Nr1 S3217Nr1 RW67 Yes RNAseq  4.1 s4d &
all_in_one_pipeline.sh S3217Nr3 S3217Nr3 PW60 Yes RNAseq  4.1 s4d &

all_in_one_pipeline.sh indel_test2 indel_test2 indel_test2 Yes RNAseq  4.1 s4c &


## 3228, IL85
all_in_one_pipeline.sh NR78 S3328Nr1 NR78 Yes RNAseq  4.1 s4as4bs4c &

all_in_one_pipeline.sh EB72 S3328Nr3 EB72 Yes RNAseq  4.1 s2s3s8as8b &
all_in_one_pipeline.sh EB72 S3328Nr3 EB72 Yes RNAseq  4.1 s4as4bs4c &

all_in_one_pipeline.sh IL85 IL85 IL85 No TCGA-CESC 4.1 s4d &

## indel xlsx test
all_in_one_pipeline.sh S2975Nr1-COAD S2975Nr1 JD52_yh Yes RNAseq  4.1 s8d &

## indel debug sequences
all_in_one_pipeline.sh S3217Nr3 S3217Nr3 PW60 Yes RNAseq  4.1 s4a &

## re-run indel for IL85 NR78 EB72
all_in_one_pipeline.sh IL85 IL85 IL85 No TCGA-CESC 4.1 s4as4bs4c &
all_in_one_pipeline.sh NR78 S3328Nr1 NR78 Yes RNAseq  4.1 s4as4bs4c &
all_in_one_pipeline.sh EB72 S3328Nr3 EB72 Yes RNAseq  4.1 s4as4bs4c &

## move results to omics
all_in_one_pipeline.sh S3217Nr1 S3217Nr1 RW67 Yes RNAseq  4.1 s4as4bs4c &
all_in_one_pipeline.sh S3217Nr3 S3217Nr3 PW60 Yes RNAseq  4.1 s4as4bs4c &

all_in_one_pipeline.sh S3005Nr1 S3005Nr1 GO53 No TCGA-LIHC 4.1 s8b &

all_in_one_pipeline.sh S2914Nr2 S2914Nr2 AbiH64 No TCGA-COAD 4.1 s8b &
all_in_one_pipeline.sh S2914Nr4 S2914Nr4 MD63 Yes RNAseq 4.1 s8b &
all_in_one_pipeline.sh S2975Nr1-COAD S2975Nr1 JD52 No TCGA-COAD 4.1 s8b &
all_in_one_pipeline.sh S2975Nr3-dkms S2975Nr3 ZK55 No TCGA-PAAD 4.1 s8b &

all_in_one_pipeline.sh trcpoc004 trcpoc004 trcpoc004 Yes RNAseq 4.1 s8b &

## indel debug2
# all_in_one_pipeline.sh S3217Nr1 S3217Nr1 RW67 Yes RNAseq 4.1 s4as4bs4c > S3217Nr1.tmp 2>&1  &
all_in_one_pipeline.sh S3217Nr3 S3217Nr3 PW60 Yes RNAseq 4.1 s4as4bs4c > PW60.tmp 2>&1  &
all_in_one_pipeline.sh NR78 S3328Nr1 NR78 Yes RNAseq 4.1 s4as4bs4c > NR78.tmp 2>&1  &
all_in_one_pipeline.sh EB72 S3328Nr3 EB72 Yes RNAseq 4.1 s4as4bs4c > EB72.tmp 2>&1  &
all_in_one_pipeline.sh IL85 IL85 IL85 No TCGA-CESC 4.1 s4as4bs4c > IL85.tmp 2>&1  &


## K26K-MK6UTZ
all_in_one_pipeline K26K-MK6UTZ 52656 HS80 Yes RNAseq 4.1a s1b
all_in_one_pipeline HS80 52656 HS80 Yes RNAseq 4.1 s2s3s8as8b
all_in_one_pipeline HS80 52656 HS80 Yes RNAseq 4.1 s4as4bs4c

## omics move ====
all_in_one_pipeline NR78 S3328Nr1 NR78 Yes RNAseq 4.1 s4cs8b 


## catch -dataset:
all_in_one_pipeline K26K-tumor 39032 K26K yes RNAseq 4.1 s1a &

all_in_one_pipeline K26K-HGXQVM sampleID K26K-HGXQVM Yes RNAseq 4.1 s2 &

## 141 debug

all_in_one_pipeline K26K-HGXQVM-metastasis22-01 sampleID K26K-HGXQVM-metastasis22-01 Yes RNAseq 4.1 s4a &
all_in_one_pipeline K26K-HGXQVM-tumor11 sampleID K26K-HGXQVM-tumor11 Yes RNAseq 4.1 s4a &

all_in_one_pipeline K26K-J1JGYX-metastasis11-01-p-buffy-coat1-01-p2 sampleID K26K-J1JGYX-metastasis11-01-p-buffy-coat1-01-p2 Yes RNAseq 4.1 s3 &

all_in_one_pipeline K26K-J1JGYX-metastasis2-01-p-buffy-coat1-01-p snv  Yes RNAseq 4.1 s4b 

## tcga 
all_in_one_pipeline K26K-HGXQVM-tumor11 sampleID K26K-HGXQVM-tumor11 No TCGA-BRCA 4.1 s3 &

all_in_one_pipeline K26K-HGXQVM-metastasis22-01 sampleID K26K-HGXQVM-metastasis22-01 No TCGA-BRCA 4.1 s3 &


##  HL50
all_in_one_pipeline HL50-3 sampleID HL50-3 No TCGA-BRCA s8b debug
all_in_one_pipeline HL50-4 sampleID HL50-3 No TCGA-BRCA s8b debug

all_in_one_pipeline K26K-MK6UTZ sampleID K26K-MK6UTZ No TCGA-BRCA s3s8as8b &
all_in_one_pipeline K26K-MK6UTZ sampleID K26K-MK6UTZ No TCGA-BRCA i4ai4bi4c &


## ts73
all_in_one_pipeline TS73 sampleID TS73 Yes RNAseq s2s3s8as8b &
all_in_one_pipeline LF89 sampleID LF89 Yes RNAseq s2s3s8as8b &


all_in_one_pipeline TS73 sampleID TS73 Yes RNAseq s8b &
all_in_one_pipeline LF89 sampleID LF89 Yes RNAseq s8b &


## IT01 vg53
all_in_one_pipeline IT01 sampleID IT01 No TCGA-THYM s1bs2s3s8as8bi4ai4bi4c
all_in_one_pipeline VG53 sampleID VG53 No TCGA-COAD s1bs2s3s8as8bi4ai4bi4c



## promise S014-1HBNPG
all_in_one_pipeline \
-r S014-1HBNPG \
-a sampleID \
-p S014-1HBNPG \
-n No \
-t RNAseq \
-s s1a \
-d debug \
-h 68_R 

all_in_one_pipeline -r S014-1HBNPG -a sampleID -p S014-1HBNPG -n No -t RNAseq -s s1a -h 68_R
all_in_one_pipeline -r S014-1HBNPG -a sampleID -p S014-1HBNPG -n No -t RNAseq -s s1a -h 69_R
all_in_one_pipeline -r S014-1HBNPG -a sampleID -p S014-1HBNPG -n No -t RNAseq -s s1a -h 70_R 
all_in_one_pipeline -r S014-1HBNPG -a sampleID -p S014-1HBNPG -n No -t RNAseq -s s1a -h 71_R 

## IT01 vg53

all_in_one_pipeline -r IT01 -a sampleID -p IT01 -n No -t TCGA-THYM -s i4c
all_in_one_pipeline -r VG53 -a sampleID -p VG53 -n No -t TCGA-COAD -s i4c


## 
all_in_one_pipeline -r S014-1HBNPG -a sampleID -p S014-1HBNPG -n No -t RNAseq -s s1a -h 68_R -d debug

all_in_one_pipeline -r S014-1HBNPG -a sampleID -p S014-1HBNPG -n No -t RNAseq -s s1a -h 68_R 


## K26K-BZC91B 
all_in_one_pipeline -r K26K-BZC91B -s s1a -h manual 
all_in_one_pipeline -r K26K-P7V4JK -s s0 


## HL50
all_in_one_pipeline -r HL50-3 -t TCGA-BRCA -s i4a -d debug
all_in_one_pipeline -r HL50-4 -t TCGA-BRCA -s s3s8as8b 


## stabpan:

all_in_one_pipeline -r EB72 -s s2 -d debug

all_in_one_pipeline -r EB72 -s s2 -e netMHCstabpan -d debug

all_in_one_pipeline -r EB72 -e netMHCstabpan -n Yes -t RNAseq  -s s3 -d debug
all_in_one_pipeline -r EB72 -e netMHCstabpan -s s3 -d debug


all_in_one_pipeline -r EB72 -e netMHCstabpan -t RNAseq  -s s8as8b
all_in_one_pipeline -r EB72 -e netMHCpan4_1  -t RNAseq  -s s8as8b


all_in_one_pipeline -r EB72 -e netMHCstabpan -s s8as8b
all_in_one_pipeline -r EB72 -e netMHCpan4_1  -s s8as8b


all_in_one_pipeline -r EB72 -s s1c -d debug



## master3
all_in_one_pipeline -r VT71 -s s1a 
all_in_one_pipeline -r OH61 -s s1a
all_in_one_pipeline -r HX66NQ -s s1a


arbpipe -r VT71 -s s0b
arbpipe -r OH61 -s s0b
arbpipe -r HX66NQ -s s0b

all_in_one_pipeline -r EB72 -s s1b -d debug
all_in_one_pipeline -r S014-1HBNPG -s s1b -h promise -d debug 


# all_in_one_pipeline -r EB72 -e netMHCpan4_1 -s s2 -d debug

# all_in_one_pipeline -r VT71 -t RNAseq -e netMHCstabpan -s s3s8as8b -d debug
# all_in_one_pipeline -r OH61 -t RNAseq -e netMHCstabpan -s s8as8b -d debug
# all_in_one_pipeline -r HX66NQ -t RNAseq -e netMHCstabpan -s s2s3s8as8b


# all_in_one_pipeline -r K26K-1GUTR3   -t RNAseq -e netMHCpan4_1 -s s0 -d debug

# all_in_one_pipeline -r K26K-1GUTR3   -t RNAseq -e netMHCpan4_1 -s s3s8as8b
# all_in_one_pipeline -r K26K-BZC91B-2 -t RNAseq -e netMHCpan4_1 -s s2s3s8as8b
# all_in_one_pipeline -r K26K-1GUTR3  -t RNAseq -s i4c -d debug


# all_in_one_pipeline -r K26K-BZC91B-2 -s indel

# arbpipe -r K26K-1GUTR3 -s s5 -d debug
# arbpipe -r K26K-BZC91B-2 -s s5 -d debug


# all_in_one_pipeline -r K26K-BZC91B-2 -s s1b -d debug

# arbpipe -r VT71 -s s1bs2s3s4s5 -d debug
# arbpipe -r OH61 -s s4s5 -d debug
# arbpipe -r HX66NQ -s s4s5 -d debug


# all_in_one_pipeline -r EB72 -s s2 -d debug -e netMHCpan4_1


# all_in_one_pipeline -r EB72 -s s2 -d debug

# =====================

# all_in_one_pipeline -r S014-1HBNPG -h promise -e netMHCpan4_1 -s i4c -d debug
# all_in_one_pipeline -r S014-1HPKHP -h promise -e netMHCpan4_1 -s i4c -d debug
# all_in_one_pipeline -r S014-2CDKKU -h promise -e netMHCpan4_1 -s i4bi4c -d debug
# all_in_one_pipeline -r S014-2UF991 -h promise -e netMHCpan4_1 -s i4a -d debug

all_in_one_pipeline -r S014-2W4HUU -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-39NDBF -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-3CZ28M -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-3DMGUU -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-3SLKPB -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-3X4VS6 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-4DBUYX -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-4LNKDY -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-531DEJ -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-5HMWPS -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-5N9GYN -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-6UWZCH -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-8B6PB4 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-9L7E3P -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-9UAKU5 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-9YE8T4 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-A9GRM7 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-BGYMAE -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-BYZX35 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-CDK4GP -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-CZPM27 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-D249ZS -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-D2EW7K -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-DDXRZC -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-DTJ3Q9 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-DX441S -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-FFAWYJ -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-GYUXWL -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-H1A6K7 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-JPFZPE -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-K6L83S -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-KPGX4P -h promise -e netMHCpan4_1 -s indel
# all_in_one_pipeline -r S014-P1ECRS -h promise -e netMHCpan4_1 -s indel  no indel vcf
all_in_one_pipeline -r S014-PQRN8X -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-QNHL6F -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-RG1UEN -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-SHQDHT -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-SJ1GLQ -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-SQVZUQ -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-TC14F4 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-U8ETE8 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-UFL1D9 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-VAWMM3 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-VLRS4K -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-VR9K44 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-WW1C6Y -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-X3Q2S2 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-X47792 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-XD89A7 -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-YRMC2D -h promise -e netMHCpan4_1 -s indel
all_in_one_pipeline -r S014-YWR6N4 -h promise -e netMHCpan4_1 -s indel








all_in_one_pipeline -r S014-VAWMM3 -h promise -e netMHCpan4_1 -s f2 -d debug
all_in_one_pipeline -r S014-QNHL6F -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-JPFZPE -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-SJ1GLQ -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-3X4VS6 -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-3DMGUU -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-CZPM27 -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-BGYMAE -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-8B6PB4 -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-6UWZCH -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-5N9GYN -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-531DEJ -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-4DBUYX -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-VLRS4K -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-DX441S -h promise -e netMHCpan4_1 -s fusion 
all_in_one_pipeline -r S014-2CDKKU -h promise -e netMHCpan4_1 -s fusion 




# all_in_one_pipeline -r S014-VAWMM3 -h promise -t RNAseq -e netMHCpan4_1 -s s3 -d debug
all_in_one_pipeline -r S014-QNHL6F -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-JPFZPE -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-SJ1GLQ -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-3X4VS6 -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-3DMGUU -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-CZPM27 -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-BGYMAE -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-8B6PB4 -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-6UWZCH -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-5N9GYN -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-531DEJ -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-4DBUYX -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
all_in_one_pipeline -r S014-VLRS4K -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b
# all_in_one_pipeline -r S014-DX441S -h promise -t RNAseq -e netMHCpan4_1 -s snv
all_in_one_pipeline -r S014-2CDKKU -h promise -t RNAseq -e netMHCpan4_1 -s s3s8as8b


## batch2 17 s2
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-3X4VS6_T1T2_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-4DBUYX_T1T2_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-531DEJ_T1T2_tumor11
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-8B6PB4_T1T2_tumor11
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-9UAKU5_T2ERW_tumor31
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-BGYMAE_T12_tumor12
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-GYUXWL_T1T2_tumor11
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-GYUXWL_T3_tumor41
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-QNHL6F_T1T2_tumor11
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-QNHL6F_T3_tumor41
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-SJ1GLQ_T1T2_tumor1
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-SQVZUQ_T1T2_tumor11
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-U8ETE8_T1T2_tumor1
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-WW1C6Y_T1T2_tumor11
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-WW1C6Y_T12_tumor12
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-X47792_T3_tumor41
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s2 -r S014-YWR6N4_T1T2_tumor11



## batch2 17 s3 rna
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s3 -r S014-3X4VS6_tumor11
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s3 -r S014-4DBUYX_tumor11
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s3 -r S014-531DEJ_tumor11
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s3 -r S014-8B6PB4_tumor11
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s3 -r S014-BGYMAE_tumor12
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s3 -r S014-QNHL6F_tumor11
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s3 -r S014-QNHL6F_tumor41

## batch2 17 s3 tcga
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s3 -r S014-9UAKU5_tumor31 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s3 -r S014-GYUXWL_tumor11
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s3 -r S014-GYUXWL_tumor41
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s3 -r S014-SJ1GLQ_tumor1
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s3 -r S014-SQVZUQ_tumor11
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s3 -r S014-U8ETE8_tumor1
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s3 -r S014-WW1C6Y_tumor11
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s3 -r S014-WW1C6Y_tumor12
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s3 -r S014-X47792_tumor41
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s3 -r S014-YWR6N4_tumor11


## batch2 indel
# all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s i4c -r S014-3X4VS6_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s i4c -r S014-4DBUYX_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s i4c -r S014-531DEJ_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s i4c -r S014-8B6PB4_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s i4c -r S014-BGYMAE_tumor12 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s i4c -r S014-QNHL6F_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s i4c -r S014-QNHL6F_tumor41 -d debug

all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s i4c -r S014-9UAKU5_tumor31 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s i4c -r S014-GYUXWL_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s i4c -r S014-GYUXWL_tumor41 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s i4c -r S014-SJ1GLQ_tumor1 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s i4c -r S014-SQVZUQ_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s i4c -r S014-U8ETE8_tumor1 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s i4c -r S014-WW1C6Y_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s i4c -r S014-WW1C6Y_tumor12 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s i4c -r S014-X47792_tumor41 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s i4c -r S014-YWR6N4_tumor11 -d debug



## batch2 fusion
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-3X4VS6_tumor11 -s f2 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-531DEJ_tumor11 -s f2 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-BGYMAE_tumor12 -s f2 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-QNHL6F_tumor11 -s f2 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-4DBUYX_tumor11 -s f2 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-8B6PB4_tumor11 -s f2 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-QNHL6F_tumor41 -s f2 -d debug


S014-3X4VS6_tumor11
S014-531DEJ_tumor11
S014-BGYMAE_tumor12
S014-QNHL6F_tumor11




## H021-S4CLSR
all_in_one_pipeline -e netMHCpan4_1 -r H021-S4CLSR -s s1b -d debug 
all_in_one_pipeline -e netMHCpan4_1 -r H021-S4CLSR-tumor -s s1b -d debug 


all_in_one_pipeline -e netMHCpan4_1 -t RNAseq -r VT71 -s s3 -d debug 



##  17 of 72 sample : fusion based prediction

all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-3X4VS6_tumor11 -s f5 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-531DEJ_tumor11 -s f4f5 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-BGYMAE_tumor12 -s f4f5 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-QNHL6F_tumor11 -s f4f5 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-4DBUYX_tumor11 -s f4f5 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-8B6PB4_tumor11 -s f4f5 -d debug
all_in_one_pipeline -h promise -e netMHCpan4_1 -r S014-QNHL6F_tumor41 -s f4f5 -d debug

## 38 of 55 sample: snv, luad-lusc re-run 

all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-1HBNPG_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-1HPKHP_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-2UF991_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-2W4HUU_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-39NDBF_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-3CZ28M_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-3SLKPB_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-4LNKDY_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-5HMWPS_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-9L7E3P_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-9UAKU5_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-9YE8T4_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-A9GRM7_tumor12 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-BYZX35_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-CDK4GP_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-D249ZS_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-D2EW7K_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-DTJ3Q9_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-FFAWYJ_tumor1 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-GYUXWL_tumor61 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-H1A6K7_tumor1 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-K6L83S_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-KPGX4P_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-P1ECRS_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-PQRN8X_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-RG1UEN_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-SHQDHT_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-SQVZUQ_tumor12 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-TC14F4_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-U8ETE8_tumor31 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-UFL1D9_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-VR9K44_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-WW1C6Y_tumor41 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-X3Q2S2_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-X47792_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-XD89A7_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-YRMC2D_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s8b -r S014-YWR6N4_tumor41 -d debug

## OH61 X issue
all_in_one_pipeline -r OH61 -t RNAseq -e netMHCstabpan -s i4a -d debug



## sr66
all_in_one_pipeline -r sr66 -t RNAseq -e netMHCstabpan -d debug

## ak64
all_in_one_pipeline -r ak64 -t RNAseq -e netMHCpan4_1 -s snv-indel-fusion


## PE67
all_in_one_pipeline -r PE67 -t RNAseq -e netMHCpan4_1 -s f1b -d debug

## 
all_in_one_pipeline -e netMHCpan4_1 -r H021-S4CLSR-tumor -t RNAseq -s s3 -d debug 
all_in_one_pipeline -e netMHCpan4_1 -r MD45 -t RNAseq -s s3 -d debug 

##
all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s s1b -t RNAseq -d debug 
all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s snv -t RNAseq

	all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s s8a -t RNAseq -d debug 
	all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s indel -t RNAseq
	all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s i4c -t RNAseq -d debug
	all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s fusion -t RNAseq
	all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s f4 -t RNAseq -d debug



## 6A8HTK
all_in_one_pipeline -r H021-6A8HTK -t RNAseq -e netMHCpan4_1 -d debug
all_in_one_pipeline -r H021-6A8HTK -t RNAseq -e netMHCpan4_1 -s s1b 
all_in_one_pipeline -r H021-6A8HTK -t RNAseq -e netMHCpan4_1 -s s1b -d debug
all_in_one_pipeline -r H021-6A8HTK -t RNAseq -e netMHCpan4_1 -s snv
all_in_one_pipeline -r H021-6A8HTK -t RNAseq -e netMHCpan4_1 -s indel
	all_in_one_pipeline -r H021-6A8HTK -t RNAseq -e netMHCpan4_1 -s indel

all_in_one_pipeline -r H021-6A8HTK -t RNAseq -e netMHCpan4_1 -s fusion
	all_in_one_pipeline -r H021-6A8HTK -t RNAseq -e netMHCpan4_1 -s f4f5 -d debug

## EO65, LS90
all_in_one_pipeline -r EO65 -t RNAseq -e netMHCpan4_1 -d debug
	all_in_one_pipeline -r EO65 -t RNAseq -e netMHCpan4_1 -s s1b -d debug
	all_in_one_pipeline -r EO65 -t RNAseq -e netMHCpan4_1 -s snv
		all_in_one_pipeline -r EO65 -t RNAseq -e netMHCpan4_1 -s s8b -d debug
	all_in_one_pipeline -r EO65 -t RNAseq -e netMHCpan4_1 -s indel
	all_in_one_pipeline -r EO65 -t RNAseq -e netMHCpan4_1 -s fusion
		all_in_one_pipeline -r EO65 -t RNAseq -e netMHCpan4_1 -s f4f5 -d debug

all_in_one_pipeline -r LS90 -t RNAseq -e netMHCpan4_1 -d debug
	all_in_one_pipeline -r LS90 -t RNAseq -e netMHCpan4_1 -s s1b
	all_in_one_pipeline -r LS90 -t RNAseq -e netMHCpan4_1 -s s1b -d debug
	all_in_one_pipeline -r LS90 -t RNAseq -e netMHCpan4_1 -s snv
	all_in_one_pipeline -r LS90 -t RNAseq -e netMHCpan4_1 -s indel
	all_in_one_pipeline -r LS90 -t RNAseq -e netMHCpan4_1 -s fusion
	all_in_one_pipeline -r LS90 -t RNAseq -e netMHCpan4_1 -s s8bi4cf4f5 -d debug

## K26K-A2EDZS_meta2
all_in_one_pipeline -r K26K-A2EDZS_meta2 -t RNAseq -e netMHCpan4_1 -d debug
	all_in_one_pipeline -r K26K-A2EDZS_meta2 -t RNAseq -e netMHCpan4_1 -s s1b
	all_in_one_pipeline -r K26K-A2EDZS_meta2 -t RNAseq -e netMHCpan4_1 -s s1b -d debug
	all_in_one_pipeline -r K26K-A2EDZS_meta2 -t RNAseq -e netMHCpan4_1 -s snv
	all_in_one_pipeline -r K26K-A2EDZS_meta2 -t RNAseq -e netMHCpan4_1 -s indel
	all_in_one_pipeline -r K26K-A2EDZS_meta2 -t RNAseq -e netMHCpan4_1 -s fusion
	all_in_one_pipeline -r K26K-A2EDZS_meta2 -t RNAseq -e netMHCpan4_1 -s s8bi4cf4f5 -d debug

## Y4X8ML
all_in_one_pipeline -r H021-Y4X8ML -t RNAseq -e netMHCpan4_1 -d debug
	all_in_one_pipeline -r H021-Y4X8ML -t RNAseq -e netMHCpan4_1 -s s1b
	all_in_one_pipeline -r H021-Y4X8ML -t RNAseq -e netMHCpan4_1 -s s1b -d debug
	all_in_one_pipeline -r H021-Y4X8ML -t RNAseq -e netMHCpan4_1 -s snv
	all_in_one_pipeline -r H021-Y4X8ML -t RNAseq -e netMHCpan4_1 -s indel
	all_in_one_pipeline -r H021-Y4X8ML -t RNAseq -e netMHCpan4_1 -s fusion
	all_in_one_pipeline -r H021-Y4X8ML -t RNAseq -e netMHCpan4_1 -s s8bi4cf4f5 -d debug

## ak46
all_in_one_pipeline -r AK46 -t RNAseq -e netMHCpan4_1 -d debug
	all_in_one_pipeline -r AK46 -t RNAseq -e netMHCpan4_1 -s s1b
	all_in_one_pipeline -r AK46 -t RNAseq -e netMHCpan4_1 -s s1b -d debug
	all_in_one_pipeline -r AK46 -t RNAseq -e netMHCpan4_1 -s snv
	all_in_one_pipeline -r AK46 -t RNAseq -e netMHCpan4_1 -s indel
	all_in_one_pipeline -r AK46 -t RNAseq -e netMHCpan4_1 -s fusion
	all_in_one_pipeline -r AK46 -t RNAseq -e netMHCpan4_1 -s s8bi4cf4f5 -d debug
## aa52
all_in_one_pipeline -r AA52 -t RNAseq -e netMHCpan4_1 -d debug
	all_in_one_pipeline -r AA52 -t RNAseq -e netMHCpan4_1 -s s1b
	all_in_one_pipeline -r AA52 -t RNAseq -e netMHCpan4_1 -s s1b -d debug
	all_in_one_pipeline -r AA52 -t RNAseq -e netMHCpan4_1 -s snv
	all_in_one_pipeline -r AA52 -t RNAseq -e netMHCpan4_1 -s indel
	all_in_one_pipeline -r AA52 -t RNAseq -e netMHCpan4_1 -s fusion
	all_in_one_pipeline -r AA52 -t RNAseq -e netMHCpan4_1 -s s8bi4cf4f5 -d debug



## K26K-AT59U7
all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s s3 -t RNAseq

all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s indel -t RNAseq
all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s fusion -t RNAseq
all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s s8bi4cf4f5 -t RNAseq
all_in_one_pipeline -e netMHCpan4_1 -r K26K-AT59U7 -s i4c -t RNAseq -d debug

## promise 21 =================

## run hla type
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-9C2UWH_tumor12 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-NHGSH5_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-E9NNVR_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-FYVYBF_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-AGLMAZ_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-XTGWF2_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-JALDBT_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-BAJUW7_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-6617PL_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-UPSW7P_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-VVXZ39_tumor11 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s1b -r S014-ZNWXFT_tumor11 -d debug


## run prediction
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-3SLKPB_tumor41 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-4LNKDY_tumor41 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-6617PL_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-8B6PB4_tumor61 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-9C2UWH_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-9C2UWH_tumor12 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-AGLMAZ_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-BAJUW7_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-E9NNVR_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-FYVYBF_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-JALDBT_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-NHGSH5_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-U8ETE8_tumor61 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-UFL1D9_tumor41 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-UPSW7P_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-VVXZ39_tumor11 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-X3Q2S2_tumor31 -d debug
all_in_one_pipeline -h promise -t RNAseq -e netMHCpan4_1 -s s0 -r S014-XTGWF2_tumor11 -d debug

all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s0 -r S014-PQRN8X_tumor41 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s0 -r S014-TC14F4_tumor41 -d debug
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s s0 -r S014-ZNWXFT_tumor11 -d debug


## 2 catch patients 
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s i4cf5 -r K26K-JS49GN_metastasis11 -d debug
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s i4cf5 -r K26K-JS49GN_metastasis21 -d debug

all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s snv -r H021-3W9ZP5_metastasis031



##########################
## K26K-P7V4JK-tumor21
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s s0 -r K26K-P7V4JK-tumor11 -d debug
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s s0 -r K26K-P7V4JK-tumor21 -d debug

all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s s1b -r K26K-P7V4JK-tumor21
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s s1b -r K26K-P7V4JK-tumor11 -d debug


all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s snv -r K26K-P7V4JK-tumor21
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s indel -r K26K-P7V4JK-tumor21
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s fusion -r K26K-P7V4JK-tumor21


all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s snv -r K26K-P7V4JK-tumor11 &
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s indel -r K26K-P7V4JK-tumor11 &
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s fusion -r K26K-P7V4JK-tumor11 &


all_in_one_pipeline -r K26K-P7V4JK-tumor21 -t RNAseq -e netMHCpan4_1 -s s8bi4cf4f5
all_in_one_pipeline -r K26K-P7V4JK-tumor11 -t RNAseq -e netMHCpan4_1 -s s8bi4cf4f5


## K1FZ3P
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s s0 -r K1FZ3P -d debug


all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s s1b -r K1FZ3P -d debug

all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s snv -r K1FZ3P &
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s indel -r K1FZ3P &
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s fusion -r K1FZ3P &



## JS53 wait for odcf response
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s snv- -r JS53 -d debug

all_in_one_pipeline -t TCGA-PRAD -e netMHCpan4_1 -s s3 -r JS53
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s indel -r JS53 &
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s fusion -r JS53 &


## loh JS53
all_in_one_pipeline -t RNAseq -e netMHCpan4_1 -s loh -r JS53 -d debug

## loh on promise
all_in_one_pipeline -h promise -t luad_lusc -e netMHCpan4_1 -s loh -r S014-9C2UWH_tumor12 -d debug


## neo promise batch3 

## neo SF67 
all_in_one_pipeline -r SF67 -t RNAseq -e netMHCpan4_1 -d debug
	all_in_one_pipeline -r SF67 -t RNAseq -e netMHCpan4_1 -s s1b
	all_in_one_pipeline -r SF67 -t RNAseq -e netMHCpan4_1 -s s1b -d debug
	all_in_one_pipeline -r SF67 -t RNAseq -e netMHCpan4_1 -s snv &
	all_in_one_pipeline -r SF67 -t RNAseq -e netMHCpan4_1 -s indel &
	all_in_one_pipeline -r SF67 -t RNAseq -e netMHCpan4_1 -s fusion &
	all_in_one_pipeline -r SF67 -t RNAseq -e netMHCpan4_1 -d debug

all_in_one_pipeline -r SF67 -t RNAseq -e netMHCpan4_1 -s cn -d debug


## 98140 : MO79
all_in_one_pipeline -r MO79 -d debug
all_in_one_pipeline -r MO79 -s s1b -d debug
all_in_one_pipeline -r MO79 -s snv 
all_in_one_pipeline -r MO79 -s s8b -d debug

all_in_one_pipeline -r MO79 -s indel
all_in_one_pipeline -r MO79 -s f4f5 -d debug
all_in_one_pipeline -r MO79 -s loh -d debug

all_in_one_pipeline -r MO79 -s cn -d debug
all_in_one_pipeline -r MO79 -s loh -d debug


all_in_one_pipeline  -r MO79 -s blav -d debug
 
## LOH of K26K-P7V4JK
all_in_one_pipeline  -r K26K-P7V4JK-tumor21 -s loh
all_in_one_pipeline  -r K26K-P7V4JK-tumor11 -s loh



## MASTER cs86
all_in_one_pipeline -r CS86 -s s1b
all_in_one_pipeline -r CH80 -s s1b

all_in_one_pipeline -r CS86 -s snv
all_in_one_pipeline -r CH80 -s snv


## K26K-P7V4JK
all_in_one_pipeline -r  -d debug


## K1FZ3P
all_in_one_pipeline -r K1FZ3P -s s3 -d debug


#############

## K26K-TQ77H3-wes-metastasis11
all_in_one_pipeline -t TCGA-BRCA -r K26K-TQ77H3-wes-metastasis11 -d debug
all_in_one_pipeline -t TCGA-BRCA -r K26K-TQ77H3-wes-metastasis11 -s s1b 
all_in_one_pipeline -t TCGA-BRCA -r K26K-TQ77H3-wes-metastasis11 -s snv 
all_in_one_pipeline -t TCGA-BRCA -r K26K-TQ77H3-wes-metastasis11 -s indel 

all_in_one_pipeline -t TCGA-BRCA -r K26K-TQ77H3-wes-metastasis11 -s s3 -d debug

all_in_one_pipeline -t TCGA-BRCA -r K26K-TQ77H3-wes-metastasis11 -s xlsx -d debug
all_in_one_pipeline -t TCGA-BRCA -r K26K-TQ77H3-wes-metastasis11 -s nz -d debug



## H021-KQQR2Q-metastasis05
all_in_one_pipeline -r H021-KQQR2Q-metastasis05 -d debug
all_in_one_pipeline -r H021-KQQR2Q-metastasis05 -s s1b 

all_in_one_pipeline -r H021-KQQR2Q-metastasis05 -s snv &
all_in_one_pipeline -r H021-KQQR2Q-metastasis05 -s indel &
all_in_one_pipeline -r H021-KQQR2Q-metastasis05 -s fusion &
all_in_one_pipeline -r H021-KQQR2Q-metastasis05 -s xlsx -d debug
all_in_one_pipeline -r H021-KQQR2Q-metastasis05 -s nz -d debug



## H021-QXVLHV-metastasis04
all_in_one_pipeline -r H021-QXVLHV-metastasis04 -d debug
all_in_one_pipeline -r H021-QXVLHV-metastasis04 -s s1b

all_in_one_pipeline -r H021-QXVLHV-metastasis04 -s snv &
all_in_one_pipeline -r H021-QXVLHV-metastasis04 -s indel &
all_in_one_pipeline -r H021-QXVLHV-metastasis04 -s fusion &



## AW58
all_in_one_pipeline -r AW58 -d debug
all_in_one_pipeline -r AW58 -s s1b

## MW78
all_in_one_pipeline -r MW78 -d debug
all_in_one_pipeline -r MW78 -s s1b

all_in_one_pipeline -r MW78 -s xlsx -d debug

## MD48
all_in_one_pipeline -t TCGA-PAAD -r MD48 -d debug
all_in_one_pipeline -t TCGA-PAAD -r MD48 -s s3 -d debug
all_in_one_pipeline -r MD48 -s xlsx -d debug

## 237
all_in_one_pipeline -r  EB61  -t TCGA-OV -d debug

all_in_one_pipeline -r P980 -d debug
all_in_one_pipeline -r GH72 -d debug
all_in_one_pipeline -r MK55 -d debug
all_in_one_pipeline -r HB81 -d debug
all_in_one_pipeline -r SVB75 -d debug
all_in_one_pipeline -r AW58 -d debug


## K26K-P7V4JK meta1
all_in_one_pipeline -r K26K-P7V4JK_meta1  -d debug


# K9665
all_in_one_pipeline -r K9665  -d debug
all_in_one_pipeline -r K9665 -s s1b
all_in_one_pipeline -r K9665 -t TCGA-COAD -s s3s8as8b
all_in_one_pipeline -r K9665 -s indel
all_in_one_pipeline -r K9665 -s nz -d debug



# K9665_original_snv
all_in_one_pipeline -r K9665_originalSnv  -d debug
all_in_one_pipeline -r K9665_originalSnv -t TCGA-COAD -s snv 
all_in_one_pipeline -r K9665_originalSnv -t TCGA-COAD -s s8a -d debug


all_in_one_pipeline -r K9665 -t TCGA-COAD -s loh -d debug

all_in_one_pipeline -r K9665_originalSnv -t TCGA-COAD -s nz -d debug


# 
all_in_one_pipeline -r LS62  -d debug
all_in_one_pipeline -r LS62 -s snv 

## 
all_in_one_pipeline -r LM86  -d debug

## BAD60
all_in_one_pipeline -r BAD60 -t luad_lusc -s snv

## poc1 

all_in_one_pipeline -r EH63_cell-line01 -s nz -d debug
all_in_one_pipeline -r EH63_metastasis2 -s nz -d debug
all_in_one_pipeline -r EH63_metastasis3 -s nz -d debug
all_in_one_pipeline -r EH63_primary -s nz -d debug
all_in_one_pipeline -r EH63_tumor -s nz -d debug
all_in_one_pipeline -r EH63_xenograft1 -s nz -d debug
all_in_one_pipeline -r EH63_xenograft2 -s nz -d debug

## 
all_in_one_pipeline -r H021-UKY79T -s s1b
all_in_one_pipeline -r H021-RK6XZV_tumor -s s1b -d debug
all_in_one_pipeline -r H021-RK6XZV_tumor02 -s s1b

all_in_one_pipeline -r H021-RK6XZV_tumor04 -d debug

## poc003
all_in_one_pipeline -r poc003 -s s1b

## 6 NCT patients

all_in_one_pipeline -r AK74 -s xlsx_nz -d debug 
all_in_one_pipeline -r GB42 -s xlsx_nz -d debug 
all_in_one_pipeline -r KV65 -s xlsx_nz -d debug 
all_in_one_pipeline -r MF67 -s xlsx_nz -d debug 
all_in_one_pipeline -r MH73 -s xlsx_nz -d debug 
all_in_one_pipeline -r TS73_tumor02 -s xlsx_nz -d debug 
all_in_one_pipeline -r K26K-P7V4JK_meta01 -s xlsx_nz -d debug 

##

all_in_one_pipeline -r RMN50 -d debug 


## 2022.11.07
all_in_one_pipeline -r AL52 -d debug 
all_in_one_pipeline -r BK70 -d debug 
all_in_one_pipeline -r JM63 -d debug 

all_in_one_pipeline -r HSK48 -d debug 


###


# find ~/omics/EH63_metastasis2 ~/omics/EH63_primary ~/omics/EH63_tumor ~/omics/H021-S4CLSR ~/omics/EH63_cell-line01 ~/omics/EH63_metastasis3 ~/omics/EH63_xenograft1 ~/omics/EH63_xenograft2 ~/omics/H021-UKY79T ~/omics/poc003 ~/omics/H021-RK6XZV_tumor ~/omics/H021-RK6XZV_tumor02 ~/omics/H021-RK6XZV_tumor04 -type f | xargs -I {} ls -l {}

all_in_one_pipeline -r EH63_metastasis2 -t TCGA-READ -s xlsx-nz 
all_in_one_pipeline -r EH63_primary -t TCGA-READ -s xlsx-nz 
all_in_one_pipeline -r EH63_tumor -t TCGA-READ -s xlsx-nz 
all_in_one_pipeline -r H021-S4CLSR -s xlsx-nz 

all_in_one_pipeline -r EH63_cell-line01 -s xlsx-nz 
all_in_one_pipeline -r EH63_metastasis3 -s xlsx-nz 
all_in_one_pipeline -r EH63_xenograft1 -s xlsx-nz 
all_in_one_pipeline -r EH63_xenograft2 -s xlsx-nz 
all_in_one_pipeline -r H021-UKY79T -s xlsx-nz 
all_in_one_pipeline -r poc003 -s xlsx-nz 
all_in_one_pipeline -r H021-RK6XZV_tumor -s xlsx-nz 
all_in_one_pipeline -r H021-RK6XZV_tumor02 -s xlsx-nz 
all_in_one_pipeline -r H021-RK6XZV_tumor04 -s xlsx-nz 

##

all_in_one_pipeline -r BT60 -t TCGA-UCEC -d debug
all_in_one_pipeline -r GH71 -d debug
all_in_one_pipeline -r BS74 -d debug
all_in_one_pipeline -r WW46 -d debug

all_in_one_pipeline -r SH55 -d debug
all_in_one_pipeline -r SF82 -d debug

##

all_in_one_pipeline -r SH55_tumor -s fusion
all_in_one_pipeline -r SH55_tumor02 -s fusion
all_in_one_pipeline -r AJ63_tumor05 -d debug


all_in_one_pipeline -r TLS80 -t TCGA-COAD -d debug
##
all_in_one_pipeline -r PC63 -s s1b -d debug
all_in_one_pipeline -r SF82 -d debug


##
all_in_one_pipeline -r SaN65 -s fusion
all_in_one_pipeline -r AL52_tumor02 -s fusion
all_in_one_pipeline -r LF89_tumor02 -s fusion

all_in_one_pipeline -r AJ63_liver300X -d debug


## 

all_in_one_pipeline -r JB73 -s s1b -d debug
all_in_one_pipeline -r HO66 -s s1b -d debug

## 
all_in_one_pipeline -r LF89_tumor02 -s xlsx -d debug
all_in_one_pipeline -r LF89 -s xlsx -d debug

## 
all_in_one_pipeline -r JH59 -d debug
all_in_one_pipeline -r AA57 -d debug

all_in_one_pipeline -r AA57 -t TCGA-PAAD -s s8a

all_in_one_pipeline -r JH59 -s nz -t TCGA-BLCA -d debug
all_in_one_pipeline -r CA59 -s nz -t TCGA-PAAD -d debug

##
all_in_one_pipeline -r AW88 -t TCGA-COAD -s snv


## 

all_in_one_pipeline -h promise -r S014-2NHL7W_tumor11_wgs_RNAseq -d debug

## 
all_in_one_pipeline  -r K9665 -s loh -d debug
all_in_one_pipeline  -r K9665_tumor2 -s xlsx -d debug



## 
all_in_one_pipeline -r BK60 -d debug

all_in_one_pipeline -r TS73_tumor03 -t TCGA-COAD -s snv



##  fusion for ca59
all_in_one_pipeline -r CA59 -s f2 -d debug
all_in_one_pipeline -r JH59 -s f2 -d debug

## ts73

all_in_one_pipeline -r TS73_tumor03 -s xlsx -d debug

all_in_one_pipeline -r TLS80_tumor2 -t TCGA-COAD -d debug

## 
all_in_one_pipeline -r uv54 -t TCGA-PRAD -d debug -s snv


## 

all_in_one_pipeline -r RK66 -s xlsx -d debug

##
all_in_one_pipeline -r VM73 -s snv

## 
all_in_one_pipeline -r ES60 -t TCGA-COAD -s s3 
all_in_one_pipeline -r ES60 -t TCGA-COAD -s s1b 


all_in_one_pipeline -r MT68 -s s3 -d debug


## 

all_in_one_pipeline -r LP91 -s s3 -t TCGA-TGCT -d debug



## 

all_in_one_pipeline -r HDK59  -s loh
all_in_one_pipeline -r NH72  -s loh
all_in_one_pipeline -r K26K-ZWJDB4 -s loh
all_in_one_pipeline -r WW46_tumor02 -t TCGA-SKCM -s loh
all_in_one_pipeline -r ET76 -t TCGA-BRCA -s loh

#
all_in_one_pipeline -r LP91 -s i4a -d debug
all_in_one_pipeline -r LP91 -s s2 -d debug
all_in_one_pipeline -r HDK59  -s i4a -d debug


## 
5 odcf process
2 vcf process
arriba on vm 
immuX
hipo34 paper
gene_name_change_list
ask sven about fusion duel name


download 115024.21/22
upload   115024.21/22; 118612


## 

all_in_one_pipeline -t RNAseq    -c pa   -r GW46  -s snv-indel-fusion-xlsx
all_in_one_pipeline -t TCGA-BRCA -c brca -r RVK75 -s snv-indel-fusion-xlsx
all_in_one_pipeline -t RNAseq    -c PL   -r SMI80 -s snv-indel-fusion-xlsx
all_in_one_pipeline -t RNAseq    -c R    -r MTK72 -s snv-indel-fusion-xlsx
all_in_one_pipeline -t RNAseq    -c sg   -r SH82  -s snv-indel-fusion-xlsx

## 

all_in_one_pipeline -t TCGA-BRCA -r RVK75 -c brca -s s2 -d debug

## 
all_in_one_pipeline -t TCGA-READ -r JG48 -d debug
all_in_one_pipeline -t TCGA-LUAD -r NKH82 -d debug


##

all_in_one_pipeline  -r HaK67 -d debug
all_in_one_pipeline  -r GS42 -d debug -s s2 -d debug


all_in_one_pipeline -r BK70 -c brca -s c1-c2 -d debug

##

all_in_one_pipeline -r KW81         -t TCGA-MESO_RNAseq -c meso -s c1 

## 
all_in_one_pipeline -r IE67         -t TCGA-COAD -c coc -s s0 -d debug
all_in_one_pipeline -r IE67         -t TCGA-COAD -c coc -s xlsx -b Yes

all_in_one_pipeline -r JAS62        -t TCGA-PAAD -c pa -s s0 -b Yes -d debug
all_in_one_pipeline -r JAS62        -t TCGA-PAAD -c pa -s snv -b Yes

all_in_one_pipeline -r MO64 -s f5 -d debug



all_in_one_pipeline -r MB65  -t TCGA-PAAD -c paad  -b yes -s s1b -d debug
all_in_one_pipeline -r MN80  -t TCGA-PAAD -c paad  -b yes -s s1b -d debug
all_in_one_pipeline -r IW62  -t TCGA-PAAD -c paad  -b yes -s s1b -d debug
all_in_one_pipeline -r MK66  -t TCGA-READ -c read  -b yes -s s1b -d debug
all_in_one_pipeline -r cDA75 -t TCGA-OV   -c ov    -b yes -s s1b -d debug
all_in_one_pipeline -r FS69  -t TCGA-LUSC -c nsclc -b yes -s s1b -d debug
## 
all_in_one_pipeline -r H021-3TZJZW -t TCGA-LUSC_RNAseq -c uvm -s s0 -d debug 

##
all_in_one_pipeline -r K50K-QEQ5DD -t TCGA-BRCA_RNAseq -c brca -s xlsx

## 

all_in_one_pipeline -r P122886 -t TCGA-COAD_RNAseq -c coc -s s0 -d debug


all_in_one_pipeline -r MR80 -t TCGA-COAD -c coc -s s2 -b Yes -d debug

all_in_one_pipeline -r YP81 -t TCGA-COAD -c coc -s s2 -b Yes -d debug

all_in_one_pipeline -r 29056 -t TCGA-KIRP -c r -s indel -b Yes
all_in_one_pipeline -r 19714 -t TCGA-PAAD -c paad -s indel -b Yes


## 
all_in_one_pipeline -r P124838 -t TCGA-PAAD -c paad -s s0 -d debug
all_in_one_pipeline -r BM44 -t TCGA-READ_RNAseq -c read -s s2 -d debug
