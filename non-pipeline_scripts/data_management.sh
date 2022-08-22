# working dir
/mnt/users/ngda/proj/paper1
# attac dir
#cd /mnt/project/Aqua-Faang/seq_data/AtlanticSalmon/ATAC

############ WGS
cd /mnt/users/ngda/ngs_data/atlantic_salmon/wgs

bcftools +tag2tag in.vcf -- -r --pl-to-gl

for i in 1..29
do
    bcftools index ssa${i}.DP10.GQ10.MS0.7.recode.vcf.gz &
done

for i in 1..29
do
    bcftools index ssa${i}.DP10.GQ10.MS0.7.recode.vcf.gz
    #bcftools +tag2tag ssa${i}.DP10.GQ10.MS0.7.recode.vcf.gz -r --pl-to-gl -Oz -o ssa${i}.DP10.GQ10.MS0.7.recode_GL.vcf.gz
done


i=29
bcftools +tag2tag ssa${i}.DP10.GQ10.MS0.7.recode.vcf.gz -- -r --pl-to-gl | bgzip > ssa${i}.recode_GL.vcf.gz
















########### ATAC SEQ
## tissues code:
	# tissueCode <- c(
	# 	"1" = "Brain",
	# 	"2" = "Gill",
	# 	"3" = "Gonad",
	# 	"4" = "HeadKidney",
	# 	"5" = "Liver",
	# 	"6" = "Muscle",
	# 	"7" = "SuppDistalIntestine",
	# 	"8" = "Sperm"
	# )
#

cd /mnt/project/Aqua-Faang/seq_data/AtlanticSalmon/ATAC
home_dir=/mnt/users/ngda/ngs_data/atlantic_salmon/atac

# 1. brain
mkdir ${home_dir}/brain
for fn in C11????1_?.fq.gz
do
    ln -s /mnt/project/Aqua-Faang/seq_data/AtlanticSalmon/ATAC/${fn} ${home_dir}/brain/${fn}
done

# 2. gill

mkdir ${home_dir}/gill
for fn in C11????2_?.fq.gz
do
    ln -s /mnt/project/Aqua-Faang/seq_data/AtlanticSalmon/ATAC/${fn} ${home_dir}/gill/${fn}
done

