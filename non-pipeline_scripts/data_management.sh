# working dir
/mnt/users/ngda/proj/paper1
# attac dir
#cd /mnt/project/Aqua-Faang/seq_data/AtlanticSalmon/ATAC

############ WGS
cd /mnt/users/ngda/ngs_data/atlantic_salmon/wgs



########### ATAC BAM
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


# 1. brain
source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/Brain/results/bwa/mergedLibrary
cd $source_data

home_dir=/mnt/users/ngda/ngs_data/atlantic_salmon/atac

mkdir ${home_dir}/brain
mkdir ${home_dir}/brain/bam


for fn in *bam*
do
    out_fn=$(echo "$fn" | sed "s/sorted.//")
	ln -s ${source_data}/${fn} ${home_dir}/brain/bam/${out_fn}
done



# 2. gill



