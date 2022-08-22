# working dir
/mnt/users/ngda/proj/paper1
# attac dir
#cd /mnt/project/Aqua-Faang/seq_data/AtlanticSalmon/ATAC
## link to home
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

