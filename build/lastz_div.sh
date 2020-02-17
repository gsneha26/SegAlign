rm run_lastz_ydrop.sh
for i in `seq 1 $1`;
do
    echo $i;
    echo "./lastz ~/WGA_GPU/data/hg38_chr1.fa ~/WGA_GPU/data/mm10_chr1.fa --segments=tmp"$i".segments --format=maf > gpu_out"$i".maf" >> run_lastz_ydrop.sh
done
