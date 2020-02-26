cd $HOME/WGA_GPU/data1
ref=$(basename -s .fa $1)
mkdir $ref
faSplit byname $ref.fa $ref/
cd $ref
for i in *;
do
    echo "faToTwoBit $i $(basename -s .fa $i).2bit" >> cmd.sh
done
lines=$(wc -l cmd.sh | cut -d " " -f 1)
parallel --jobs=$lines < cmd.sh
rm cmd.sh

cd $HOME/WGA_GPU/data1
query=$(basename -s .fa $2)
mkdir $query
faSplit byname $query.fa $query/
cd $query
for i in *;
do
    echo "faToTwoBit $i $(basename -s .fa $i).2bit" >> cmd.sh
done
lines=$(wc -l cmd.sh | cut -d " " -f 1)
parallel --jobs=$lines < cmd.sh
rm cmd.sh
