#! /bin/bash
# takes in file name and generates alt/az images
cd "${0%/*}"
source /Users/sjosh/opt/anaconda3/etc/profile.d/conda.sh
conda activate img_2 

python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/gen_debris_disk_single.py animation >& gentest.out
for i in {0..10}
do 
	let alt=-2*i 
	python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 animation $alt 90 1 100
done 