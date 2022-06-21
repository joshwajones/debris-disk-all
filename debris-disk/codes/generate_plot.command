#! /bin/bash
# takes in file name and generates alt/az images
cd "${0%/*}"
read name
source /Users/sjosh/opt/anaconda3/etc/profile.d/conda.sh
conda activate img_2 
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/gen_debris_disk_single.py $name
for az in 0 45 90 135 180 
do
	for alt in 90 45 10 05 00 
	do 
		Echo "Generating image, ALT: $alt, AZ: $az"
		python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name $alt $az 1
	done
done