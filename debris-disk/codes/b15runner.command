#! /bin/bash
# takes in file name and generates images at 0, 5, 10, 45, 90 degree altitudes 
cd "${0%/*}"
read name
source /Users/sjosh/opt/anaconda3/etc/profile.d/conda.sh
conda activate img_2 
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/gen_debris_disk_single.py $name
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 00 00 1
 python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 90 00 1
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 05 00 1
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 10 00 1
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 45 00 1