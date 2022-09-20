#! /bin/bash
# takes in file name and generates images at 3, 5, 10, 90 altitudes
cd "${0%/*}"
read name
source /Users/sjosh/opt/anaconda3/etc/profile.d/conda.sh
conda activate img_2 

python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/gen_debris_disk_single.py $name

Echo "Generating image: Alt -05, Az 90"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name -05 90 1 100 

 

Echo "Generating image: Alt -03, Az 90"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name -03 90 1 100

Echo "Generating image: Alt -90, Az 90"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name -90 90 1 100 


Echo "Generating image: Alt -10, Az 90"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name -10 90 1 100 

Echo "Generating image: Alt 00, Az 90"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 00 90 1 100

Echo "Generating image: Alt -45, Az 90"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name -45 90 1 100
