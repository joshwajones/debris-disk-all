#! /bin/bash
# takes in file name and generates alt/az images
cd "${0%/*}"
read name
source /Users/sjosh/opt/anaconda3/etc/profile.d/conda.sh
conda activate img_2 
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/gen_debris_disk_single.py $name
Echo "Generating image, ALT: 90, AZ: 00"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 90 00 1
Echo "Generating image, ALT: 45, AZ: 00"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 45 00 1
Echo "Generating image, ALT: 10, AZ: 00"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 10 00 1
Echo "Generating image, ALT: 05, AZ: 00"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 05 00 1
Echo "Generating image, ALT: 00, AZ: 00"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 00 00 1
Echo "Generating image, ALT: 90, AZ: 45"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 90 45 1
Echo "Generating image, ALT: 45, AZ: 45"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 45 45 1
Echo "Generating image, ALT: 10, AZ: 45"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 10 45 1
Echo "Generating image, ALT: 05, AZ: 45"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 05 45 1
Echo "Generating image, ALT: 00, AZ: 45"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 00 45 1
Echo "Generating image, ALT: 90, AZ: 90"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 90 90 1
Echo "Generating image, ALT: 45, AZ: 90"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 45 90 1
Echo "Generating image, ALT: 10, AZ: 90"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 10 90 1
Echo "Generating image, ALT: 05, AZ: 90"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 05 90 1
Echo "Generating image, ALT: 00, AZ: 90"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 00 90 1
Echo "Generating image, ALT: 90, AZ: 135"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 90 135 1
Echo "Generating image, ALT: 00, AZ: 135"
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 45 135 1
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 10 135 1
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 05 135 1
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 00 135 1
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 90 180 1
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 45 180 1
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 10 180 1
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 05 180 1
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 $name 00 180 1
