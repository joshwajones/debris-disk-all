#! /bin/bash
cd "${0%/*}"
source /Users/sjosh/opt/anaconda3/etc/profile.d/conda.sh
conda activate img_2 
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/gen_debris_disk_single.py background_24
python /Users/sjosh/pycharmprojects/research/img_2/debris-disk/codes/save_altaz.py 400 10 background_24 90 00 1