#num=("2" "3" "4" "5" "6" "7" "8")
num=( "15" "25")
#num=("6")
#for num_wtag_int in ${num[@]}
#do
#num_wtag=0.$num_wtag_int
#name_wtag=0p$num_wtag_int
idx=-1; unset img;
num_wtag=0.6
num_rtag=0.6
name_wtag=0p6
name_rtag=0p6
python Makeplots_WWW_W.py --channel had --REGION SR1 --wtag ${num_wtag} --rtag ${num_rtag} --name MJJ__${name_wtag}__${name_rtag}
python Makeplots_WWW_W.py --channel had --REGION SR2 --wtag ${num_wtag} --rtag ${num_rtag} --name MJJ__${name_wtag}__${name_rtag}
python Makeplots_WWW_W.py --channel had --REGION SR3 --wtag ${num_wtag} --rtag ${num_rtag} --name MJJ__${name_wtag}__${name_rtag}


num_wtag=0.3
num_rtag=0.6
name_wtag=0p3
name_rtag=0p6




python Makeplots_WWW_W.py --channel had --REGION SR4 --wtag ${num_wtag} --rtag ${num_rtag} --name MJJ__${name_wtag}__${name_rtag}




num_wtag=0.6
num_rtag=0.6
name_wtag=0p6
name_rtag=0p6



python Makeplots_WWW_W.py --channel had --REGION SR5 --wtag ${num_wtag} --rtag ${num_rtag} --name MJJ__${name_wtag}__${name_rtag}


num_wtag=0.6
num_rtag=0.6
name_wtag=0p6
name_rtag=0p6




#idx=-1; unset img;
#((idx++)); img[idx]='had_SR2_/MJJ__'${name_wtag}__${name_rtag}'.png';
((idx++)); img[idx]='had_SR1_/MJJ__'${name_wtag}__${name_rtag}'.png';
((idx++)); img[idx]='had_SR2_/MJJ__'${name_wtag}__${name_rtag}'.png';
((idx++)); img[idx]='had_SR3_/MJJ__'${name_wtag}__${name_rtag}'.png';

num_wtag=0.3
num_rtag=0.6
name_wtag=0p3
name_rtag=0p6


((idx++)); img[idx]='had_SR4_/MJJ__'${name_wtag}__${name_rtag}'.png';


num_wtag=0.6
num_rtag=0.6
name_wtag=0p6
name_rtag=0p6



((idx++)); img[idx]='had_SR5_/MJJ__'${name_wtag}__${name_rtag}'.png';
mkdir -p R_fixed_W_changes9
montage -mode concatenate -tile 3x2 ${img[*]} R_fixed_W_changes9/M_alljets__${name_wtag}__${name_rtag}.png;
#done


