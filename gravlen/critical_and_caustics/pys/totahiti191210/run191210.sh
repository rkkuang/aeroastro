raynum=30000
imgsize=2048
massratio=(0.1 0.5)
base=1
X=(0.1 0.2 0.3 0.5 0.7 0.9 1 1.2)
#1.5
for q in ${massratio[@]}
do
for i in ${X[@]}
do
    echo $i
    if (( $(echo "$i < $base" |bc -l) )); then
        python3 two_pmnegq_inverse_ray_shoot_191210.py 4.5 $q $i $((imgsize)) $((raynum)) $i slow &
    else
        # echo $i "> 1"
        python3 two_pmnegq_inverse_ray_shoot_191210.py 4.5 $q $i $((imgsize)) $((raynum)) 1 slow &
    fi

done
python3 two_pmnegq_inverse_ray_shoot_191210.py 4.5 $q 1.5 $((imgsize)) $((raynum)) 1 slow
done


massratio=(0.9)
for q in ${massratio[@]}
do
for i in ${X[@]}
do
    echo $i
    if (( $(echo "$i < $base" |bc -l) )); then
        python3 two_pmnegq_inverse_ray_shoot_191210.py 6 $q $i $((imgsize)) $((raynum)) $i slow &
    else
        python3 two_pmnegq_inverse_ray_shoot_191210.py 6 $q $i $((imgsize)) $((raynum)) 1 slow &
    fi
done
python3 two_pmnegq_inverse_ray_shoot_191210.py 6 $q 1.5 $((imgsize)) $((raynum)) 1 slow
done







# for i in ${X[@]}
# do
#     echo $i
#     if [ $i -le  0.3]
#     then
#     python3 two_pmnegq_inverse_ray_shoot_191206.py 4.5 ${massratio} $i $((imgsize)) $((raynum)) 0.35 slow
#     elif [ $i -le  0.7]
#     python3 two_pmnegq_inverse_ray_shoot_191206.py 4.5 ${massratio} $i $((imgsize)) $((raynum)) 0.5 slow
# else 
#     python3 two_pmnegq_inverse_ray_shoot_191206.py 4.5 ${massratio} $i $((imgsize)) $((raynum)) $i slow
# done


# if [ $1 -gt $2 ]
# then echo “$1>$2”
# elif ....; then
# else echo “$2>$1”
# fi
#-gt是大于的意思
#-lt是小于
#-eq是等于
#-ne是不等于
#-ge是大于等于
#le是小于等于
