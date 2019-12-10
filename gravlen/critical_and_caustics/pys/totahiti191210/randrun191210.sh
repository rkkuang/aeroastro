xylim=3
Nrows=5
posscale=1
imgsz=2048
raynum=30
kscale=0.6
msrang0=-1
msrang1=1
# run 6 * 8 = 48 times
for (( i=0; i<1; ++i)); do
for (( j=0; j<7; ++j)); do
python3 n_pointM_lens_inverse_ray_shoot.py $((xylim)) $((Nrows)) $((posscale)) $((imgsz)) $((raynum)) $kscale slow $((msrang0)) $((msrang1)) &
done
python3 n_pointM_lens_inverse_ray_shoot.py $((xylim)) $((Nrows)) $((posscale)) $((imgsz)) $((raynum)) $kscale slow $((msrang0)) $((msrang1))
done