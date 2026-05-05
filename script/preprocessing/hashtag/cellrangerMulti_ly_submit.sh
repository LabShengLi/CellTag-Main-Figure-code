
sbatch -p epyc-64 --job-name=cr_multi_O_vitro cellrangerMulti_ly.sh O_vitro

exit 0

cd /project2/sli68423_1316/users/yang/workspace/U01_aim2/scripts/yang/preprocessing/cellranger/Cell_Ranger_Multi

sbatch -p largemem --job-name=cr_multi_OO cellrangerMulti_ly.sh OO
sbatch -p largemem --job-name=cr_multi_OY cellrangerMulti_ly.sh OY
sbatch -p largemem --job-name=cr_multi_Y_vitro cellrangerMulti_ly.sh Y_Vitro
sbatch -p largemem --job-name=cr_multi_YY cellrangerMulti_ly.sh YY
sbatch -p largemem --job-name=cr_multi_YO  cellrangerMulti_ly.sh YO

