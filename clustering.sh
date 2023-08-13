#! /bin/bash

rm -rf /data/underworld/kbas/03_data/raw_dif &
rm -rf /data/underworld/kbas/03_data/source_dif &
rm -rf /data/underworld/kbas/03_data/derivatives_dif &
rm -rf /data/underworld/kbas/03_data/processed_dif &
wait
cp -R /data/underworld/kbas/03_data/raw/. /data/underworld/kbas/03_data/raw_dif &
cp -R /data/underworld/kbas/03_data/source/. /data/underworld/kbas/03_data/source_dif &
cp -R /data/underworld/kbas/03_data/derivatives/. /data/underworld/kbas/03_data/derivatives_dif &
mkdir /data/underworld/kbas/03_data/processed_dif
cd /data/underworld/kbas/03_data/processed_dif
mkdir 112111 128221 130519 170192 176117 208010 211787 214685 232237 308597 324038 330406 346878
cd /data/underworld/kbas/clustering
wait 

matlab -batch "mpm_registration_to_b0_multi" &

matlab -batch "between_subject_transformation" &

wait

python /data/underworld/kbas/clustering/estimate.py &

matlab -batch "test" &

wait

matlab -batch "parameter_map_registration_to_b0_multi_v3" &

source ./data/underworld/kbas/clustering/tractography_multisubject &

wait

matlab -batch "labelling_diffusion_space" &

matlab -batch /data/underworld/kbas/clustering/deformation_reslicing_3.m &

wait

matlab -batch "parameter_map_masking"

matlab -batch /data/underworld/kbas/clustering/klara_em_multi_dif_qmri.m
