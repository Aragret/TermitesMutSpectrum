#!/bin/bash

declare -a sp_names=("Blaberus_craniifer" "Blattella_germanica" "Blattella_germanica" "Cryptocercus_primarius" "Cryptocercus_punctulatus" "Geoscapheus_robustus" "Macropanesthia_lithgowae" "Macropanesthia_rhinoceros" "Panchlora_nivea" "Panesthia_cribrata" "Panesthia_sloanei" "Panesthia_tryoni_tryoni" "Paratemnopteryx_howarthi" "Paratemnopteryx_stonei" "Periplaneta_americana" "Periplaneta_australasiae" "Pycnoscelus_surinamensis" "Salganea_aequaliterspinosa" "Salganea_guentheri" "Salganea_raggei")

for sp_name in ${sp_names[@]}; do

PATH_TO_IN="/home/glebo/Documents/lab/TermitesMutSpectrum/Body/2Derived/Cockroaches/binom_subs_per_sp/$sp_name.csv"
PATH_TO_OUT="/home/glebo/Documents/lab/TermitesMutSpectrum/Body/3Results/binom_results/$sp_name.compared.csv" 
python3 binom.py $PATH_TO_IN $PATH_TO_OUT

done



