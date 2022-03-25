#!/bin/bash

declare -a sp_names=("Nasutitermes_similis"
"Odontotermes_hainanensis"
"Reticulitermes_virginicus"
"Microhodotermes_viator"
"Macrotermes_malaccensis"
"Nasutitermes_graveolus"
"Macrotermes_gilvus"
"Reticulitermes_aculabialis"
"Zootermopsis_nevadensis"
"Nasutitermes_exitiosus"
"Nasutitermes_corniger"
"Macrotermes_falciger"
"Termes_fatalis"
"Tumulitermes_pastinator"
"Reticulitermes_santonensis"
"Labiotermes_labralis"
"Macrotermes_carbonarius"
"Nasutitermes_triodiae"
"Reticulitermes_flaviceps"
"Nasutitermes_octopilis"
"Reticulitermes_chinensis"
"Coptotermes_lacteus"
"Macrotermes_subhyalinus"
"Pseudacanthotermes_spiniger"
"Macrotermes_natalensis"
"Reticulitermes_labralis"
"Reticulitermes_leptomandibularis"
"Mastotermes_darwiniensis"
"Microcerotermes_crassus"
"Cavitermes_tuberosus"
"Reticulitermes_hageni"
"Zootermopsis_angusticollis"
"Reticulitermes_flavipes"
"Coptotermes_formosanus"
"Embiratermes_brevinasus"
"Nasutitermes_longipennis"
"Embiratermes_neotenicus")

for sp_name in ${sp_names[@]}; do

PATH_TO_IN="/home/glebo/Documents/lab/TermitesMutSpectrum/Body/2Derived/Cockroaches/binom_subs_per_sp/$sp_name.csv"
PATH_TO_OUT="/home/glebo/Documents/lab/TermitesMutSpectrum/Body/3Results/binom_results/$sp_name.compared.csv" 
python3 binom.py $PATH_TO_IN $PATH_TO_OUT

done



