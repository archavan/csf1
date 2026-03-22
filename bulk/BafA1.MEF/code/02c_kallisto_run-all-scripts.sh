#!/usr/bin/env bash

cd ~/csf1/bulk/BafA1.MEF/code/02b_kallisto_scripts

sbatch BafA1.MEF.0h_1.sh
sbatch BafA1.MEF.0h_2.sh
sbatch BafA1.MEF.4h_1.sh
sbatch BafA1.MEF.4h_2.sh
sbatch BafA1.MEF.8h_1.sh
sbatch BafA1.MEF.8h_2.sh
sbatch BafA1.MEF.12h_1.sh
sbatch BafA1.MEF.12h_2.sh
sbatch BafA1.MEF.24h_1.sh
sbatch BafA1.MEF.24h_2.sh
sbatch BafA1.MEF.48h_1.sh
sbatch BafA1.MEF.48h_2.sh
