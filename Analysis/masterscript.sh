#!/bin/bash
myArray=(2a_createinst.R 2b_instwinkler.R 2bb_reswinkler.R 2c_NAFLD5e6.R 5_harmonise_and_results.R 6a_createplots.R 7b_createtables.R 8_textresultssection.R)
for str in ${myArray[@]}; do
chmod u+x ./$str
done
echo 'Initializing 2a_createinst.R' && ./2a_createinst.R && echo 'Initializing 2b_instwinkler.R' && ./2b_instwinkler.R && echo 'Initializing 2bb_reswinkler.R' && ./2bb_reswinkler.R && echo 'Initializing 2c_NAFLD5e6.R' && ./2c_NAFLD5e6.R && echo 'Initializing 5_harmonise_and_results.R' && ./5_harmonise_and_results.R && echo 'Initializing 6a_createplots.R' && ./6a_createplots.R && echo 'Initializing 7b_createtables.R' && ./7b_createtables.R && echo 'Initializing 8_textresultssection.R' && ./8_textresultssection.R && echo 'The master script finished without errors'
