# 
# This script executes script GF_standalone (compiles and execute GF) and then compare with base results.
#

# MONAN
#
# Author: Denis Eiras
#
# E-mail: denis.eiras@inpe.br  
#
# Date: 2023/03/30
#
# Version: TODO -> USAR VERSÃO DA RELEASE, QUANDO FOR GERADA
#
#
# Full description
# ================
# - executes GF_standalone.sh script
# - run binay comparison between grads .gra files
# - run grads script to compare all variables of all levels of files that fails binary comparison.
# 
#
#
# Instructions
# ~~~~~~~~~~~~
# 
# - First of all (once):
#   - The base results are reference files generated before for comparison. Parametrize then in "DIR_REF" variable. 
#
# - Now you can execute this script:
#   $ ./compare_local_modifications.sh
#
#
# History
# =======
# 2023/03/60 - Creation - deniseiras
# 2023/05/11 - Update - ncolumns; using script compare_diff_and_contourn_all_variables
#
#
# Licence
# =======
#
# <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the  Free  Software  Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it  will be useful, but
# ** WITHOUT  ANY  WARRANTY **;  without  even  the   implied   warranty  of
# **MERCHANTABILITY** or **FITNESS FOR A  PARTICULAR PURPOSE**.  See  the, GNU
# GNU General Public License for more details.
#
# You should have received a copy  of the GNU General  Public  License
# along with this program.  If not, see [GNU Public    License](https://www.gnu.org/licenses/gpl-3.0.html).
#

DIR_TESTING="../dataout"
DIR_REF="../refs"

echo -e "\n\nCompiling and executing GF"
echo -e "==============================\n"

./GF_standalone.sh
echo -e "\n\nComparing binary diferences"
echo -e "===========================\n"


diff -qr $DIR_TESTING $DIR_REF > diff.txt
sed -i '/^Somente/d' diff.txt
if [ -s diff.txt ]; then

  while read line
  do
    echo
    echo -e "\n\nComparing using grads"
    echo -e "=====================\n"
    file_path_test=$(echo $line | awk '{print $3}')
    file_path_test=$(echo $file_path_test | sed 's/gra/ctl/g')

    file_path_ref=$(echo $line | awk '{print $5}')
    file_path_ref=$(echo $file_path_ref | sed 's/gra/ctl/g')

    file_name=$(basename $file_path_test)
    file_base_name="${file_name%.*}"
    
        
    GRADS_PARAMS='grads -lc "compare_diff_and_contourn_all_variables.gs '$file_path_test' '$file_path_ref' '$file_base_name' "'
    echo $GRADS_PARAMS
    eval $GRADS_PARAMS
    
  done < diff.txt

else
  echo -e "\nCongratulations, no diferences found!!!"
fi




