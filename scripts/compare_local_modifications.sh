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
# - run grads diff comparison of variables outt1, outt2 and outt3
# - run scripts
#
#
# Instructions
# ~~~~~~~~~~~~
# 
# - First of all (once):
#
#   - The base results are reference files generated before for comparison. Must be at ../dataout__compare_local_modifications. 
#     So, execute:
#          
#     $ ./GF_Standalone.sh
#     $ copy -rf ../dataout ../dataout__compare_local_modifications
#
# - Now you can execute this script:
#
#   $ ./compare_local_modifications.sh
#
#
# History
# =======
# 2023/03/60 - Creation - deniseiras
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

echo -e "\n\nCompiling and executing GF"
echo -e "==============================\n"

./GF_standalone.sh
echo -e "\n\nComparing binary diferences"
echo -e "===========================\n"
diff -qr ../dataout/ref_g.gra ../dataout__compare_local_modifications/ref_g.gra 

echo
echo -e "\n\nComparing using grads"
echo -e "=====================\n"

grads -lc "compare_diff_and_contourn_outt_variables.gs"