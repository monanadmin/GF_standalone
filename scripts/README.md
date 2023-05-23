# Scripts auxiliares para execução do gf

- GF_standalone.sh					                - Script que compila e executa o codigo GF_StandAlone localmente

- run.csh						                        - Executa o gf com uma configuração definida no script e compara com a referência

- compare_local_modifications.sh            - Executa o script GF_standalone e compara os resultados com uma base.
  - o script usa os arquivos:
    - compare_diff_and_contourn_all_variables.gs

- run_gf_times_with_hpctoolkit_or_gprof.sh  - executa o gf x vezes instrumentado com o hpc_toolkit ou grof (comentar / descomentar a parte do script para selecionar)
  - o script usa os arquivos: 
    - env_hpctoolkit.sh
    - submit_hpctoolkit_from_script_run_gf_times.sbatch

- submit.sh                                 - submete um job com o gf
