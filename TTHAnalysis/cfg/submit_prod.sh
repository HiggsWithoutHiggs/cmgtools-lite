for Y in 2018 ; do
     nanopy_batch.py -o ${Y} -r /store/cmst3/group/wmass/secret/NanoTrees_HWH_2lskim_170920/${Y} run_ttH_fromNanoAOD_cfg.py --option analysis=main  --option year=$Y --option selectComponents=HWHSIGNAL -b 'run_condor_simple.sh -t 1200 -a group_u_CMST3.all  ./batchScript.sh' -B;  echo sleep;  : sleep 240m
done
