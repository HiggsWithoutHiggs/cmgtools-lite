for Y in 2018 ; do
     nanopy_batch.py -o ${Y} -r /store/group/phys_higgs/emanuele/higgs_wo_higgs/NanoTrees_HWH_100720_v9/${Y} run_HwH_fromNanoAOD_cfg.py --option analysis=main  --option year=$Y -b 'run_condor_simple.sh -t 1200 -a group_u_CMST3.all  ./batchScript.sh' -B;  echo sleep;  : sleep 240m
done
