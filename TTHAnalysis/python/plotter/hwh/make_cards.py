import os, sys
nCores=8
#submit = '''sbatch -c %d -p short  --wrap '{command}' '''%nCores
submit = '{command}' 


ORIGIN="/eos/cms/store/cmst3/group/wmass/secret/"; 

if len(sys.argv) < 4: 
    print 'Sytaxis is %s [outputdir] [year] [region] [other]'%sys.argv[0]
    raise RuntimeError 
OUTNAME=sys.argv[1]
YEAR=sys.argv[2]
REGION=sys.argv[3]
OTHER=sys.argv[4:] if len(sys.argv) > 4 else ''

if   YEAR in '2016': LUMI="35.9"
elif YEAR in '2017': LUMI="41.4"
#elif YEAR in '2018': LUMI="59.7"
elif YEAR in '2018': LUMI="137.0" # using 2018 MC/data as proxy for tot Run 2
else:
    raise RuntimeError("Wrong year %s"%YEAR)


#print "Normalizing to {LUMI}/fb".format(LUMI=LUMI);
OPTIONS=" --tree NanoAOD --s2v -j {J} -l {LUMI} -f --WA prescaleFromSkim --split-factor=-1 ".format(LUMI=LUMI,J=nCores)
os.system("test -d cards/{OUTNAME} || mkdir -p cards/{OUTNAME}".format(OUTNAME=OUTNAME))
OPTIONS="{OPTIONS} --od cards/{OUTNAME} ".format(OPTIONS=OPTIONS, OUTNAME=OUTNAME)
T2L="-P {ORIGIN}/NanoTrees_HWH_2lskim_170920/{YEAR} --Fs '{{P}}/3_recleaner_vMarcAllVars' --Fs '{{P}}/5_evtVars_v3' ".format(ORIGIN=ORIGIN, YEAR=YEAR)
T3L=T2L
T4L=T2L

SYSTS="--unc hwh/systsUnc.txt --amc "
MCAOPTION=""
#MCAOPTION="-splitdecays"
ASIMOV="--asimov 's+b'"
#ASIMOV="" 
SCRIPT= "makeShapeCardsNew.py"
PROMPTSUB=""
#PROMPTSUB="--plotgroup data_fakes+=.*_promptsub"

if 'unblind' in OTHER:
    ASIMOV=""

print "We are using the asimov dataset"
OPTIONS="{OPTIONS} -L ttH-multilepton/functionsTTH.cc -L hwh/functionsHWH.cc {PROMPTSUB} --neg   --threshold 0.01 {ASIMOV}  ".format(OPTIONS=OPTIONS,PROMPTSUB=PROMPTSUB,ASIMOV=ASIMOV) # neg necessary for subsequent rebin #  
CATPOSTFIX=""

FUNCTION_2L="hwh_catIndex_2lss_MVA(LepClean_Recl1_pdgId,LepClean_Recl2_pdgId,DNN_2lss_predictions_tjvv,DNN_2lss_predictions_ttV,DNN_2lss_predictions_tt,DNN_2lss_predictions_other)"
FUNCTION_3L="hwh_catIndex_3l_MVA(DNN_2lss_predictions_tjvv,DNN_2lss_predictions_ttV,DNN_2lss_predictions_tt,DNN_2lss_predictions_other,LepClean_Recl1_pdgId,LepClean_Recl2_pdgId,LepClean_Recl3_pdgId,mZ1)"

ONEBIN="1 1,0.5,1.5"
MCASUFFIX="mc"

DOFILE = ""

if REGION == "2lss":
    OPT_2L='{T2L} {OPTIONS} -W "puWeight"'.format(T2L=T2L, OPTIONS=OPTIONS)
    CATPOSTFIX=""

    CATBINS="[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5,35.5,36.5,37.5,38.5,39.5,40.5,41.5,42.5,43.5,44.5,45.5,46.5,47.5,48.5,49.5,50.5,51.5,52.5,53.5,54.5,55.5,56.5,57.5,58.5,59.5,60.5,61.5,62.5,63.5,64.5,65.5,66.5,67.5,68.5,69.5,70.5,71.5,72.5,73.5,74.5,75.5,76.5,77.5,78.5,79.5,80.5,81.5,82.5,83.5,84.5,85.5,86.5,87.5,88.5,89.5,90.5,91.5,92.5,93.5,94.5,95.5,96.5,97.5,98.5,99.5,100.5,101.5,102.5,103.5,104.5,105.5,106.5,107.5,108.5,109.5,110.5,111.5,112.5,113.5,114.5,115.5,116.5,117.5,118.5,119.5,120.5]"
    RANGES = '[1,6,14,20,24,37,45,64,75,88,99,114,121]'
    NAMES  = ','.join( '%s_%s'%(x,YEAR) for x in 'ee_tjVVnode,ee_Othernode,ee_ttVnode,ee_ttnode,em_tjVVnode,em_Othernode,em_ttVnode,em_ttnode,mm_tjVVnode,mm_Othernode,mm_ttVnode,mm_ttnode'.split(','))
            
    TORUN='''python {SCRIPT} {DOFILE} hwh/mca-2l-{MCASUFFIX}{MCAOPTION}.txt hwh/2lss_tight.txt "{FUNCTION_2L}" "{CATBINS}" {SYSTS} {OPT_2L} --binname hwh_2lss --year {YEAR} --categorize-by-ranges "{RANGES}" "{NAMES}" '''.format(SCRIPT=SCRIPT, DOFILE=DOFILE, MCASUFFIX=MCASUFFIX, MCAOPTION=MCAOPTION, FUNCTION_2L=FUNCTION_2L, CATBINS=CATBINS, SYSTS=SYSTS, OPT_2L=OPT_2L,YEAR=YEAR,RANGES=RANGES,NAMES=NAMES)
    print submit.format(command=TORUN)
            

if REGION == "3l":
    OPT_3L='{T3L} {OPTIONS} -W "puWeight"'.format(T3L=T3L, OPTIONS=OPTIONS)
    CATPOSTFIX=""
    CATBINS="[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5,35.5,36.5,37.5,38.5,39.5,40.5,41.5,42.5,43.5,44.5,45.5,46.5,47.5,48.5,49.5,50.5,51.5,52.5,53.5,54.5,55.5,56.5,57.5,58.5,59.5,60.5,61.5,62.5,63.5,64.5,65.5,66.5,67.5,68.5,69.5,70.5]"
    RANGES = '[1,8,15,22,29,36,43,50,57,64,71]'
    NAMES  = ','.join( '%s_%s'%(x,YEAR) for x in ['tjVV_nr','tjVV_r','eee_Othernode','eee_ttnode','eem_Othernode','eem_ttnode','emm_Othernode','emm_ttnode','mmm_Othernode','mmm_ttnode'])
    TORUN = 'python {SCRIPT} {DOFILE} hwh/mca-3l-{MCASUFFIX}{MCAOPTION}.txt hwh/3l_tight.txt "{FUNCTION_3L}" "{CATBINS}" {SYSTS} {OPT_3L} --binname hwh_3l --year {YEAR} --categorize-by-ranges "{RANGES}" "{NAMES}"'.format(SCRIPT=SCRIPT, DOFILE=DOFILE,MCASUFFIX=MCASUFFIX,MCAOPTION=MCAOPTION,FUNCTION_3L=FUNCTION_3L,CATBINS=CATBINS,YEAR=YEAR, SYSTS=SYSTS, OPT_3L=OPT_3L,RANGES=RANGES, NAMES=NAMES)
    print submit.format(command=TORUN)

