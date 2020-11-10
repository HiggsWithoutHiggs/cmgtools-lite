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
T2L="-P {ORIGIN}/NanoTrees_HWH_2lskim_170920/{YEAR} --Fs '{{P}}/3_recleaner_vMarcAllVars' --Fs '{{P}}/5_evtVars_v6' ".format(ORIGIN=ORIGIN, YEAR=YEAR)
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

FUNCTION_2L="hwh_catIndex_2lss_MVA(LepClean_Recl1_pdgId,LepClean_Recl2_pdgId,DNN_2lss_predictions_tjvv,DNN_2lss_predictions_tt,DNN_2lss_predictions_other)"
FUNCTION_3L="hwh_catIndex_3l_MVA(DNN_3l_predictions_tjvv,DNN_3l_predictions_tt,DNN_3l_predictions_other,LepClean_Recl1_pdgId,LepClean_Recl2_pdgId,LepClean_Recl3_pdgId,mZ1)"
FUNCTION_4L=''' "hwh_catIndex_4l_MVA(DNN_4l_predictions_tjvv,DNN_4l_predictions_tt,DNN_4l_predictions_other)" "[0.5,1.5,2.5]" '''
ONEBIN="1 1,0.5,1.5"
MCASUFFIX="mc"

DOFILE = ""

if REGION == "2lss":
    OPT_2L='{T2L} {OPTIONS} -W "puWeight"'.format(T2L=T2L, OPTIONS=OPTIONS)
    CATPOSTFIX=""

    CATBINS="[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5,35.5,36.5,37.5,38.5,39.5,40.5,41.5,42.5,43.5,44.5,45.5,46.5,47.5,48.5,49.5,50.5,51.5,52.5,53.5,54.5,55.5,56.5,57.5,58.5,59.5,60.5,61.5,62.5,63.5,64.5,65.5,66.5,67.5,68.5,69.5,70.5,71.5]"
    RANGES = '[1,6,12,25,44,57,72]'
    NAMES  = ','.join( '%s_%s'%(x,YEAR) for x in 'ee_tjVVnode,ee_Othernode,em_tjVVnode,em_Othernode,mm_tjVVnode,mm_Othernode'.split(','))
            
    TORUN='''python {SCRIPT} {DOFILE} hwh/mca-2l-{MCASUFFIX}{MCAOPTION}.txt hwh/2lss_tight.txt "{FUNCTION_2L}" "{CATBINS}" {SYSTS} {OPT_2L} --binname hwh_2lss --year {YEAR} --categorize-by-ranges "{RANGES}" "{NAMES}" '''.format(SCRIPT=SCRIPT, DOFILE=DOFILE, MCASUFFIX=MCASUFFIX, MCAOPTION=MCAOPTION, FUNCTION_2L=FUNCTION_2L, CATBINS=CATBINS, SYSTS=SYSTS, OPT_2L=OPT_2L,YEAR=YEAR,RANGES=RANGES,NAMES=NAMES)
    print submit.format(command=TORUN)
            

if REGION == "3l":
    OPT_3L='{T3L} {OPTIONS} -W "puWeight"'.format(T3L=T3L, OPTIONS=OPTIONS)
    CATPOSTFIX=""
    CATBINS="[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5,35.5,36.5,37.5,38.5,39.5,40.5,41.5,42.5]"
    RANGES = '[1,8,15,22,29,36,43]'
    NAMES  = ','.join( '%s_%s'%(x,YEAR) for x in ['tjVV_nr','tjVV_r','eee_Bkgnode','eem_Bkgnode','emm_Bkgnode','mmm_Bkgnode'])
    TORUN = 'python {SCRIPT} {DOFILE} hwh/mca-3l-{MCASUFFIX}{MCAOPTION}.txt hwh/3l_tight.txt "{FUNCTION_3L}" "{CATBINS}" {SYSTS} {OPT_3L} --binname hwh_3l --year {YEAR} --categorize-by-ranges "{RANGES}" "{NAMES}"'.format(SCRIPT=SCRIPT, DOFILE=DOFILE,MCASUFFIX=MCASUFFIX,MCAOPTION=MCAOPTION,FUNCTION_3L=FUNCTION_3L,CATBINS=CATBINS,YEAR=YEAR, SYSTS=SYSTS, OPT_3L=OPT_3L,RANGES=RANGES, NAMES=NAMES)
    print submit.format(command=TORUN)


if REGION=="4l": 
    OPT_4L='{T4L} {OPTIONS} -W "puWeight"'.format(T4L=T4L,OPTIONS=OPTIONS)
    CATPOSTFIX=""
    TORUN="python {SCRIPT} {DOFILE} hwh/mca-4l-{MCASUFFIX}{MCAOPTION}.txt hwh/4l_tight.txt {FUNCTION_4L} {SYSTS} {OPT_4L} --binname hwh_4l --year {YEAR} ".format(SCRIPT=SCRIPT, DOFILE=DOFILE,MCASUFFIX=MCASUFFIX,MCAOPTION=MCAOPTION,FUNCTION_4L=FUNCTION_4L,SYSTS=SYSTS,OPT_4L=OPT_4L,YEAR=YEAR)
    print submit.format(command=TORUN)
