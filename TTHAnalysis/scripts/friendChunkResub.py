#!/usr/bin/env python
# USAGE: 
# 1. if you have already the output of friendChunkCheck.sh -z fdir > zombies.txt
# ./scripts/friendChunkResub.py frienddir maintreedir zombies.txt -N 500000
# 2. if you don't have it
#  ./scripts/friendChunkResub.py frienddir maintreedir --run-checker -N 500000

import sys,os,re

MODULES_DATA = ["recleaner_step1","recleaner_step2_data","triggerSequence"]
MODULES_MC = ["recleaner_step1","recleaner_step2_mc","mcMatch_seq","higgsDecay","triggerSequence"]

if __name__ == "__main__":
    
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] dir_with_friend_trees dir_with_trees [checkedfile.txt]")
    parser.add_option("-c",      "--run-checker", dest="runChecker",  action='store_true', default=False, help="Run the script friendChunkCheck.sh");
    parser.add_option("-N", "--events",  dest="chunkSize", type="int",    default=1000000, help="Default chunk size when splitting trees");
    parser.add_option(      "--data",  dest="data", action="store_true", default=False, help="Run the data modules (default is MC modules)");
    parser.add_option("-p", "--pretend", dest="pretend",   action="store_true", default=False, help="Don't run anything");
    parser.add_option("-q", "--queue",   dest="queue",     type="string", default=None, help="Run jobs on lxbatch queue or condor instead of locally");
    parser.add_option("--sub", "--subfile", dest="subfile", type="string", default="resub_condor.sub", help="Subfile for condor (default: condor.sub)");
    parser.add_option("--maxruntime", "--time",  dest="maxruntime", type="int", default=360, help="Condor job wall clock time in minutes (default: 6h)");
    (options, args) = parser.parse_args()

    fdir = args[0]
    maindir = args[1]

    MODULES = MODULES_DATA if options.data else MODULES_MC

    if options.runChecker:
        print "Running friendChunkCheck.py now on ",fdir,". Will take time..."
        if len(args)==3:
            print "Can't run the checker if you pass the output file of a previous check (for safety)"
            sys.exit(1)
        pyCheckScript = '{cmssw}/src/CMGTools/WMass/python/postprocessing/scripts/friendChunkCheck.py -z {frienddir} > tmpcheck.txt'.format(cmssw=os.environ['CMSSW_BASE'],
                                                                                                                                              frienddir=fdir)
        if os.path.isfile('tmpcheck.txt'):
            print "File tmpcheck.txt exists. It means you could be running friendChunkCheck.sh. If not, remove it, and run it again."
            sys.exit()
        else:
            os.system(pyCheckScript)
        print "Done. The list of files is in tmpcheck.txt."

    condor_file_name = options.subfile
    if options.queue == "condor":
        condor_file = open(condor_file_name,'w')
        condor_file.write('''Universe = vanilla
Executable = lxbatch_runner.sh
Log        = {jd}/logs/resubjob$(ProcId).log
Output     = {jd}/logs/resubjob$(ProcId).out
Error      = {jd}/logs/resubjob$(ProcId).error

use_x509userproxy = True
getenv = True
request_memory = 2000
+MaxRuntime = {rt}\n
'''.format(jd=fdir, rt=options.maxruntime * 60, here=os.environ['PWD'] ) )
        
    tmpfile = 'tmpcheck.txt' if len(args)<3 else args[2]
    txtfile=open(tmpfile,'r')
    for line in txtfile:
        l = line.rstrip()
        if l.startswith('#') or l.startswith('DONE') or l.startswith('check file'): continue

        if not l.endswith('OK'):
            base = os.path.basename(l)
            tokens = base.split('.')
            if len(tokens)<2: continue
            dataset = '_'.join(tokens[0].split('_')[:-1])
            chunk = tokens[1].split('chunk')[-1]
            #print "# resubmitting dataset = ",dataset," chunk ",chunk
            
            cmd = "python prepareEventVariablesFriendTree.py -j 0 -t NanoAOD --compression ZLIB:3 -I CMGTools.TTHAnalysis.tools.nanoAOD.ttH_modules {modules} {maintreedir} {frienddir} -N {nevents} -d {dataset} -c {chunk}".format(maintreedir=maindir,frienddir=fdir,nevents=options.chunkSize,dataset=dataset,chunk=chunk,modules=','.join(MODULES))
            if options.queue == 'condor':
                condor_file.write('Arguments = {pwd} {cmssw} {cmd} \nqueue 1 \n\n'.format(pwd=os.getcwd(),cmssw = os.environ['CMSSW_BASE'],cmd=cmd))
            else:
                print cmd

    if options.queue == 'condor':
        condor_file.close()
        print "condor submission file is ", condor_file_name
