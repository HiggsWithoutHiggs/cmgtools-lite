#!/usr/bin/env python
import sys,os.path,re
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

class CheckOneFriend:
    def __init__(self,inputdir,dataset,frienddir="friends",friendTree="Friends",
                 skipFriend=False):
        self.file = "/".join([inputdir,dataset])
        self.friendfile = "/".join([frienddir,dataset]).replace('.root','_Friend.root')
        self.tree = 'Events'
        self.ftree = friendTree
        self.skipFriend = skipFriend
    def check(self,verbose=0):
        if not os.path.isfile(self.file): 
            print "File %s does not exist!" % self.file
            return False
        tf = ROOT.TFile.Open(self.file)
        if not tf:
            print "ERROR! Unable to open file %s" % (self.file)
            return False
        if tf.TestBit(ROOT.TFile.kRecovered):
            print "WARNING! Attempt was made to recover file %s" % (self.file)
            return False
        tree = tf.Get(self.tree)
        if not tree:
            print "ERROR! Base tree in file %s is missing" % (self.file)
            return False
        tentries = tree.GetEntries()
        tf.Close()
        if self.skipFriend:
            if verbose>0:
                print "%s OK." % self.file
            return True
        else:
            if not os.path.isfile(self.friendfile): 
                print "Friend file %s does not exist!" % self.friendfile
                return False
            ff = ROOT.TFile.Open(self.friendfile)
            if not ff:
                print "ERROR! Unable to open file %s" % (self.friendfile)
                return False
            if ff.TestBit(ROOT.TFile.kRecovered):
                print "WARNING! Attempt was made to recover file %s" % (self.friendfile)
                return False
            ftree = ff.Get(self.ftree)
            if not ftree:
                print "ERROR! Friend tree in file %s is missing" % (self.friendfile)
                return False
            fentries = ftree.GetEntries()
            ff.Close()
            if tentries != fentries:
                print "ERROR! Friend tree for dataset %s has %d entries, while tree has %d entries!" % (self.file,fentries,tentries)
                return False
            else:
                if verbose>0:
                    print "%s OK." % self.file
                return True


MODULES_DATA = ["recleaner_step1","recleaner_step2_data","triggerSequence"]
MODULES_MC = ["recleaner_step1","recleaner_step2_mc","mcMatch_seq","higgsDecay","triggerSequence"]

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] dir_with_trees dir_with_friends")
    parser.add_option(      "--sp", dest="selectProcess",  type="string", default=None, help="Process only datasets that match this regexp (or comma-separated list of regexp)");
    parser.add_option(      "--xp", dest="excludeProcess",  type="string", default=None, help="Process only datasets that do not match this regexp (or comma-separated list of regexp)");
    parser.add_option("--skipFriend", dest="skipFriend", action="store_true", default=False, help="Check only main trees (useful if for some reason you don't have friends yet)");
    parser.add_option(      "--data",  dest="data", action="store_true", default=False, help="Run the data modules (default is MC modules)");
    parser.add_option("-N", "--events",  dest="chunkSize", type="int",    default=1000000, help="Default chunk size when splitting trees");    
    (options, args) = parser.parse_args()

    if len(args)<2:
        print "Usage: python checkMergedFriends.py <dir_with_trees> <dir_with_friends> [options].\n"
        exit(1)

    MODULES = MODULES_DATA if options.data else MODULES_MC

    sel_processes = []
    if options.selectProcess != None:
        sel_processes = options.selectProcess.split(',')
    excl_processes = []
    if options.excludeProcess != None:
        excl_processes = options.excludeProcess.split(',')

    inputdir = args[0]
    frienddir = args[1]
    badDir = []
    for root,dirs,files in os.walk(inputdir):
        for f in files:
            if f.endswith(".txt"): continue  # skip txt files with selection and saved branch (when skimming ntuples)
            if options.selectProcess and not any(re.match(proc,f) for proc in sel_processes): continue
            if options.excludeProcess and any(re.match(proc,f) for proc in excl_processes): continue
            print "Checking dataset %s..." % f
            cf = CheckOneFriend(inputdir,f, frienddir=frienddir,skipFriend=options.skipFriend)
            if not cf.check(1):
                badDir.append(f)
        break
    print "-"*20
    print "Summary"
    print "-"*20
    if len(badDir):
        print "There are %d bad freiends. They are:" % len(badDir)
        for i in badDir:
            print i
        print "This is the command to resub."
        datasetMatch = ' -D '.join([ds.replace('.root','') for ds in badDir])
        cmd = "python prepareEventVariablesFriendTree.py -j 0 -t NanoAOD --compression ZLIB:3 -I CMGTools.TTHAnalysis.tools.nanoAOD.ttH_modules {modules} {maintreedir} {frienddir} -N {nevents} -q condor -D {dm} --sub resub.sub".format(maintreedir=inputdir,frienddir=frienddir,nevents=options.chunkSize,dm=datasetMatch,modules=','.join(MODULES))
        print cmd
    else:
        print "All files are OK"

