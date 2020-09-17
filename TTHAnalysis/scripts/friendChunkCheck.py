#!/usr/bin/env python

import sys,os,re
import ROOT

if __name__ == "__main__":
    
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] dir_with_friend_trees dir_with_trees [checkedfile.txt]")
    parser.add_option("-z", "--zombies",  dest="zombies", action='store_true', default=False, help="Check zombies");
    parser.add_option("-D", "--dm", "--dataset-match", dest="datasetMatches",  type="string", default=[], action="append", help="Process only this dataset (or dataset if specified multiple times): REGEXP");
    parser.add_option("--xD", "--de", "--dataset-exclude", dest="datasetExcludes",  type="string", default=[], action="append", help="Exclude these dataset (or dataset if specified multiple times): REGEXP");
    (options, args) = parser.parse_args()

    fdir = args[0]
    allfiles = os.listdir(fdir)

    chunks = {}
    samples = []
    for f in allfiles:
        if f=='logs': continue
        tokens = (os.path.basename(f)).split('.')
        if len(tokens)<2: continue
        s = tokens[0]
        c = tokens[1]
        if not re.match('chunk\d+',c): continue
        if options.datasetMatches != []:
            found = False
            for dm in  options.datasetMatches:
                if re.match(dm,s): found = True
            if not found: continue
        if options.datasetExcludes != []:
            found = False
            for dm in  options.datasetExcludes:
                if re.match(dm,s): found = True
            if found: continue
        print "check file ",f
        if s not in samples: samples.append(s)
        if s not in chunks:
            chunks[s] = [c]
        else:
            chunks[s].append(c)

    goodfiles = {}
    for s,c in chunks.iteritems():
        c.sort(key=lambda x: int(x.split('chunk')[-1]))
        maxc = int(c[-1].split('chunk')[-1])
        #print "Sample ",s," has ",len(c)," chunks: ",c,". Max chunk = ",maxc
        for ic in xrange(maxc+1):
            ftest = '{fdir}/{sample}.chunk{c}.root'.format(fdir=fdir,sample=s,c=ic)
            if not os.path.isfile(ftest):
                print ftest, " # not present"
            else:
                if os.stat(ftest)<1000:
                    print ftest, " # has zero size"
                else:
                    if s not in goodfiles:
                        goodfiles[s] = [ftest]
                    else:
                        goodfiles[s].append(ftest)
    
    if options.zombies:
        print "# Testing for zombies or not correctly closed files"
        for s,chunks in goodfiles.iteritems():
            good = True
            for f in chunks:
                tf = ROOT.TFile(f)
                if not tf: good=False
                else:
                    if tf.IsZombie() or tf.TestBit(ROOT.TFile.kRecovered): good=False
                print '{f} # {qual}'.format(f=f,qual="OK" if good else "zombie")
    
    print "DONE."