from __future__ import print_function

import ROOT as r 
import numpy as np
import pickle,math,os,tqdm
from keras.utils import np_utils
from multiprocessing import Pool

testCut   = lambda ev: ev.event%5==0
trainCut  = lambda ev: ev.event%5!=0

def getVar(ev, l, s):
    return getattr(ev,'LepClean_Recl_%s'%s)[l]


commonFeatureList = [
    "lep1_pt"      ,
    "lep1_eta"         ,
    "lep1_phi"         ,
    "lep1_charge"      ,

    "lep2_pt"      ,
    "lep2_eta"         ,
    "lep2_phi"         ,

    "maxeta"           ,
    "Dilep_pdgId"      , 

    "jet1_pt"          ,
    "jet1_eta"         ,
    "jet1_phi"         ,    

    "jet2_pt"          ,
    "jet2_eta"         ,
    "jet2_phi"         ,

    "jet3_pt"          ,
    "jet3_eta"         ,
    "jet3_phi"         ,
    
    "jetFwd1_pt"       ,
    "jetFwd1_eta"      ,

    "n_presel_jetFwd"  ,
    "n_presel_jet"     ,

    "nTopjet"          ,
    "topjet1_pt"       ,
    "topjet1_eta"      ,
    "topjet1_phi"      ,
    "topjet1_Ttag"     ,
    
    "nBJetLoose"       ,
    "nBJetMedium"      ,
]


features = {

    "lep1_pt"      : lambda ev : ev.LepClean_Recl_pt[0] if getattr(ev,'nLepClean_Recl') >= 1 else -9,
    "lep1_eta"         : lambda ev : ev.LepClean_Recl_eta[0] if getattr(ev,'nLepClean_Recl') >= 1 else 0,
    "lep1_phi"         : lambda ev : ev.LepClean_Recl_phi[0] if getattr(ev,'nLepClean_Recl') >= 1 else -9,
    "lep1_charge"      : lambda ev : ev.LepClean_Recl_charge[0] if getattr(ev,'nLepClean_Recl') >= 1 else -9,

    "lep2_pt"      : lambda ev : ev.LepClean_Recl_pt[1] if getattr(ev,'nLepClean_Recl') >= 2 else -9,
    "lep2_eta"         : lambda ev : ev.LepClean_Recl_eta[1] if getattr(ev,'nLepClean_Recl') >= 2 else -9,
    "lep2_phi"         : lambda ev : ev.LepClean_Recl_phi[1] if getattr(ev,'nLepClean_Recl') >= 2 else -9,

    "lep3_pt"      : lambda ev : ev.LepClean_Recl_pt[2] if getattr(ev,'nLepClean_Recl') >= 3 else -9,
    "lep3_eta"         : lambda ev : ev.LepClean_Recl_eta[2] if getattr(ev,'nLepClean_Recl') >= 3 else -9,
    "lep3_phi"         : lambda ev : ev.LepClean_Recl_phi[2] if getattr(ev,'nLepClean_Recl') >= 3 else -9,

    "maxeta"           : lambda ev : max( [abs(ev.LepClean_Recl_eta[0]), abs(ev.LepClean_Recl_eta[1])]) if getattr(ev,'nLepClean_Recl') >=2 else -9,
    "Dilep_pdgId"      : lambda ev : (28 - abs(ev.LepClean_Recl_pdgId[0]) - abs(ev.LepClean_Recl_pdgId[1]))/2 if getattr(ev,'nLepClean_Recl') >=2 else -9,
    
    "jet1_pt"          : lambda ev : getattr(ev,'JetSel_Recl_pt')[0] if getattr(ev,'nJet25_Recl') > 0 else -9,
    "jet1_eta"         : lambda ev : abs(ev.JetSel_Recl_eta[0]) if getattr(ev,'nJet25_Recl') > 0 else 9,
    "jet1_phi"         : lambda ev : ev.JetSel_Recl_phi[0] if getattr(ev,'nJet25_Recl') > 0 else -9,
    
    "jet2_pt"          : lambda ev : getattr(ev,'JetSel_Recl_pt')[1] if getattr(ev,'nJet25_Recl') > 1 else -9,
    "jet2_eta"         : lambda ev : abs(ev.JetSel_Recl_eta[1]) if getattr(ev,'nJet25_Recl') > 1 else 9,
    "jet2_phi"         : lambda ev : ev.JetSel_Recl_phi[1] if getattr(ev,'nJet25_Recl') >= 2 else -9,

    "jet3_pt"          : lambda ev : ev.JetSel_Recl_pt[2] if getattr(ev,'nJet25_Recl') >= 3 else -9,
    "jet3_eta"         : lambda ev : ev.JetSel_Recl_eta[2] if getattr(ev,'nJet25_Recl') >= 3 else -9,
    "jet3_phi"         : lambda ev : ev.JetSel_Recl_phi[2] if getattr(ev,'nJet25_Recl') >= 3 else -9,
    
    "jetFwd1_pt"       : lambda ev : getattr(ev,'FwdJet1_pt_Recl') if getattr(ev,'nFwdJet_Recl') else -9,    
    "jetFwd1_eta"      : lambda ev : abs(getattr(ev,'FwdJet1_eta_Recl')) if getattr(ev,'nFwdJet_Recl') else 9,

    "n_presel_jetFwd"  : lambda ev : getattr(ev,'nFwdJet_Recl'),
    "n_presel_jet"     : lambda ev : getattr(ev,'nJet25_Recl'),

    "nTopjet"          : lambda ev : getattr(ev,'nTJetSel_Recl'),
    "topjet1_pt"       : lambda ev : getattr(ev,'TJetSel_Recl_pt')[0] if getattr(ev,'nTJetSel_Recl') > 0 else -9,
    "topjet1_eta"      : lambda ev : getattr(ev,'TJetSel_Recl_eta')[0] if getattr(ev,'nTJetSel_Recl') > 0 else -9,
    "topjet1_phi"      : lambda ev : getattr(ev,'TJetSel_Recl_phi')[0] if getattr(ev,'nTJetSel_Recl') > 0 else -9,
    "topjet1_Ttag"     : lambda ev : getattr(ev,'TJetSel_Recl_deepTag_TvsQCD')[0] if getattr(ev,'nTJetSel_Recl') > 0 else -9,

    "nBJetLoose"       : lambda ev : getattr(ev,'nBJetLoose25_Recl'),
    "nBJetMedium"      : lambda ev : getattr(ev,'nBJetMedium25_Recl'),

    }

cuts = {
    '2lss'       : {'sigll' : lambda ev: ev.GenV1DecayMode>1 and ev.GenV2DecayMode>1 and ev.nLepClean_Recl>=2 and ev.LepClean_Recl_pt[0]>50 and ev.LepClean_Recl_pt[1]>30 and ev.nBJetMedium25_Recl>=1, 
                    'ttv'   : lambda ev: ev.nLepClean_Recl>=2 and ev.LepClean_Recl_pt[0]>50 and ev.LepClean_Recl_pt[1]>30 and ev.nBJetMedium25_Recl>=1,
                    'tt'    : lambda ev: ev.nLepClean_Recl>=2 and ev.LepClean_Recl_pt[0]>50 and ev.LepClean_Recl_pt[1]>30 and ev.nBJetMedium25_Recl>=1,
                    'other' : lambda ev: ev.nLepClean_Recl>=2 and ev.LepClean_Recl_pt[0]>50 and ev.LepClean_Recl_pt[1]>30 and ev.nBJetMedium25_Recl>=1},
    '3l'         : {'sigll' : lambda ev: ev.GenV1DecayMode>1 and ev.GenV2DecayMode>1 and ev.nLepClean_Recl >= 3, 'ttv' : lambda ev: ev.nLepClean_Recl >= 3, 'tt' : lambda ev: ev.nLepClean_Recl >= 3, 'other' : lambda ev: ev.nLepClean_Recl >= 3}
}

classes = {
    'sigll'      : { 'cut': cuts['2lss']['sigll'], 'lst_train' : [], 'lst_test' : [] , 'lst_y_train' : [], 'lst_y_test' : [] },
    'ttv'        : { 'cut': cuts['2lss']['ttv'],   'lst_train' : [], 'lst_test' : [] , 'lst_y_train' : [], 'lst_y_test' : [] },
    'tt'         : { 'cut': cuts['2lss']['tt'],    'lst_train' : [], 'lst_test' : [] , 'lst_y_train' : [], 'lst_y_test' : [] },
    'other'      : { 'cut': cuts['2lss']['other'], 'lst_train' : [], 'lst_test' : [] , 'lst_y_train' : [], 'lst_y_test' : [] },
}

sampleDir='/eos/cms/store/cmst3/group/wmass/secret/NanoTrees_HWH_2lskim_170920/2018/'

sigSamplesWpWp = []
sigSamplesWZ = []
ttvSamples = []
ttSamples  = []
othSamples = []

sigSamplesWpWp.extend( ['TJWpWp_SM_2018.root','TJWpWp_0p8_2018.root','TJWpWm_SM_2018.root','TJWpWm_0p8_2018.root'] )
sigSamplesWZ.extend( ['TJWZ_SM_2018.root'] )
ttvSamples.extend( ['TTWToLNu_fxfx.root']+['TTZToLLNuNu_amc_part%d.root'%i for i in range(1,2)] )
ttSamples.extend( ['TTJets_SingleLeptonFromT.root','TTJets_SingleLeptonFromTbar.root']+['TTJets_DiLepton_part%d.root'%i for i in range(1,2)])
othSamples.extend( ['WWTo2L2Nu.root','WZTo3LNu_fxfx.root','ZZTo4L.root'] )


def toNumpy(featureList,maxEntries,task):
    print('starting', task)
    fil, typs = task
    path = os.path.dirname(fil)
    fname = os.path.basename(fil)
    friendpath = path+'/3_recleaner_vMarcAllVars/'
    print('List of features for', featureList + eval('featureList'))
    tfile = r.TFile(fil); ttree = tfile.Events; ttree.AddFriend('Friends',friendpath+fname.replace('.root','_Friend.root'))
    results = {}
    for ty in typs: 
        results[ty + '_test']  = []
        results[ty + '_train'] = []

    print("start looping on file ",fil)
    for iev,ev in enumerate(ttree):
        if iev%1000==0: print('Processing event ',iev,'...')
        if iev>maxEntries: break
        tstr = 'test' if testCut(ev) else 'train'
        for ty in typs:
            if classes[ty]['cut'](ev):
                results[ty+'_'+tstr].append([ features[s](ev) for s in (featureList) ])
    tfile.Close()
    print('finishing', task)
    return results




if __name__ == "__main__":

    from functools import partial

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options]")
    parser.add_option("--max-entries",   dest="maxEntries", default=1000000000, type="int", help="Max entries to process in each tree")
    parser.add_option("-o", "--outfile", dest="outfile", type="string", default="vars.pkl", help="Output pickle file (default: vars.pkl)");
    parser.add_option("-c", "--channel", dest="channel", type="string", default="2lss", help="Final state: 2lss or 3l (default: 2lss)");
    (options, args) = parser.parse_args()

    tasks = []
    print('Setting up the tasks')
    sigSamples = sigSamplesWpWp if options.channel=='2lss' else sigSamplesWZ
    for samp in sigSamples:
        tasks.append( (sampleDir+'/'+samp, ['sigll']) )

    for samp in ttvSamples:
        tasks.append( (sampleDir+'/'+samp, ['ttv']) )

    for samp in ttSamples:
        tasks.append( (sampleDir+'/'+samp, ['tt']) )
    
    for samp in othSamples:
        tasks.append( (sampleDir+'/'+samp, ['other']) )
    
    print('Numpy conversion. It will take time...')
    print("max entries = ",options.maxEntries)

    featureList = commonFeatureList
    if options.channel=='3l':
        featureList += ["lep3_pt","lep3_eta","lep3_phi"]
        for cl,vals in classes.iteritems():
            vals['cut'] = cuts[options.channel][cl]
            
    ## lxplus seems to have 10 cores/each
    p =  Pool(min(10,len(sigSamples+ttvSamples+ttSamples+othSamples)))
    func = partial(toNumpy,featureList,options.maxEntries)
    results = list(tqdm.tqdm(p.imap(func, tasks), total=len(tasks)))

    print('Now putting everything together')

    types = ['sigll', 'ttv', 'tt', 'other']
    for result in results: 
        for ty in types:
            if ty+'_train' in result:
                classes[ty]['lst_train'].extend( result[ty+'_train'])
                classes[ty]['lst_test' ].extend( result[ty+'_test'])

            
    print('Setting the indices')
    toDump = {} 
    for i, ty in enumerate(types):
        classes[ty]['lst_train'  ] = np.asarray(classes[ty]['lst_train'])
        classes[ty]['lst_y_train'] = i*np.ones((classes[ty]['lst_train'].shape[0],1))
        classes[ty]['lst_test'   ] = np.asarray(classes[ty]['lst_test'])
        classes[ty]['lst_y_test' ] = i*np.ones((classes[ty]['lst_test'].shape[0],1))

    train_x = np.concatenate( tuple( [classes[ty]['lst_train'] for ty in types] ), axis=0)
    train_y = np_utils.to_categorical( np.concatenate( tuple( [classes[ty]['lst_y_train'] for ty in types] ), axis=0), len(classes))
    test_x = np.concatenate( tuple( [classes[ty]['lst_test'] for ty in types] ), axis=0)
    test_y = np_utils.to_categorical( np.concatenate( tuple( [classes[ty]['lst_y_test'] for ty in types] ), axis=0), len(classes))
    toDump['train_x'] = train_x
    toDump['train_y'] = train_y
    toDump['test_x' ] = test_x
    toDump['test_y' ] = test_y

    ### dump to file
    print ('dump to ',options.outfile,' now...')
    pickle_out = open(options.outfile,'wb')
    pickle.dump( toDump, pickle_out)
    pickle_out.close()
