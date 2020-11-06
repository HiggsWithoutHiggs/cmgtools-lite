from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module 
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection as NanoAODCollection 

from CMGTools.TTHAnalysis.treeReAnalyzer import Collection as CMGCollection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import declareOutput, writeOutput
from CMGTools.TTHAnalysis.tools.tfTool import TFTool
import os 
from copy import deepcopy

class finalDNN_3l(Module):
    def __init__(self, variations=[], doSystJEC=False, fillInputs=False):
        self.outVars = []
        self._MVAs   = [] 
        self.fillInputs = fillInputs
        varorder = ["lep1_pt","lep1_eta","lep1_phi","lep1_charge","lep2_pt","lep2_eta","lep2_phi","maxeta","Dilep_pdgId","jet1_pt","jet1_eta","jet1_phi","jet2_pt","jet2_eta","jet2_phi","jet3_pt","jet3_eta","jet3_phi","jetFwd1_pt","jetFwd1_eta","n_presel_jetFwd","n_presel_jet","nTopjet","topjet1_pt","topjet1_eta","topjet1_phi","topjet1_Ttag","nBJetLoose","nBJetMedium","lep3_pt","lep3_eta","lep3_phi"]
        cats_3l = ['predictions_tjvv','predictions_tt','predictions_other']

        if fillInputs:
            self.outVars.extend(varorder+['nEvent'])
        self.inputHelper = self.getVarsForVariation('')
        self.inputHelper['nEvent'] = lambda ev : ev.event

        self.systsJEC = {0:"",\
                         1:"_jesTotalCorrUp"  , -1:"_jesTotalCorrDown",\
                         2:"_jesTotalUnCorrUp", -2: "_jesTotalUnCorrDown",\
                         3:"_jerUp", -3: "_jerDown",\
                     } if doSystJEC else {0:""}
        if len(variations): 
            self.systsJEC = {0:""}
            for i,var in enumerate(variations):
                self.systsJEC[i+1]   ="_%sUp"%var
                self.systsJEC[-(i+1)]="_%sDown"%var

        for var in self.systsJEC: 
            self._MVAs.append( TFTool('DNN_3l%s'%self.systsJEC[var], os.environ['CMSSW_BASE'] + '/src/CMGTools/TTHAnalysis/data/kinMVA/hwh/model_3l_3classes_tightpresel.pb',
                               self.getVarsForVariation(self.systsJEC[var]), cats_3l, varorder))

            self.outVars.extend( ['DNN_3l%s_'%self.systsJEC[var] + x for x in cats_3l])


        vars_3l_unclUp = deepcopy(self.getVarsForVariation(''))
        vars_3l_unclUp["metLD"            ] =  lambda ev : (ev.MET_pt_unclustEnUp if ev.year != 2017 else ev.METFixEE2017_pt_unclustEnUp) *0.6 + ev.mhtJet25_Recl*0.4
        vars_3l_unclUp["mT_lep1"          ] =  lambda ev : ev.MT_met_lep1_unclustEnUp
        vars_3l_unclUp["mT_lep2"          ] =  lambda ev : ev.MT_met_lep2_unclustEnUp
        self.outVars.extend( ['DNN_3l_unclUp_' + x for x in cats_3l])

        vars_3l_unclDown = deepcopy(self.getVarsForVariation(''))
        vars_3l_unclDown["metLD"            ] =  lambda ev : (ev.MET_pt_unclustEnDown if ev.year != 2017 else ev.METFixEE2017_pt_unclustEnDown) *0.6 + ev.mhtJet25_Recl*0.4
        vars_3l_unclDown["mT_lep1"          ] =  lambda ev : ev.MT_met_lep1_unclustEnDown
        vars_3l_unclDown["mT_lep2"          ] =  lambda ev : ev.MT_met_lep2_unclustEnDown
        self.outVars.extend( ['DNN_3l_unclDown_' + x for x in cats_3l])

        worker_3l_unclUp        = TFTool('DNN_3l_unclUp', os.environ['CMSSW_BASE'] + '/src/CMGTools/TTHAnalysis/data/kinMVA/hwh/model_3l_3classes_tightpresel.pb',
                                           vars_3l_unclUp, cats_3l, varorder)
        worker_3l_unclDown      = TFTool('DNN_3l_unclDown', os.environ['CMSSW_BASE'] + '/src/CMGTools/TTHAnalysis/data/kinMVA/hwh/model_3l_3classes_tightpresel.pb',
                                           vars_3l_unclDown, cats_3l, varorder)
        
        self._MVAs.extend( [worker_3l_unclUp, worker_3l_unclDown])



    def getVarsForVariation(self, var ): 
        return {     "lep1_pt"          : lambda ev : ev.LepClean_Recl_pt[0] if ev.nLepFO_Recl >= 1 else -9,
                     "lep1_eta"         : lambda ev : ev.LepClean_Recl_eta[0] if ev.nLepFO_Recl >= 1 else 0,
                     "lep1_phi"         : lambda ev : ev.LepClean_Recl_phi[0] if ev.nLepFO_Recl >= 1 else -9,
                     "lep1_charge"      : lambda ev : ev.LepClean_Recl_charge[0] if ev.nLepFO_Recl >= 1 else -9,
                     
                     "lep2_pt"          : lambda ev : ev.LepClean_Recl_pt[1] if ev.nLepFO_Recl >= 2 else -9,
                     "lep2_eta"         : lambda ev : ev.LepClean_Recl_eta[1] if ev.nLepFO_Recl >= 2 else -9,
                     "lep2_phi"         : lambda ev : ev.LepClean_Recl_phi[1] if ev.nLepFO_Recl >= 2 else -9,

                     "lep3_pt"          : lambda ev : ev.LepClean_Recl_pt[2] if ev.nLepFO_Recl >= 3 else -9,
                     "lep3_eta"         : lambda ev : ev.LepClean_Recl_eta[2] if ev.nLepFO_Recl >= 3 else -9,
                     "lep3_phi"         : lambda ev : ev.LepClean_Recl_phi[2] if ev.nLepFO_Recl >= 3 else -9,
                     
                     "maxeta"           : lambda ev : max( [abs(ev.LepClean_Recl_eta[0]), abs(ev.LepClean_Recl_eta[1])]),
                     "Dilep_pdgId"      : lambda ev : (28 - abs(ev.LepClean_Recl_pdgId[0]) - abs(ev.LepClean_Recl_pdgId[1]))/2,
                     
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


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        print self.outVars
        declareOutput(self, wrappedOutputTree, self.outVars)
        
    def analyze(self,event):
        myvars = [event.iLepFO_Recl[0],event.iLepFO_Recl[1],event.iLepFO_Recl[2]]
        ret = []
        if self.fillInputs and self.inputHelper:
            for var in self.inputHelper:
                ret.append( (var, self.inputHelper[var](event)))
        for worker in self._MVAs:
            name = worker.name
            if ( not hasattr(event,"nJet25_jerUp_Recl") and not hasattr(event, "nJet25_jesBBEC1_yearDown_Recl")) and ('_jes' in name or  '_jer' in name or '_uncl' in name): continue # using jer bc components wont change
            ret.extend( [(x,y) for x,y in worker(event).iteritems()])
        writeOutput(self, dict(ret))
        return True
