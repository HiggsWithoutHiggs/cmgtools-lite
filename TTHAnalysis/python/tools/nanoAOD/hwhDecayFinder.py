from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 


class hwhDecayFinder( Module ): 
    def __init__(self, lepCollection="LepGood", genCollection="GenPart"):
        self.lepCollection = lepCollection
        self.genCollection = genCollection

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch('GenTopDecayMode','I')
        self.out.branch('GenV1DecayMode','I')
        self.out.branch('GenV2DecayMode','I')
        self.out.branch('nLepGood','I')
        self.out.branch('LepGood_genAncestor','I', lenVar='nLepGood')

    def getTopDecay(self, idx, genParts):
        decayproducts = [] 
        for g in genParts:
            if g.genPartIdxMother == idx:
                if abs(g.pdgId) == 6 or abs(g.pdgId) == 24:
                    return self.getTopDecay(genParts.index(g),genParts)
                else:
                    decayproducts.append( abs(g.pdgId) )
        if len(decayproducts) > 2: print "More than two top decay products?", decayproducts
        if not len(decayproducts): print 'No top decay products found...'
        return decayproducts

    def getPromptVDecay(self, idx, genParts):
        decayproducts = [] 
        for g in genParts:
            if g.genPartIdxMother == idx: 
                if (abs(g.pdgId) in [23,24]):
                    return self.getPromptVDecay( genParts.index(g),genParts)
                else: 
                    decayproducts.append( abs(g.pdgId) )
        if len(decayproducts) > 2: print "More than two V decay products?", decayproducts
        if not len(decayproducts): print 'No V decay products found...'
        return decayproducts

    def getLeptonAncestor(self, idx, genParts):
        ancestor = -1
        if idx>0: 
            for g in genParts:
                if genParts.index(g) == idx:
                    if g.genPartIdxMother > 0:
                        return self.getLeptonAncestor( g.genPartIdxMother, genParts)
                    elif g.genPartIdxMother==0 and (g.statusFlags>>7)&1==1:
                        ancestor = g.pdgId
                    else:
                        break
            #if ancestor == -1: print "No lepton ancestor found?"
        return ancestor

    def analyze(self, event):
        leps = [ l for l in Collection(event, self.lepCollection) ]
        genParts = [ g for g in Collection(event, self.genCollection) ]
        promptVDecays = []
        topDecay = None
        ancestors = []

        # look for the TVV decay chain
        for g in genParts:
            if abs(g.pdgId) in [23,24] and g.genPartIdxMother==0:
                promptVDecay = self.getPromptVDecay(genParts.index(g), genParts)
                promptVDecays.append(promptVDecay)
            elif abs(g.pdgId) == 6 and g.genPartIdxMother==0:
                topDecay = self.getTopDecay(genParts.index(g), genParts)
        
        if not topDecay:  
            self.out.fillBranch('GenTopDecayMode',-1);
        elif   topDecay == [11,12] or topDecay == [12,11] : 
            self.out.fillBranch('GenTopDecayMode',11);
        elif   topDecay == [13,14] or topDecay == [14,13] : 
            self.out.fillBranch('GenTopDecayMode',13);
        elif   topDecay == [15,16] or topDecay == [16,15] : 
            self.out.fillBranch('GenTopDecayMode',15);
        else:  
            self.out.fillBranch('GenTopDecayMode',1);

        for v,vDecay in enumerate(promptVDecays):
            if   vDecay == [11,12] or vDecay == [12,11] :
                self.out.fillBranch('GenV%dDecayMode'%(v+1),11)
            elif vDecay == [13,14] or vDecay == [14,13] : 
                self.out.fillBranch('GenV%dDecayMode'%(v+1),13);
            elif vDecay == [15,16] or vDecay == [16,15] : 
                self.out.fillBranch('GenV%dDecayMode'%(v+1),15);
            else:
                self.out.fillBranch('GenV%dDecayMode'%(v+1),1);

        # look for the ancestors of the matched leptons
        for l in leps:
            idxGen = l.genPartIdx
            ancestor = self.getLeptonAncestor(idxGen, genParts)
            ancestors.append(ancestor)
        self.out.fillBranch("LepGood_genAncestor", ancestors)

        return True

        
