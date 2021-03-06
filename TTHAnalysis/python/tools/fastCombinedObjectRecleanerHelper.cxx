#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <DataFormats/Math/interface/deltaR.h>
#include <CMGTools/TTHAnalysis/interface/CollectionSkimmer.h>
#include "CMGTools/TTHAnalysis/interface/CombinedObjectTags.h"
#include "DataFormats/Math/interface/LorentzVector.h"

struct JetSumCalculatorOutput {
  int thr;
  float htJetj;
  float mhtJet;
  int nBJetLoose;
  int nBJetMedium; 
  int nJet;
  int nFwdJet;
  float fwd1_pt;
  float fwd1_eta;
};

class fastCombinedObjectRecleanerHelper {
public:
  typedef TTreeReaderValue<unsigned>   ruint;
  typedef TTreeReaderValue<int>   rint;
  typedef TTreeReaderArray<float> rfloats;
  typedef TTreeReaderArray<int> rints;
  class rcount {
      public:
          rcount() : signed_(NULL), unsigned_(NULL) {}
          rcount(rint *src) : signed_(src), unsigned_(NULL) {}
          rcount(ruint *src) : signed_(NULL), unsigned_(src) {}
          rcount & operator=(rint *src) { signed_ = src; return *this; }  
          rcount & operator=(ruint *src) { unsigned_ = src; return *this; }  
          int operator*() const { return signed_ ? **signed_ : int(**unsigned_); }
      private:
          rint * signed_;
          ruint * unsigned_;
  };
  
  fastCombinedObjectRecleanerHelper(CollectionSkimmer &clean_leps, CollectionSkimmer &clean_taus, CollectionSkimmer &clean_jets, CollectionSkimmer &clean_tjets, CollectionSkimmer &clean_vjets, CollectionSkimmer &clean_fjets, CollectionSkimmer &clean_ljets, bool cleanJetsWithFOTaus, float bTagL, float bTagM, bool cleanWithRef=false) : clean_leps_(clean_leps), clean_taus_(clean_taus), clean_jets_(clean_jets), clean_tjets_(clean_tjets), clean_vjets_(clean_vjets), clean_fjets_(clean_fjets), clean_ljets_(clean_ljets), deltaR2cut(0.16), cleanJetsWithFOTaus_(cleanJetsWithFOTaus), bTagL_(bTagL), bTagM_(bTagM), cleanWithRef_(cleanWithRef), deltaR2cut_taus(0.09), deltaR2cut_fatjets(0.64) {
    _cl.reset(new std::vector<int>);
    _ct.reset(new std::vector<int>);
    _cj.reset(new std::vector<int>);
    _ctj.reset(new std::vector<int>);
    _cvj.reset(new std::vector<int>);
    _cfj.reset(new std::vector<int>);
    _clj.reset(new std::vector<int>);
}
  
  void setLeptons(rint *nLep, rfloats* lepPt, rfloats *lepEta, rfloats *lepPhi) {
    nLep_ = nLep; Lep_pt_ = lepPt; Lep_eta_ = lepEta; Lep_phi_ = lepPhi;
    if (!nLep || !lepPt || !lepEta || !lepPhi) { std::cout << "ERROR: fastCombinedObjectRecleanerHelper initialized setLeptons with a null reader" << std::endl; }
  }
  void setLeptons(ruint *nLep, rfloats* lepPt, rfloats *lepEta, rfloats *lepPhi, rints *lepJet) {
    nLep_ = nLep; Lep_pt_ = lepPt; Lep_eta_ = lepEta; Lep_phi_ = lepPhi; Lep_jet_ = lepJet;
    if (!nLep || !lepPt || !lepEta || !lepPhi || !lepJet) { std::cout << "ERROR: fastCombinedObjectRecleanerHelper initialized setLeptons with a null reader" << std::endl; }
  }
  void setTaus(rint *nTau, rfloats *tauPt, rfloats *tauEta, rfloats *tauPhi) {
    nTau_ = nTau; Tau_pt_ = tauPt; Tau_eta_ = tauEta; Tau_phi_ = tauPhi;
    if (!nTau || !tauPt || !tauEta || !tauPhi) { std::cout << "ERROR: fastCombinedObjectRecleanerHelper initialized setTaus with a null reader" << std::endl; }
  }
  void setTaus(ruint *nTau, rfloats *tauPt, rfloats *tauEta, rfloats *tauPhi, rints *tauJet) {
    nTau_ = nTau; Tau_pt_ = tauPt; Tau_eta_ = tauEta; Tau_phi_ = tauPhi; Tau_jet_ = tauJet;
    if (!nTau || !tauPt || !tauEta || !tauPhi || !tauJet) { std::cout << "ERROR: fastCombinedObjectRecleanerHelper initialized setTaus with a null reader" << std::endl; }
  }
  void setJets(ruint *nJet, rfloats *jetPt, rfloats *jetEta, rfloats *jetPhi, rfloats *jetbtagCSV, vector<rfloats*> jetpt) {
    nJet_ = nJet; Jet_pt_ = jetPt; Jet_eta_ = jetEta; Jet_phi_ = jetPhi; Jet_btagCSV_ = jetbtagCSV; 
    Jet_corr_   = jetpt;
  }
  void setFatJets(ruint *nJet, rfloats *jetPt, rfloats *jetEta, rfloats *jetPhi, rfloats *jetbtagCSV, rfloats *jetttag, rfloats *jetwtag, rfloats *jetztag) {
    nFJet_ = nJet; FJet_pt_ = jetPt; FJet_eta_ = jetEta; FJet_phi_ = jetPhi; FJet_btagCSV_ = jetbtagCSV; FJet_deepTag_TvsQCD_ = jetttag; FJet_deepTag_WvsQCD_ = jetwtag; FJet_deepTag_ZvsQCD_ = jetztag; 
  }

  void addJetPt(int pt){
    _jetptcuts.insert(pt);
  }
  void setFwdPt(float fwdJetPt1, float fwdJetPt2){
    fwdJetPt1_= fwdJetPt1;
    fwdJetPt2_= fwdJetPt2;
  }


  typedef math::PtEtaPhiMLorentzVectorD ptvec;
  typedef math::XYZTLorentzVectorD crvec;

  std::vector<JetSumCalculatorOutput> GetJetSums(int variation = 0){

    std::vector<JetSumCalculatorOutput> output;
    
    crvec _mht(0,0,0,0);
    
    for (int i=0; i<*nLep_; i++) {
      if (!sel_leps[i]) continue;
      crvec lep(ptvec((*Lep_pt_)[i],0,(*Lep_phi_)[i],0));
      _mht = _mht - lep;
    }
    if (cleanJetsWithFOTaus_) {
      for (auto i : *_ct){
	crvec tau(ptvec((*Tau_pt_)[i],0,(*Tau_phi_)[i],0));
	_mht = _mht - tau;
      }
    }
    
    for (auto thr : _jetptcuts){
      auto mht = _mht;
      JetSumCalculatorOutput sums;
      sums.thr = float(thr);
      sums.htJetj = 0;
      sums.nBJetLoose = 0;
      sums.nBJetMedium = 0;
      sums.nJet = 0;
      sums.nFwdJet = 0;
      sums.fwd1_pt = 0;
      sums.fwd1_eta = 0;
      int var = -1;
      if (variation < 0)
	var = 2*abs(variation) - 1;
      else
	var = 2*variation-2;


      for (auto j : *_cj){
	float pt = (*Jet_pt_)[j];
	if (variation != 0) 
	  pt = (*(Jet_corr_.at(var)))[j];
	float abseta = fabs((*Jet_eta_)[j]) ;
	if (abseta > 2.7 && abseta < 3 ){
	  if (pt  > fwdJetPt2_){
	    sums.nFwdJet++;
	    if (pt > sums.fwd1_pt){
	      sums.fwd1_pt = pt; sums.fwd1_eta = (*Jet_eta_)[j];
	    }
	  }
	  continue;
	}
	else if (abseta > 2.4 && abseta < 5){
	  if (pt > fwdJetPt1_){
	    sums.nFwdJet++;
	    if (pt > sums.fwd1_pt){
	      sums.fwd1_pt = pt; sums.fwd1_eta = (*Jet_eta_)[j];
	    }
	  }
	  continue;
	}
	else if(abseta > 2.4) continue;
	if (pt<=thr) continue;
	float phi = (*Jet_phi_)[j];
	float csv = (*Jet_btagCSV_)[j];
	sums.htJetj += pt;
	crvec jp4(ptvec(pt,0,phi,0));
	mht = mht - jp4;
	if (csv>bTagL_) sums.nBJetLoose += 1;
	if (csv>bTagM_) sums.nBJetMedium += 1;
	sums.nJet += 1;
      }

      sums.mhtJet = mht.Pt();
      output.push_back(sums);
    }

    return output;
  }
  
  void clear() {
    sel_leps.reset(new bool[*nLep_]);
    sel_leps_extrafortau.reset(new bool[*nLep_]);
    sel_taus.reset(new bool[*nTau_]);
    sel_jets.reset(new bool[*nJet_]);
    std::fill_n(sel_leps.get(),*nLep_,false);
    std::fill_n(sel_leps_extrafortau.get(),*nLep_,false);
    std::fill_n(sel_taus.get(),*nTau_,false);
    std::fill_n(sel_jets.get(),*nJet_,false);
  }
  void selectLepton(uint i, bool what=true) {sel_leps.get()[i]=what;}
  void selectLeptonExtraForTau(uint i, bool what=true) {sel_leps_extrafortau.get()[i]=what;}
  void selectTau(uint i, bool what=true) {sel_taus.get()[i]=what;}
  void selectJet(uint i, bool what=true) {sel_jets.get()[i]=what;}
  void loadTags(CombinedObjectTags *tags, bool cleanTausWithLooseLeptons, float wPL=0, float wPM=0){
    //std::cout << __LINE__ << " in loadtags" << std::endl;
    bTagL_ = wPL; bTagM_ = wPM;
    //std::copy(tags->lepsC.get(),tags->lepsC.get()+*nLep_,sel_leps.get()); //marc
    std::copy(tags->lepsT.get(),tags->lepsT.get()+*nLep_,sel_leps.get());
    if (cleanTausWithLooseLeptons) std::copy(tags->lepsL.get(),tags->lepsL.get()+*nLep_,sel_leps_extrafortau.get());
    std::copy(tags->tausF.get(),tags->tausF.get()+*nTau_,sel_taus.get());
    std::copy(tags->jetsS.get(),tags->jetsS.get()+*nJet_,sel_jets.get());
  }

  void setDR(float f) {deltaR2cut = f*f;}
  void setDRFatJets(float f) {deltaR2cut_fatjets = f*f;}

  std::vector<std::vector<int>* > run() {

//std::cout<< "=================="<<std::endl;
    clean_taus_.clear();
    clean_leps_.clear();
    clean_jets_.clear();
    clean_tjets_.clear();
    clean_vjets_.clear();
    clean_fjets_.clear();
    clean_ljets_.clear();

//std::cout<<__LINE__<<std::endl;
    _cl->clear();
    _ct->clear();
    _cj->clear();
    _ctj->clear();
    _cvj->clear();
    _cfj->clear();
    _clj->clear();

    std::vector<int> other_fatjets;
    std::vector<int> bad_leptons;
    bad_leptons.clear();

//std::cout<< "total number of leptons " << *nLep_ << std::endl;
//std::cout<< "total number of fat jets " << *nFJet_ << std::endl;
// beginning marc

    std::vector<float> vetos_eta;
    std::vector<float> vetos_phi;
    std::vector<int>   vetos_indices;
    for (int iL = 0, nL = *nLep_; iL < nL; ++iL) if (sel_leps[iL]) {vetos_eta.push_back((*Lep_eta_)[iL]); vetos_phi.push_back((*Lep_phi_)[iL]); if (Lep_jet_) vetos_indices.push_back((*Lep_jet_)[iL]);}
    for (int iT = 0, nT = *nTau_; iT < nT; ++iT) if (sel_taus[iT]) {vetos_eta.push_back((*Tau_eta_)[iT]); vetos_phi.push_back((*Tau_phi_)[iT]); if (Tau_jet_) vetos_indices.push_back((*Tau_jet_)[iT]);}

    std::unique_ptr<bool[]> good;   good.reset(new bool[*nJet_]);    std::fill_n(good.get(),*nJet_,true);
    std::unique_ptr<bool[]> goodfj; goodfj.reset(new bool[*nFJet_]); std::fill_n(goodfj.get(),*nFJet_,true);

    for (int iFJ = 0, nFJ = *nFJet_; iFJ < nFJ; ++iFJ) {
      float tagger_t = (*FJet_deepTag_TvsQCD_)[iFJ];
      if (tagger_t>0.436) { // these are the good top fatjets, WP iwth 2.5% mistag rate
        clean_tjets_.push_back(iFJ);
        _ctj->push_back(iFJ);

        // now clean the leptons from these good fat topjets

        //  mark here the ones that aren't good!
        for (int iL = 0, nL = *nLep_; iL < nL; ++iL) {
          if (!(sel_leps.get()[iL])) { // || sel_leps_extrafortau.get()[iL])) { 
            bad_leptons.push_back(iL); 
            continue; // not sure that second statement after the or is...
          }
          if (reco::deltaR2((*Lep_eta_)[iL], (*Lep_phi_)[iL], (*FJet_eta_)[iFJ], (*FJet_phi_)[iFJ]) < deltaR2cut_fatjets) {
            bad_leptons.push_back(iL);
            //std::cout << "found a bad lepton with index " << iL << std::endl;
            // clean_leps_.push_back(iL);
            // _cl.push_back(iL);
          }
        }
      } 
      else {
        other_fatjets.push_back(iFJ);
      }

    } // cleaned the top fat jets from the leptons


//std::cout<<__LINE__<<std::endl;
    //std::cout << " this is the size of bad_leptons " << bad_leptons.size() << std::endl;
    //std::cout << " this is the size of tp jets " << _ctj->size() << std::endl;

    // put all clean leptons into a vector ...
    for (int iL = 0, nL = *nLep_; iL < nL; ++iL) {

      bool lepIsGood = true;
      if (!(sel_leps.get()[iL] )) lepIsGood = false; //|| sel_leps_extrafortau.get()[iL])) continue; // not sure that second statement after the or is...

      for (std::vector<int>::iterator iBL = bad_leptons.begin(); iBL != bad_leptons.end(); ++iBL) {
        if (*iBL == iL) lepIsGood = false;
      }

      if (lepIsGood) {_cl->push_back(iL); clean_leps_.push_back(iL);}
    }

    //std::cout << " this is the size of good_leptons " << _cl->size() << std::endl;

//std::cout<<__LINE__<<std::endl;

    
    // loop on the rest of the fatjets
    for (std::vector<int>::iterator iFJ = other_fatjets.begin(); iFJ != other_fatjets.end(); iFJ++){
      float tagger_w = (*FJet_deepTag_WvsQCD_)[*iFJ];
      float tagger_z = (*FJet_deepTag_ZvsQCD_)[*iFJ];

      bool jetIsGood = true;

      for(std::vector<int>::iterator iL = _cl->begin(); iL != _cl->end(); ++iL) {
        float tmp_dr =  reco::deltaR2((*Lep_eta_)[*iL], (*Lep_phi_)[*iL], (*FJet_eta_)[*iFJ], (*FJet_phi_)[*iFJ]);

        if (tmp_dr < deltaR2cut_fatjets) jetIsGood = false;
      }

      if (!jetIsGood) continue;

      if( (tagger_w > 0.762 || tagger_z > 0.479) ) {
        clean_vjets_.push_back(*iFJ); // tag-values for 2.5% mistag
        _cvj->push_back(*iFJ); // tag-values for 2.5% mistag
      }
      else {
        clean_fjets_.push_back(*iFJ);
        _cfj->push_back(*iFJ);
      }
    }

    for (int iJ = 0, nJ = *nJet_; iJ < nJ; ++iJ) {
      if ( good[iJ] && sel_jets[iJ]) {
        _cj->push_back(iJ);
        if(fabs((*Jet_eta_)[iJ]) < 2.4) clean_jets_.push_back(iJ); // only to count fwd jets. marc: i don't get it...
        bool lightjet = true;

        for (int iFJ = 0, nFJ = *nFJet_; iFJ < nFJ; ++iFJ) { // clean all the jets mathcing a fatjet to count as a light jet
          float dr2 = reco::deltaR2((*Jet_eta_)[iJ],(*Jet_phi_)[iJ],(*FJet_eta_)[iFJ],(*FJet_phi_)[iFJ]);
          if (goodfj[iFJ] && dr2<deltaR2cut_fatjets)  
            lightjet = false;
        }

        for(std::vector<int>::iterator iL = _cl->begin(); iL != _cl->end(); ++iL) { // clean all the regular jets with the leptons
          float dr2 = reco::deltaR2((*Jet_eta_)[iJ],(*Jet_phi_)[iJ],(*Lep_eta_)[*iL],(*Lep_phi_)[*iL]);
          if (dr2 < deltaR2cut)  lightjet = false;
        }

        if (lightjet ) {
          _clj->push_back(iJ);
          clean_ljets_.push_back(iJ);
        }
      }
    }

//std::cout<<__LINE__<<std::endl;
    std::vector<std::vector<int>* > ret;
    ret.push_back(_cl .get());
    ret.push_back(_ct .get());
    ret.push_back(_cj .get());
    ret.push_back(_ctj.get());
    ret.push_back(_cvj.get());
    ret.push_back(_cfj.get());
    ret.push_back(_clj.get());
//std::cout<<__LINE__<<std::endl;
    return ret;
  }

// end marc


private:
  std::unique_ptr<bool[]> sel_leps, sel_leps_extrafortau, sel_taus, sel_jets;
  CollectionSkimmer &clean_leps_, &clean_taus_, &clean_jets_, &clean_tjets_, &clean_vjets_, &clean_fjets_, &clean_ljets_;
  rcount nLep_, nTau_, nJet_, nFJet_, nLJet_;
  rfloats *Lep_pt_, *Lep_eta_, *Lep_phi_;
  rfloats *Tau_pt_, *Tau_eta_, *Tau_phi_;
  rfloats *Jet_pt_, *Jet_phi_, *Jet_eta_, *Jet_btagCSV_;
  rfloats *FJet_pt_, *FJet_phi_, *FJet_eta_, *FJet_btagCSV_, *FJet_deepTag_TvsQCD_, *FJet_deepTag_WvsQCD_, *FJet_deepTag_ZvsQCD_;
  vector<rfloats*> Jet_corr_;
  rints    *Lep_jet_, *Tau_jet_;
  float deltaR2cut;
  float deltaR2cut_taus;
  float deltaR2cut_fatjets;
  std::set<int> _jetptcuts;
  std::unique_ptr<std::vector<int> > _cl;
  std::unique_ptr<std::vector<int> > _ct;
  std::unique_ptr<std::vector<int> > _cj;
  std::unique_ptr<std::vector<int> > _ctj;
  std::unique_ptr<std::vector<int> > _cvj;
  std::unique_ptr<std::vector<int> > _cfj;
  std::unique_ptr<std::vector<int> > _clj;
  bool cleanJetsWithFOTaus_;
  float bTagL_,bTagM_;
  bool cleanWithRef_;
  float fwdJetPt1_, fwdJetPt2_;
};
