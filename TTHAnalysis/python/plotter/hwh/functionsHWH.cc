#include "TFile.h"
#include "TH2.h"
#include "TH2Poly.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>

int hwh_catIndex_2lss(int LepGood1_pdgId, int LepGood2_pdgId, float tjvv, float ttv, float tt, float other)
{

//2lss_ee_tjvv   1
//2lss_ee_other  2
//2lss_ee_ttv    3
//2lss_ee_tt     4
//2lss_em_tjvv   5
//2lss_em_other  6
//2lss_em_ttv    7
//2lss_em_tt     8
//2lss_mm_tjvv   9
//2lss_mm_other  10
//2lss_mm_ttv    11
//2lss_mm_tt     12
  int flch = 0;
  int procch = 0;

  if (abs(LepGood1_pdgId)+abs(LepGood2_pdgId) == 22)
    flch = 0;
  else if (abs(LepGood1_pdgId)+abs(LepGood2_pdgId) == 24)
    flch = 1;
  else if (abs(LepGood1_pdgId)+abs(LepGood2_pdgId) == 26)
    flch = 2;
  else
    cout << "[2lss]: It shouldnt be here. pdgids are " << abs(LepGood1_pdgId) << " " << abs(LepGood2_pdgId)  << endl;

  if (tjvv >= ttv && tjvv >= tt && tjvv >= other)
    procch = 0;
  else if (other >= tjvv && other >= ttv && other >= tt)
    procch = 1;
  else if (ttv >= tjvv && ttv >= other && ttv >= tt)
    procch = 2;
  else if (tt >= tjvv && tt >= other && tt >= ttv)
    procch = 3;
  else 
    cout << "[2lss]: It shouldnt be here. DNN scores are " << tjvv << " " << other << " " << ttv << " " << tt << endl;
      
  return flch*4+procch+1;

}

std::map<TString,int> hwhbins2lss = {{"ee_tjVVnode",5},{"ee_Othernode",8},{"ee_ttVnode",6},{"ee_ttnode",4},
				  {"em_tjVVnode",13},{"em_Othernode",8},{"em_ttVnode",19},{"em_ttnode",11},
				  {"mm_tjVVnode",13},{"mm_Othernode",11},{"mm_ttVnode",15},{"mm_ttnode",7}};
std::vector<TString> hwhbin2lsslabels = {
  "ee_tjVVnode","ee_Othernode","ee_ttVnode","ee_ttnode",
  "em_tjVVnode","em_Othernode","em_ttVnode","em_ttnode",
  "mm_tjVVnode","mm_Othernode","mm_ttVnode","mm_ttnode"
};

std::map<TString, TH1F*> hwhbinHistos2lss;
std::map<TString, int> hwhbins2lsscumul;
TFile* f2lssHwhbins;

int hwh_catIndex_2lss_MVA(int LepGood1_pdgId, int LepGood2_pdgId, float tjvv, float ttv, float tt, float other)
{
  if (!f2lssHwhbins){
    int offset = 0;
    f2lssHwhbins = TFile::Open("../../data/kinMVA/DNNBin_hwh_2lss.root");
    for (auto & la : hwhbin2lsslabels){
      int bins = hwhbins2lss[la];
      hwhbinHistos2lss[la] = (TH1F*) f2lssHwhbins->Get(Form("%s_2018_Map_nBin%d", la.Data(), bins));
      hwhbins2lsscumul[la] = offset;
      offset += bins;
    }
  }
  
  int idx = hwh_catIndex_2lss(LepGood1_pdgId, LepGood2_pdgId, tjvv,ttv, tt,other); 
  
  TString binLabel = hwhbin2lsslabels[idx-1];
  
  float mvavar = 0;
  if (tjvv >= ttv && tjvv >= tt && tjvv >= other)
    mvavar = tjvv;
  else if (other >= tjvv && other >= ttv && other >= tt)
    mvavar =other;
  else if (ttv >= tjvv && ttv >= other && ttv >= tt)
    mvavar = ttv;
  else if (tt >= tjvv && tt >= other && tt >= ttv)
    mvavar = tt;
  else 
    cout << "It shouldnt be here" << endl;
  return hwhbinHistos2lss[binLabel]->FindBin( mvavar ) + hwhbins2lsscumul[binLabel];

}


// for plots

int hwh_2lss_node( float tjvv, float ttv, float tt, float other ){

  int procch = 0;
  if (tjvv >= ttv && tjvv >= tt && tjvv >= other)
    procch = 0;
  else if (other >= tjvv && other >= ttv && other >= tt)
    procch = 1;
  else if (ttv >= tjvv && ttv >= other && ttv >= tt)
    procch = 2;
  else if (tt >= tjvv && tt >= other && tt >= ttv)
    procch = 3;
  else 
    cout << "[2lss]: It shouldnt be here. DNN scores are " << tjvv << " " << other << " " << ttv << " " << tt << endl;

  return procch;
}


std::vector<TString> hwhbin2lsslabels_plots = {
  "ee_tjVVnode" , "em_tjVVnode" ,  "mm_tjVVnode", 
  "ee_Othernode", "em_Othernode",  "mm_Othernode",
  "ee_ttVnode" , "em_ttVnode" ,  "mm_ttVnode",
  "ee_ttnode" , "em_ttnode" ,  "mm_ttnode",

};

std::map<TString, TH1F*> hwhbinHistos2lss_plots;
TFile* f2lssHwhbins_plots;


int hwh_catIndex_2lss_plots(int LepGood1_pdgId, int LepGood2_pdgId, float tjvv, float ttv, float tt, float other)
{

  if (!f2lssHwhbins_plots){
    f2lssHwhbins_plots = TFile::Open("../../data/kinMVA/DNNBin_v3_xmas.root");
    for (auto & la : hwhbin2lsslabels_plots){
      int bins = hwhbins2lss[la];
      hwhbinHistos2lss_plots[la] = (TH1F*) f2lssHwhbins_plots->Get(Form("%s_2018_Map_nBin%d", la.Data(), bins));
    }
  }

  int idx = hwh_catIndex_2lss(LepGood1_pdgId, LepGood2_pdgId, tjvv,ttv, tt,other); 
  TString binLabel = hwhbin2lsslabels[idx-1];
  int offset=0;
  int node = hwh_2lss_node(tjvv, ttv,tt, other);
  if (abs(LepGood1_pdgId*LepGood2_pdgId) == 143){
    if (node == 0) offset = 5;
    else if (node == 1) offset = 8;
    else if (node == 2) offset = 6;
    else offset = 4;
  }
  if (abs(LepGood1_pdgId*LepGood2_pdgId) == 169){
    if (node == 0) offset = 5+13;
    else if (node == 1) offset = 8+8;
    else if (node == 2) offset = 6+19;
    else offset = 4+11;
  }

  float mvavar = 0;
  if (tjvv >= ttv && tjvv >= tt && tjvv >= other)
    mvavar = tjvv;
  else if (other >= tjvv && other >= ttv && other >= tt)
    mvavar =other;
  else if (ttv >= tjvv && ttv >= other && ttv >= tt)
    mvavar = ttv;
  else if (tt >= tjvv && tt >= other && tt >= ttv)
    mvavar = tt;
  else 
    cout << "It shouldnt be here" << endl;

  return hwhbinHistos2lss_plots[binLabel]->FindBin( mvavar ) + offset;
    
}


int hwh_catIndex_3l(float tjvv, float tt, float other, int lep1_pdgId, int lep2_pdgId, int lep3_pdgId, float mZ1 )
{

  int sumpdgId = abs(lep1_pdgId)+abs(lep2_pdgId)+abs(lep3_pdgId);
  int flch = 0;
  int procch = 0;

  if (tjvv >= tt && tjvv >= other) {
    if (abs(mZ1-91.2) > 10)
      return 1; // VV_nr
    else
      return 2; // VV_r
  }
  else {
    if (sumpdgId == 33){ // eee 
      flch = 0;
    }
    else if (sumpdgId == 35){ // eem
      flch = 1;
    }
    else if (sumpdgId == 37){ // emm
      flch = 2;
    }
    else if (sumpdgId == 39){ // mmm
      flch = 3;
    }
    else
      cout << "[3l]: It shouldnt be here. pdgids are " << abs(lep1_pdgId) << " " <<  abs(lep2_pdgId) << " " << abs(lep3_pdgId)  << endl;
    return flch+3;
  }

  cout << "[hwh_catIndex_3l]: It should not be here" << endl;
  return -1;

}

std::vector<TString> hwhbin3llabels = {
  "tjVV_nr","tjVV_r","eee_Bkgnode","eem_Bkgnode","emm_Bkgnode","mmm_Bkgnode" 
};

std::map<TString, TH1F*> hwhbinHistos3l;
std::map<TString, int> hwhbins3lcumul;
TFile* f3lHwhBins;

int hwh_catIndex_3l_MVA(float tjvv, float tt, float other, int lep1_pdgId, int lep2_pdgId, int lep3_pdgId, float mZ1 )
{

  if (!f3lHwhBins){
    // for the moment take for all the same binning: [0,0.43,0.47,0.5,0.55,0.61,0.71,1.0]
    f3lHwhBins=TFile::Open("../../data/kinMVA/binning_3l.root");
    int count=0;
    for (auto label : hwhbin3llabels){
      hwhbinHistos3l[label] = (TH1F*) f3lHwhBins->Get("tH_bl"); //(label);
      hwhbins3lcumul[label] = count;
      count += hwhbinHistos3l[label]->GetNbinsX();
    }
  }
  TString binLabel = hwhbin3llabels[hwh_catIndex_3l(tjvv,tt,other,lep1_pdgId,lep2_pdgId,lep3_pdgId,mZ1)-1];
  float mvas[] = { tjvv, tt, other };
  float mvavar = *std::max_element( mvas, mvas+3 );
  return hwhbinHistos3l[binLabel]->FindBin( mvavar ) + hwhbins3lcumul[binLabel];
  
  cout << "[hwh_catIndex_3l_MVA]: It should not be here "<< tjvv << " " << tt << " " << other << endl;
  return -1;

}


int hwh_catIndex_3l_node(float tjvv, float tt, float other ){
  int procch = 0;
  // for the time being, dnn confuses other and tt. So merge these two cats
  if (tjvv >= tt && tjvv >= other)
    procch = 0;
  else
    procch = 1;
  return procch;
}

float hwh_mva_4l(float tjvv, float tt, float other ){
  return tjvv / (tjvv+tt+other);
}

int hwh_catIndex_4l_MVA(float tjvv, float tt, float other, float cut=0.4){
  if (hwh_mva_4l(tjvv,tt,other) < cut) return 1;
  else return 2;
}
