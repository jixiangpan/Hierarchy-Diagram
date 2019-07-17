///
/// Created by Xiangpan Ji at BNL, 07/10/2019
///
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<set>
#include<string>
using namespace std;

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"

#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "TROOT.h"
#include "TString.h"
#include "TMath.h"
#include "TGaxis.h"
#include "THStack.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"
#include "TExec.h"
#include "TDatabasePDG.h"

map<int, TString> map_pdg2name;

TGraph *graph_all_TreeDiagram;
map<int,int> graph_map_pdg;
void draw_text_TreeDiagram()
{
  Int_t i,n;
  Double_t x,y;
  TLatex *l;

  // TDatabasePDG *database_pdg = new TDatabasePDG();
  
  TGraph *g = graph_all_TreeDiagram;
  n = g->GetN();
  for (i=0; i<n; i++) {
    g->GetPoint(i,x,y);
    int pdg = graph_map_pdg[i];
    // if( database_pdg->GetParticle(pdg)!=NULL ) 
    //   cout<<" Particle name "<<pdg<<" "<<database_pdg->GetParticle(pdg)->GetName()<<endl;

    if( map_pdg2name.find(pdg)!=map_pdg2name.end() ) {
      TString name = map_pdg2name[pdg];
      l = new TLatex(x,y+0.2, name);
    }
    else {
      if(pdg>2e9) {
	l = new TLatex(x,y+0.2,"pseudo");
      }
      else {
	l = new TLatex(x,y+0.2,TString::Format("%d",pdg));
      }
    }
    
    // l = new TLatex(x,y+0.2,Form("%d",pdg));
    // l->SetTextSize(0.025);
    l->SetTextSize(0.03);
    l->SetTextFont(42);
    l->SetTextAlign(21);
    l->Paint();
  }
}

/////////////////////////////////////////////////////////////////////

/// gen: generation

int final_yValue_Loop_func = 0;
void Loop_func(set<int> set_root_trk, map<int,vector<int> > *map_trk_daughters, map<int,int> *map_parentrkid, int current_trk, map<int,int> *par_upper, map<int,int> *par_index, map<int,int> *out_map_final_yValue, map<int,int> *out_map_trk_generation )
{
  ///      set_root_trk: set of root trks
  /// map_trk_daughters: data, track and daughters
  ///       current_trk: current trk
  ///         par_upper: number of daughters of current trk
  ///         par_index: loop index of the daughters
  
  int user_trk = current_trk;
  if( set_root_trk.find(current_trk)!=set_root_trk.end() ) {// root_trk
    (*out_map_trk_generation)[user_trk]++;
    par_upper->clear();
    par_index->clear();
  }

  (*par_upper)[user_trk] = (*map_trk_daughters)[user_trk].size();
  (*par_index)[user_trk]++;

  if( (*map_trk_daughters)[user_trk].size()==0 ) {
    final_yValue_Loop_func++;
    (*out_map_final_yValue)[user_trk] = final_yValue_Loop_func;

    if( set_root_trk.find(user_trk)==set_root_trk.end() ) {
      int parentrkid = (*map_parentrkid)[user_trk];
      int parent_generation = (*out_map_trk_generation)[parentrkid];
      (*out_map_trk_generation)[user_trk] = parent_generation+1;
    }

    int current_generation = (*out_map_trk_generation)[user_trk];
    cout<<TString::Format(" ------> final %4d, generation %4d, yValue %4d", user_trk, current_generation, final_yValue_Loop_func)<<endl;
  }
  else {
    for( (*par_index)[user_trk]=0; (*par_index)[user_trk]<(*par_upper)[user_trk]; (*par_index)[user_trk]++ ) {
      current_trk = (*map_trk_daughters)[user_trk].at( (*par_index)[user_trk] );
      
      if( set_root_trk.find(current_trk)==set_root_trk.end() ) {
	int parentrkid = (*map_parentrkid)[current_trk];
	int parent_generation = (*out_map_trk_generation)[parentrkid];
	(*out_map_trk_generation)[current_trk] = parent_generation+1;
      }
      
      Loop_func(set_root_trk, map_trk_daughters, map_parentrkid, current_trk, par_upper, par_index, out_map_final_yValue, out_map_trk_generation );
    }
  }
  
}

/////////////////////////////////////////////////////////////////////
//////////////////////////// MAIN ///////////////////////////////////
/////////////////////////////////////////////////////////////////////

void read_fa_TreeDiagram_color(int outflag, int file_bgn, int file_end)
{
  cout<<endl<<" Hello, World!"<<endl<<endl;
  
  gROOT->ProcessLine(".x DrawOption.cc");
  TString roostr = "";

  double SizeMarker = 1.4;
  
  const int proton  = 2212;
  const int neutron = 2112;

  int color_proton = kBlue;
  int color_neutron = kRed;
  
  map_pdg2name[11] = "e";
  map_pdg2name[12] = "#nu_{e}";
  map_pdg2name[13] = "#mu^{-}";
  map_pdg2name[14] = "#nu_{#mu}";
  map_pdg2name[111] = "#pi^{0}";
  map_pdg2name[211] = "#pi^{+}";
  map_pdg2name[-211] = "#pi^{-}";
  map_pdg2name[2212] = "p";
  map_pdg2name[2112] = "n";
  map_pdg2name[2224] = "#Delta^{++}";
  map_pdg2name[2214] = "#Delta^{+}";
  map_pdg2name[2114] = "#Delta^{0}";
  map_pdg2name[1114] = "#Delta^{-}";
  map_pdg2name[3122] = "#Lambda";
  map_pdg2name[3222] = "#Sigma^{+}";
  map_pdg2name[3212] = "#Sigma^{0}";
  map_pdg2name[3112] = "#Sigma^{-}";
  map_pdg2name[3334] = "#Omega^{-}";
  map_pdg2name[1000180400] = "{}^{40}Ar";
  map_pdg2name[1000180390] = "{}^{39}Ar";
  map_pdg2name[1000170390] = "{}^{39}Cl";

  map_pdg2name[3334] = "#Omega^{-}";
  
  ////////////////////////////// input files to be analysed

  TString str_inputlist = "info_edit_temp.txt";
  
  ifstream read_list(str_inputlist, ios::in);
  if(!read_list) {
    cerr<<"Error: No "<<str_inputlist<<endl;
    exit(1);
  }

  /// Initializing data
  vector<int>     *mcp_pdg = 0;// MC particles
  vector<int>     *mcp_trkid = 0;// Track_ID
  vector<int>     *mcp_parentrkid = 0;// Parent's Track_ID
  Int_t           nu_PdgCode;
  Int_t           nu_CCNC;
  Int_t           gt_Gscatter;

  TBranch        *b_nu_PdgCode;   //!
  TBranch        *b_nu_CCNC;   //!
  TBranch        *b_gt_Gscatter;   //!
  
  TChain *tree = new TChain("GetInfo/fTree");
  //tree->Add("read_fa_data.root");
    
  int line_file = 0;
  for(int ifile=1; ifile<=50000; ifile++) {
    line_file++;
    read_list>>roostr;

    if( (line_file<file_bgn) || (line_file>file_end) ) continue;
    ifstream read_check(roostr, ios::in); if(!read_check) continue;
    
    cout<<TString::Format(" ---> %4d   ", ifile)<<roostr<<endl;
    tree->Add( roostr );
  }
  
  // Set branch addresses and branch pointers
  tree->SetBranchAddress("nu_PdgCode", &nu_PdgCode, &b_nu_PdgCode);
  tree->SetBranchAddress("mcp_pdg",        &mcp_pdg);
  tree->SetBranchAddress("mcp_trkid",      &mcp_trkid);
  tree->SetBranchAddress("mcp_parentrkid", &mcp_parentrkid);
  tree->SetBranchAddress("nu_CCNC",        &nu_CCNC, &b_nu_CCNC);
  tree->SetBranchAddress("gt_Gscatter",    &gt_Gscatter, &b_gt_Gscatter);
  
  int entries = tree->GetEntries();
  cout<<endl<<" Entries: "<<entries<<endl<<endl;

  map<int,int> map_trkid;
  map<int,int> map_parentrkid;
  map<int,int> map_pdg;

  /////
  const int flag_numu = 14;
  const int flag_nue  = 12;
  const int flag_numu_bar = -14;
  const int flag_nue_bar  = -12;

  int used_nu = flag_nue;

  int flag_find_interaction = false;
  
  for(int ientry=0; ientry<entries; ientry++) {
    tree->GetEntry( ientry );
    
    if( nu_PdgCode!=used_nu ) continue;
    if( nu_CCNC!=0 ) continue;// CC-0, NC-1
    if( gt_Gscatter!=1 ) continue;
    
    flag_find_interaction = true;
      
    // http://internal.dunescience.org/doxygen/classgenie_1_1ScatteringType.html  #details
    // gt_Gscatter:
    // case(kScQuasiElastic) :            return "QES";       1
    // case(kScSingleKaon) :              return "1Kaon";     2
    // case(kScDeepInelastic) :           return "DIS";       3
    // case(kScResonant) :                return "RES";       4
    // case(kScCoherent) :                return "COH";       5
    // case(kScDiffractive) :             return "DFR";       6
    // case(kScNuElectronElastic) :       return "NuEEL";     7
    // case(kScInverseMuDecay) :          return "IMD";       8
    // case(kScAMNuGamma) :               return "AMNuGamma"; 9
    // case(kScMEC) :                     return "MEC";       10
    // case(kScCoherentElas) :            return "COHEl";     11
    // case(kScInverseBetaDecay) :        return "IBD";       12
    // case(kScGlashowResonance) :        return "GLR";       13
    // case(kScIMDAnnihilation) :         return "IMDAnh";    14
    // case(kScDarkMatterElastic) :       return "DME";       101
    // case(kScDarkMatterDeepInelastic) : return "DMDIS";     102
    // default :                          return "Unknown";   0  
    
    for(unsigned int idx=0; idx<mcp_pdg->size(); idx++) {
      int trkid      = mcp_trkid->at(idx);
      int parentrkid = mcp_parentrkid->at(idx);
      int pdg        = mcp_pdg->at(idx);
      map_trkid[trkid]      = trkid;
      map_parentrkid[trkid] = parentrkid;
      map_pdg[trkid]        = pdg;
    }

    if( flag_find_interaction ) break;
    
  }

  map<int,int>::iterator it_map_trkid;
  map<int, vector<int> > map_trk_daughters;
  set<int> set_root_trk;
  set<int> set_final_trk;
  for( it_map_trkid=map_trkid.begin(); it_map_trkid!=map_trkid.end(); it_map_trkid++ ) {
    int trkid = it_map_trkid->first;
    int parentrkid = map_parentrkid[trkid];
    if( parentrkid<0 ) {
      set_root_trk.insert(trkid);
      continue;
    }
    map_trk_daughters[parentrkid].push_back(trkid);
  }
  
  for( it_map_trkid=map_trkid.begin(); it_map_trkid!=map_trkid.end(); it_map_trkid++ ) {
    int trkid = it_map_trkid->first;
    if( map_trk_daughters[trkid].size()==0 ) set_final_trk.insert(trkid);
  }
  
  /// Tree Diagram
  set<int>::iterator it_set_root_trk;
  map<int, int> par_upper;
  map<int, int> par_index;
  map<int, int> map_final_yValue;
  map<int, int> map_trk_generation;
  map<int, double> map_trk_yValue;
  
  cout<<endl;
  for( it_set_root_trk=set_root_trk.begin(); it_set_root_trk!=set_root_trk.end(); it_set_root_trk++ ) {
    int root_trk = (*it_set_root_trk);
    Loop_func( set_root_trk, &map_trk_daughters, &map_parentrkid, root_trk, &par_upper, &par_index, &map_final_yValue, &map_trk_generation);
  }
  cout<<endl;

  int max_generation = 0;
  int max_yValue     = 0;
  set<int>::iterator it_set_final_trk;
  for( it_set_final_trk=set_final_trk.begin(); it_set_final_trk!=set_final_trk.end(); it_set_final_trk++ ) {
    int trkid = (*it_set_final_trk);
    int generation = map_trk_generation[trkid];
    int int_yVlaue = map_final_yValue[trkid];
    if( max_generation<generation ) max_generation = generation;
    if( max_yValue<int_yVlaue ) max_yValue = int_yVlaue;
  }

  for(int iG=max_generation; iG>=1; iG--) {
    for( it_map_trkid=map_trkid.begin(); it_map_trkid!=map_trkid.end(); it_map_trkid++ ) {
      int current_trkid = it_map_trkid->first;
      if( map_trk_generation[current_trkid]==iG ) {
	int Ndaugther = map_trk_daughters[current_trkid].size();
	if( Ndaugther==0 ) {
	  map_trk_yValue[current_trkid] = map_final_yValue[current_trkid] *1.;
	}
	else {
	  for(int iD=0; iD<Ndaugther; iD++) {
	    int daughter_trk = map_trk_daughters[current_trkid].at(iD);
	    map_trk_yValue[current_trkid] += map_trk_yValue[daughter_trk];
	  }
	  map_trk_yValue[current_trkid] = map_trk_yValue[current_trkid]/Ndaugther;
	}
      }
    }
  }

  /// check
  cout<<endl;
  for( it_set_final_trk=set_final_trk.begin(); it_set_final_trk!=set_final_trk.end(); it_set_final_trk++ ) {
      int trkid = (*it_set_final_trk);
      int pdg = map_pdg[trkid];
      int parentrkid = map_parentrkid[trkid];
      cout<<TString::Format(" %2d (%d)", trkid, pdg);
      while( parentrkid>=0 ) {
	int current_trkid = parentrkid;
	int current_pdg   = map_pdg[current_trkid];	 
	int current_parent_trkid = map_parentrkid[current_trkid];
	int current_parent_pdg   = map_pdg[current_parent_trkid];
	parentrkid = current_parent_trkid;	
	cout<<TString::Format(" ---> %d (%d)", current_trkid, current_pdg);
      }
      cout<<endl;
  }
  cout<<endl;

  // cout<<endl;
  // for( it_map_trkid=map_trkid.begin(); it_map_trkid!=map_trkid.end(); it_map_trkid++ ) {
  //   int trkid = it_map_trkid->first;
  //   double generationValue = map_trk_generation[trkid] * 1.;
  //   double yValue = map_trk_yValue[trkid];
  //   cout<<TString::Format( "trk/G: %4d %2.0f", trkid, generationValue)<<endl;
  // }
  // cout<<endl;
  
  // TGraph *graph_all_TreeDiagram = new TGraph();
  graph_all_TreeDiagram = new TGraph();
  graph_all_TreeDiagram->SetName("graph_all_TreeDiagram");
  int line_graph_all_TreeDiagram = 0;
  // map<int,int> graph_map_pdg;
  for( it_map_trkid=map_trkid.begin(); it_map_trkid!=map_trkid.end(); it_map_trkid++ ) {
    int trkid = it_map_trkid->first;
    double generationValue = map_trk_generation[trkid] * 1.;
    double yValue = map_trk_yValue[trkid];

    line_graph_all_TreeDiagram++;
    graph_all_TreeDiagram->SetPoint(line_graph_all_TreeDiagram-1, generationValue, yValue);
    graph_map_pdg[line_graph_all_TreeDiagram-1] = map_pdg[trkid];
  }
  

  roostr = "canv_graph_TreeDiagram";
  TCanvas *canv_graph_TreeDiagram = new TCanvas(roostr, roostr, 900, 650);
  canv_graph_TreeDiagram->SetLeftMargin(0.14);  
  canv_graph_TreeDiagram->SetRightMargin(0.05);
  canv_graph_TreeDiagram->SetTopMargin(0.09);
  canv_graph_TreeDiagram->SetBottomMargin(0.18);

  roostr = "h2_basic_graph_TreeDiagram";
  TH2D *h2_basic_graph_TreeDiagram = new TH2D(roostr, roostr,
					      max_generation+1, -0.5, max_generation+0.5,
					      max_yValue+1+1, -0.5, max_yValue+0.5+1);
  h2_basic_graph_TreeDiagram->Draw();
  h2_basic_graph_TreeDiagram->SetStats(0);
  h2_basic_graph_TreeDiagram->SetXTitle("Generation");
  h2_basic_graph_TreeDiagram->GetXaxis()->CenterTitle();
  h2_basic_graph_TreeDiagram->GetYaxis()->CenterTitle();
  h2_basic_graph_TreeDiagram->GetXaxis()->SetNdivisions(508);
  h2_basic_graph_TreeDiagram->GetYaxis()->SetNdivisions(508);

  TExec *ex_TreeDiagram = new TExec("ex_TreeDiagram","draw_text_TreeDiagram();");
  graph_all_TreeDiagram->GetListOfFunctions()->Add(ex_TreeDiagram);
  graph_all_TreeDiagram->SetMarkerStyle(20);
  graph_all_TreeDiagram->SetMarkerColor(kBlack);
  graph_all_TreeDiagram->SetMarkerSize(SizeMarker);
  graph_all_TreeDiagram->Draw("p");

  ///
  map<int, TGraph*> map_graph_trks;
  for( it_map_trkid=map_trkid.begin(); it_map_trkid!=map_trkid.end(); it_map_trkid++ ) {
    int trk = it_map_trkid->first;
    int parent_trk = map_parentrkid[trk];

    if( parent_trk<0 ) continue;
    
    double trk_H = map_trk_generation[trk];
    double trk_Y = map_trk_yValue[trk];

    double parent_trk_H = map_trk_generation[parent_trk];
    double parent_trk_Y = map_trk_yValue[parent_trk];
    
    roostr = TString::Format("map_graph_trks_%3d", trk);
    map_graph_trks[trk] = new TGraph();
    map_graph_trks[trk]->SetName(roostr);
    map_graph_trks[trk]->SetPoint(0, parent_trk_H, parent_trk_Y);
    map_graph_trks[trk]->SetPoint(1, trk_H, trk_Y);
    map_graph_trks[trk]->Draw("pl");
  }

  ///
  map<int, TGraph*> map_graph_particles;
  for( it_map_trkid=map_trkid.begin(); it_map_trkid!=map_trkid.end(); it_map_trkid++ ) {
    int trk = it_map_trkid->first;
    int parent_trk = map_parentrkid[trk];
    int pdg = map_pdg[trk];
    
    double trk_H = map_trk_generation[trk];
    double trk_Y = map_trk_yValue[trk];

    if( pdg==proton ) {
      roostr = TString::Format("map_graph_particles_%3d", trk);
      map_graph_particles[trk] = new TGraph();
      map_graph_particles[trk]->SetName(roostr);
      map_graph_particles[trk]->SetPoint(0, trk_H, trk_Y);
      map_graph_particles[trk]->Draw("p");
      map_graph_particles[trk]->SetMarkerStyle(20);
      map_graph_particles[trk]->SetMarkerSize(SizeMarker);
      map_graph_particles[trk]->SetMarkerColor(color_proton);
    }
  
    if( pdg==proton ) {
      roostr = TString::Format("map_graph_particles_%3d", trk);
      map_graph_particles[trk] = new TGraph();
      map_graph_particles[trk]->SetName(roostr);
      map_graph_particles[trk]->SetPoint(0, trk_H, trk_Y);
      map_graph_particles[trk]->Draw("p");
      map_graph_particles[trk]->SetMarkerStyle(20);
      map_graph_particles[trk]->SetMarkerSize(SizeMarker);
      map_graph_particles[trk]->SetMarkerColor(color_proton);
    }
    
    if( pdg==neutron ) {
      roostr = TString::Format("map_graph_particles_%3d", trk);
      map_graph_particles[trk] = new TGraph();
      map_graph_particles[trk]->SetName(roostr);
      map_graph_particles[trk]->SetPoint(0, trk_H, trk_Y);
      map_graph_particles[trk]->Draw("p");
      map_graph_particles[trk]->SetMarkerStyle(20);
      map_graph_particles[trk]->SetMarkerSize(SizeMarker);
      map_graph_particles[trk]->SetMarkerColor(color_neutron);
    }
    
  }

  
}

