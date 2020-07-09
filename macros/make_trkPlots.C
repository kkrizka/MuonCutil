#include "TFile.h"
#include "TTree.h"
#include "TEfficiency.h"
#include "TH1F.h"

/*
 =====================================================================================================================

 Simple macro to plot tracks efficiency. It produces:
    - Tracks efficiency vs pT
    - Tracks efficiency vs theta
    - Tracks efficiency vs phi
    - pT resoltution vs pT 
    - pT resoltution vs thetaT 

 If required, set debug = 1, plot also the base variables: pt, phi, theta, chi2 and nhit 
 If required, set saveonfile = 1, saved also the plots in png files.  

 Use: root -l make_trkPlots.C

   Copyright (C) 2020 Massimo Casarsa
 
 =====================================================================================================================
*/

// save plots on file?
const bool saveonfile = 1;
// plot also the base variables: pt, phi, theta, chi2 and nhit
const bool debug = 1;


TFile* infile;
TTree* T_eff;
TTree* T_trk;

TEfficiency* eff_pt;
TEfficiency* eff_phi;
TEfficiency* eff_theta;

TH1F* h_pt_reco;
TH1F* h_phi_reco;
TH1F* h_theta_reco;
TH1F* h_chi2_reco;
TH1F* h_nhit_reco;

const unsigned int Nbins_pt = 10;
TH1F* h_pt_res[Nbins_pt];

const unsigned int Nbins_theta = 20;
TH1F* h_theta_res[Nbins_theta];

TH1F h_pt_res_tot("h_pt_res_tot",";p_{T} [GeV];#sigma(#Deltap_{T}/p_{T}^{2}) [GeV^{-1}]",Nbins_pt,0.,10.);
TH1F h_theta_res_tot("h_theta_res_tot",";#theta [#circ];#sigma(#Deltap_{T}/p_{T}^{2}) [GeV^{-1}]",Nbins_theta,0.,90.);

TCanvas* c[10];
TCanvas* basevar[5];

// =====================================================================================================================

void make_trkPlots(TString filename="histograms.root"){

  // --- Book histograms:
  
  eff_pt = new TEfficiency("eff_pt",";p_{T} [GeV];#epsilon_{trk}",40,0.,10.);
  eff_phi = new TEfficiency("eff_phi",";#phi [#circ];#epsilon_{trk}",180,-180.,180.);
  eff_theta = new TEfficiency("eff_theta",";#theta [#circ];#epsilon_{trk}",60,0.,180.);

  h_pt_reco = new TH1F("h_pt_reco",";p_{T} [GeV]",100,0.,10.) ;
  h_phi_reco = new TH1F("h_phi_reco",";#phi [#circ]",120,-180.,180.) ;
  h_theta_reco = new TH1F("h_theta_reco",";#theta [#circ]",100,20.,150.) ;
  h_chi2_reco = new TH1F("h_chi2_reco",";#chi^{2}/ndf",100,0.,20.) ;
  h_nhit_reco = new TH1F("h_nhit_reco",";N_{hit}",20,0.,20.) ;

  for ( unsigned int ih=0; ih<h_pt_res_tot.GetNbinsX(); ++ih ) {
    TString hname = Form("h_pt_res_%d",ih);
    h_pt_res[ih] = new TH1F(hname,"",100,-0.01,0.01);
  }
  
  for ( unsigned int ih=0; ih<h_theta_res_tot.GetNbinsX(); ++ih ) {
    TString hname = Form("h_theta_res_%d",ih);
    h_theta_res[ih] = new TH1F(hname,"",100,-0.01,0.01) ;
  }
  

  // --- Open the root file and load the trees:

  infile = new TFile(filename);


  // Track reconstruction efficiency

  T_eff = (TTree*) infile->Get("MyClicEfficiencyCalculator/simplifiedEfficiencyTree");

  double pt_mc, phi_mc, theta_mc;
  int nhits;
  bool is_reco;
  T_eff->SetBranchAddress("m_pt", &pt_mc);
  T_eff->SetBranchAddress("m_theta", &theta_mc);
  T_eff->SetBranchAddress("m_phi", &phi_mc);
  T_eff->SetBranchAddress("m_nHits", &nhits);
  T_eff->SetBranchAddress("m_reconstructed", &is_reco);

  for(int ientry=0; ientry<T_eff->GetEntries(); ++ientry){

    T_eff->GetEntry(ientry);

    eff_phi->Fill(is_reco,phi_mc);
    eff_theta->Fill(is_reco,theta_mc);

    if (theta_mc<8. || theta_mc>172. ) continue;

    eff_pt->Fill(is_reco,pt_mc);

  }

    
  // Track parameters plots

  T_trk = (TTree*) infile->Get("MyClicEfficiencyCalculator/perfTree");
  
  std::vector<double> *pt_true=0;
  std::vector<double> *phi_true=0;
  std::vector<double> *theta_true=0;
  std::vector<double> *pt_reco=0;
  std::vector<double> *phi_reco=0;
  std::vector<double> *theta_reco=0;
  std::vector<double> *chi2_reco=0;
  std::vector<double> *nhit_reco=0;
  std::vector<double> *purity_reco=0;
  T_trk->SetBranchAddress("truePt", &pt_true);
  T_trk->SetBranchAddress("truePhi", &phi_true);
  T_trk->SetBranchAddress("trueTheta", &theta_true);
  T_trk->SetBranchAddress("recoPt", &pt_reco);
  T_trk->SetBranchAddress("recoPhi", &phi_reco);
  T_trk->SetBranchAddress("recoTheta", &theta_reco);
  T_trk->SetBranchAddress("recoChi2OverNDF", &chi2_reco);
  T_trk->SetBranchAddress("recoNhits", &nhit_reco);
  T_trk->SetBranchAddress("recoPurity", &purity_reco);

  for(unsigned int ientry=0; ientry<T_trk->GetEntries(); ++ientry){

    T_trk->GetEntry(ientry);

    for (unsigned int itrk=0; itrk<pt_reco->size(); ++itrk){

      h_pt_reco->Fill(pt_reco->at(itrk));
      h_phi_reco->Fill(phi_reco->at(itrk));
      h_theta_reco->Fill(theta_reco->at(itrk));
      h_chi2_reco->Fill(chi2_reco->at(itrk));
      h_nhit_reco->Fill(nhit_reco->at(itrk));

      // pt resolution
      
      double res_pt = (pt_reco->at(itrk)-pt_true->at(itrk))/(pt_true->at(itrk)*pt_true->at(itrk));

      int pt_bin = (int) ( pt_true->at(itrk)/h_pt_res_tot.GetBinWidth(1) );
      if ( pt_bin < 0 )
        pt_bin = 0;

      double theta = theta_true->at(itrk);
      if ( theta_true->at(itrk) > 90. )
	      theta = 180.-theta_true->at(itrk);
      int theta_bin = (int) ( theta/h_theta_res_tot.GetBinWidth(1) );

      h_pt_res[pt_bin]->Fill(res_pt);
      h_theta_res[theta_bin]->Fill(res_pt);
	
    }
      
  }

  // --- Fit the pt and theta bins and make the resolution plots:
  TF1* fitFunc = new TF1("fitFunc", "gaus");

  for ( unsigned int ih=0; ih<h_pt_res_tot.GetNbinsX(); ++ih ) {

    fitFunc->SetParameters(1.,h_pt_res[ih]->GetMean(), h_pt_res[ih]->GetRMS());
    
    h_pt_res[ih]->Fit("fitFunc","0E");

    h_pt_res_tot.SetBinContent(ih+1,fitFunc->GetParameter(2));
    h_pt_res_tot.SetBinError(ih+1,fitFunc->GetParError(2));
  
  }

  for ( unsigned int ih=0; ih<h_theta_res_tot.GetNbinsX(); ++ih ) {

    fitFunc->SetParameters(1.,h_theta_res[ih]->GetMean(), h_theta_res[ih]->GetRMS());
    
    h_theta_res[ih]->Fit("fitFunc","0E");

    h_theta_res_tot.SetBinContent(ih+1,fitFunc->GetParameter(2));
    h_theta_res_tot.SetBinError(ih+1,fitFunc->GetParError(2));
    
  }


  c[0] = new TCanvas("c_0","pT efficiency",800,600);
  eff_pt->Draw();
  c[0]->Update();
  // forces Y range to [0,1]
  eff_pt->GetPaintedGraph()->GetYaxis()->SetRangeUser(0., 1.05);
  if ( saveonfile )
    c[0]->Print("Eff_pt.png");

  c[1] = new TCanvas("c_1","phi efficiency",800,600);
  eff_phi->Draw();
  if ( saveonfile )
    c[1]->Print("Eff_phi.png");  
  
  c[2] = new TCanvas("c_2","theta efficiency",800,600);
  eff_theta->Draw();
  if ( saveonfile )
    c[2]->Print("Eff_theta.png");

  c[3] = new TCanvas("c_3","pT resolution vs pT",800,600);
  h_pt_res_tot.Draw();
  if ( saveonfile )
    c[3]->Print("Res_vs_pt.png");

  c[4] = new TCanvas("c_4","pT resolution vs theta",800,600);
  h_theta_res_tot.Draw();
  if ( saveonfile )
    c[4]->Print("Res_vs_theta.png");

  if ( debug ) {
    basevar[0] = new TCanvas("basevar_0","pT",800,600);
    h_pt_reco->Draw();
    if ( saveonfile )
      basevar[0]->Print("pt_reco.png");

    basevar[1] = new TCanvas("basevar_1","phi",800,600);
    h_phi_reco->Draw();
    if ( saveonfile )
      basevar[0]->Print("phi_reco.png");

    basevar[2] = new TCanvas("basevar_2","theta",800,600);
    h_theta_reco->Draw();
    if ( saveonfile )
      basevar[2]->Print("theta_reco.png");

    basevar[3] = new TCanvas("basevar_3","chi2",800,600);
    h_chi2_reco->Draw();
    if ( saveonfile )
      basevar[3]->Print("chi2_reco.png");

    basevar[4] = new TCanvas("basevar_4","nhit",800,600);
    h_nhit_reco->Draw();
    if ( saveonfile )
      basevar[4]->Print("nhit_reco.png");

  }

}
