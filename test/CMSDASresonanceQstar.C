#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include <iostream>

void CMSDASresonanceQstar() {

  TFile* fileQstar = new TFile("histos_Qstar.root");

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  // get event counts
  TH1F *eventCountQstar = (TH1F*)fileQstar->Get("dijetAna/hEventCount");
  // get dijet mass distribution
  TH1F *hCorrDijetMassQstar = (TH1F*)fileQstar->Get("dijetAna/hCorrDijetMass");

  float xsQstar = 0.8333E-01; // Qstar xs in pb
  float intLumi = 100.;       // integrated lumi in /pb

  // scale simulated Qstar dijet mass distibution to 100 /pb
  hCorrDijetMassQstar->Scale((xsQstar*intLumi)/eventCountQstar->GetBinContent(4)); // Qstar xsection includes |eta|<2.5 and |Delta eta|<1.3 cuts. Hence, we use event counts after the full event selection

  TH1D* hDiffXsecQstar=(TH1D*)hCorrDijetMassQstar->Clone("hDiffXsecQstar");
  hDiffXsecQstar->Reset();
  
  // compute the differential cross-section (bin content/bin width/luminosity)
  for(int i=1; i<=hCorrDijetMassQstar->GetNbinsX(); ++i) 
  {
    double width=hCorrDijetMassQstar->GetBinWidth(i);
    double n=hCorrDijetMassQstar->GetBinContent(i);
    double err=hCorrDijetMassQstar->GetBinError(i);

    hDiffXsecQstar->SetBinContent(i, n/width/intLumi);
    hDiffXsecQstar->SetBinError(i, err/width/intLumi);
  }

  // make a home for the pretty plots
  TCanvas *c = new TCanvas("c", "",1000,800);
  hDiffXsecQstar->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  hDiffXsecQstar->GetYaxis()->SetTitle("d#sigma/dm [pb/GeV]");
  hDiffXsecQstar->SetTitleOffset(1.4,"Y");
  hDiffXsecQstar->SetTitle("Differential Cross Section for m_{Qstar}=4 TeV");
  hDiffXsecQstar->SetLineWidth(2);
  hDiffXsecQstar->Draw("hist");
  
  c->SaveAs("hDiffXsecQstar.png");
  
  //---------------------------------------------------------------------
  TFile* file = new TFile("test/hDijetMassPseudo.root");
  
  TH1D* hDijetMass=(TH1D*)file->Get("hDijetMassPseudo");
  TH1D* hDiffCounts=(TH1D*)hDijetMass->Clone("hDiffCounts");
  hDiffCounts->Reset();

  // compute the differential counts (bin content/bin width)
  for(int i=1; i<=hDijetMass->GetNbinsX(); ++i) 
  {
    double width=hDijetMass->GetBinWidth(i);
    double n=hDijetMass->GetBinContent(i);
    double err=TMath::Sqrt(n);

    hDiffCounts->SetBinContent(i, n/width);
    hDiffCounts->SetBinError(i, err/width);
  }
  
  hDiffCounts->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  hDiffCounts->GetYaxis()->SetTitle("dN/dm [1/GeV]");
  hDiffCounts->SetTitle("Differential Counts");
  hDiffCounts->SetMinimum(1e-4);

  // declare a function to fit, over some range, and fit it
  TF1 *fit = new TF1("fit", "[0]*pow(1-x/13000.0,[1])/pow(x/13000.,[2]+[3]*log(x/13000.))",1118,6099);
  fit->SetParameter(0,2.94554e-04);
  fit->SetParameter(1,5.57678e+00);
  fit->SetParameter(2,7.09806e+00);
  fit->SetParameter(3,3.48634e-01);
  
  hDiffCounts->Fit("fit","RLI");
  
  hDiffCounts->Draw();


  TH1D* hDiffCountsQstar=(TH1D*)hCorrDijetMassQstar->Clone("hDiffCountsQstar");
  hDiffCountsQstar->Reset();
  
  // compute the differential counts (bin content/bin width)
  for(int i=1; i<=hCorrDijetMassQstar->GetNbinsX(); ++i) 
  {
    double width=hCorrDijetMassQstar->GetBinWidth(i);
    double n=hCorrDijetMassQstar->GetBinContent(i);
    double err=hCorrDijetMassQstar->GetBinError(i);

    hDiffCountsQstar->SetBinContent(i, n/width);
    hDiffCountsQstar->SetBinError(i, err/width);
  }
  
  hDiffCountsQstar->SetLineWidth(2);
  hDiffCountsQstar->Draw("histsame");
  
  c->SetLogy();

  c->SaveAs("hDiffCountsWithFit_Qstar.png");
  
  //---------------------------------------------------------------------
  TH1D* hResiduals = (TH1D*)hDijetMass->Clone("hResiduals");
  hResiduals->Reset();
  hResiduals->SetTitle("(Data - Fit)/Fit with contribution from 4 TeV Qstar;Dijet mass [GeV]");
  
  TH1D* hResidualsQstar=(TH1D*)hCorrDijetMassQstar->Clone("hResidualsQstar");
  hResidualsQstar->Reset();

  // fill the histogram of the data minus the fit integral
  for(int bin=1; bin<=hDijetMass->GetNbinsX(); ++bin) 
  {
    double data_val = hDijetMass->GetBinContent(bin);
    double err_val  = TMath::Sqrt(data_val);
    double m_low    = hDijetMass->GetBinLowEdge(bin);
    double m_high   = m_low + hDijetMass->GetBinWidth(bin);
    double fit_val  = fit->Integral(m_low,m_high);
    // skip bins with no data value
    if (data_val != 0.0) {
      hResiduals->SetBinContent(bin, (data_val - fit_val)/fit_val );
      hResiduals->SetBinError(bin, err_val/fit_val );
    }

    double res_val = hCorrDijetMassQstar->GetBinContent(bin);
    if (res_val != 0.0) {
      hResidualsQstar->SetBinContent(bin, res_val/fit_val );
    }
  }

  hResiduals->SetMinimum(-1.);
  hResiduals->SetMaximum(4.);
  hResiduals->GetXaxis()->SetRangeUser(900,6100);
  hResiduals->Draw();

  hResidualsQstar->SetLineWidth(2);
  hResidualsQstar->Draw("histsame");
  
  TLine *line = new TLine(1118.,0.,6099.,0.);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw("same");
  
  c->SetLogy(0);
  
  c->SaveAs("hResiduals_Qstar.png");
  
  //---------------------------------------------------------------------
  TFile *fileQCD = new TFile("test/histos_FullStats.root");
  
  TH1F *eventCountQCD = (TH1F*)fileQCD->Get("dijetAna/hEventCount");
  TH1F *hCorrDijetMassQCD = (TH1F*)fileQCD->Get("dijetAna/hCorrDijetMass");

  float xsQCD = 2.435e+09; // total xs in pb

  // scale simulated QCD dijet mass distibution to 100 /pb
  hCorrDijetMassQCD->Scale((xsQCD*intLumi)/eventCountQCD->GetBinContent(2)); // QCD xsection does not include |eta|<2.5 and |Delta eta|<1.3 cuts. Hence, we use total event counts

  TH1D* hDiffCountsQCD=(TH1D*)hCorrDijetMassQstar->Clone("hDiffCountsQCD");
  hDiffCountsQCD->Reset();
  
  // compute the differential counts (bin content/bin width)
  for(int i=1; i<=hCorrDijetMassQCD->GetNbinsX(); ++i) 
  {
    double width=hCorrDijetMassQCD->GetBinWidth(i);
    double n=hCorrDijetMassQCD->GetBinContent(i);
    double err=hCorrDijetMassQCD->GetBinError(i);

    hDiffCountsQCD->SetBinContent(i, n/width);
    hDiffCountsQCD->SetBinError(i, err/width);
  }
  
  hDiffCountsQCD->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  hDiffCountsQCD->GetYaxis()->SetTitle("dN/dm [1/GeV]");
  hDiffCountsQCD->SetTitle("Differential Counts");
  hDiffCountsQCD->SetFillColor(42);
  hDiffCountsQCD->Draw("hist");

  hDiffCounts->Draw("same");
  
  c->SetLogy();

  c->SaveAs("hDiffCountsWithFit_QCD.png");
}