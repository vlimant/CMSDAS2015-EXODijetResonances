#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include <iostream>

void CMSDASmacro() {

  TFile* file = new TFile("test/hDijetMassPseudo.root");

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  //------------------Dijet mass spectrum---------------------------

  // get pointers to some histograms we will want
  TH1D* hDijetMass=(TH1D*)file->Get("hDijetMassPseudo");
  TH1D* hDiffXsec=(TH1D*)hDijetMass->Clone("hDiffXsec");
  hDiffXsec->Reset();

  float intLumi = 100.; // integrated lumi in /pb

  // compute the differential cross-section (bin content/bin width/luminosity)
  for(int i=1; i<=hDijetMass->GetNbinsX(); ++i) 
  {
    double width=hDijetMass->GetBinWidth(i);
    double n=hDijetMass->GetBinContent(i);
    double err=TMath::Sqrt(n);

    hDiffXsec->SetBinContent(i, n/width/intLumi);
    hDiffXsec->SetBinError(i, err/width/intLumi);
  }

  // make a home for the pretty plots
  TCanvas *c = new TCanvas("c", "",1000,800);
  hDiffXsec->GetXaxis()->SetTitle("Dijet Mass [GeV]");
  hDiffXsec->GetYaxis()->SetTitle("d#sigma/dm [pb/GeV]");
  hDiffXsec->SetTitle("Differential Cross Section");
  hDiffXsec->Draw();

  c->SetLogy();
  
  c->SaveAs("hDiffXsec.png");
  
  //---------------------------------------------------------------------
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
  
  c->SaveAs("hDiffCountsWithFit.png");
  
  //---------------------------------------------------------------------
  // now make a histogram for the values of the fit
  TH1D* hResiduals = (TH1D*)hDijetMass->Clone("hResiduals");
  hResiduals->Reset();
  hResiduals->SetTitle("(Data - Fit)/Fit;Dijet mass [GeV]");

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
  }

  hResiduals->SetMinimum(-1.);
  hResiduals->SetMaximum(4.);
  hResiduals->GetXaxis()->SetRangeUser(900,6100);
  hResiduals->Draw();
  
  TLine *line = new TLine(1118.,0.,6099.,0.);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw("same");
  
  c->SetLogy(0);
  
  c->SaveAs("hResiduals.png");

  //---------------------------------------------------------------------
  // calculate the pulls
  TH1D* hPulls = (TH1D*)hDijetMass->Clone("hPulls");
  hPulls->Reset();
  hPulls->SetTitle("(Data - Fit)/Error;Dijet mass [GeV]");

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
      hPulls->SetBinContent(bin, (data_val - fit_val)/err_val );
      hPulls->SetBinError(bin, err_val/err_val );
    }
  }

  hPulls->SetMinimum(-3.);
  hPulls->SetMaximum(3.);
  hPulls->GetXaxis()->SetRangeUser(900,6100);
  hPulls->Draw();
  
  line->Draw("same");
  
  c->SetLogy(0);
  
  c->SaveAs("hPulls.png");

}
