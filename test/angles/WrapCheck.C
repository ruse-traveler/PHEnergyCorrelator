#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TString.h>



// function to wrap an angle
double Wrap(double phi) {

  if (phi < 0.)             phi += TMath::TwoPi();
  if (phi > TMath::TwoPi()) phi -= TMath::TwoPi();
  return phi;

}



// function to set a color/marker
void Style(int color, int marker, TH1D* hist) {

  hist -> SetLineColor(color);
  hist -> SetMarkerColor(color);
  hist -> SetMarkerStyle(marker);
  return;

}



// mc test
void WrapCheck(const std::size_t nmc = 100000, const TString out = "wraptest.root") {

  // for getting phi spin
  TRandom3* rando = new TRandom3();

  // sine function to generate input distribution
  TF1* fSinIn  = new TF1("fSinIn",  "sin(x)+1.",  -1.*TMath::Pi(), TMath::Pi());
  TF1* fSinIn2 = new TF1("fSinIn2", "sin(2.*x)+1.", -1.*TMath::Pi(), TMath::Pi());

  // turn on histogram errors
  TH1::SetDefaultSumw2(true);

  // unwrapped histograms
  TH1D* hSinIn   = new TH1D("hSinIn",   "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hSinIn2  = new TH1D("hSinIn2",  "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiSpi  = new TH1D("hPhiSpi",  "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiSpi2 = new TH1D("hPhiSpi2", "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiHad  = new TH1D("hPhiHad",  "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiHad2 = new TH1D("hPhiHad2", "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiDif  = new TH1D("hPhiDif",  "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiDif2 = new TH1D("hPhiDif2", "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());

  // unwrapped histograms
  TH1D* hSinInW   = new TH1D("hSinInWrap",   "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hSinIn2W  = new TH1D("hSinIn2Wrap",  "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiSpiW  = new TH1D("hPhiSpiWrap",  "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiSpi2W = new TH1D("hPhiSpi2Wrap", "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiHadW  = new TH1D("hPhiHadWrap",  "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiHad2W = new TH1D("hPHiHad2Wrap", "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiDifW  = new TH1D("hPhiDifWrap",  "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());
  TH1D* hPhiDif2W = new TH1D("hPhiDif2Wrap", "", 360, -2.*TMath::TwoPi(), 2.*TMath::TwoPi());

  // mc loop
  for (std::size_t imc = 0; imc < nmc; ++imc) {

    // sine(x) test -----------------------------------------------------------

    // randomly sample from sine, uniform spin
    //   --> goal is to constrain spin - (2*had) to a sine function
    double sin  = fSinIn -> GetRandom(-1.*TMath::Pi(), 1.*TMath::Pi());
    double spin = rando  -> Uniform(-1.*TMath::Pi(), 1.*TMath::Pi());
    double had  = spin - sin;

    // now calculate unwrapped difference
    double dif  = spin - had;

    // fill unwrapped histograms
    hSinIn  -> Fill(sin);
    hPhiSpi -> Fill(spin);
    hPhiHad -> Fill(had);
    hPhiDif -> Fill(dif);

    // now wrap angles
    double sinw  = Wrap(sin);
    double spinw = Wrap(spin);
    double hadw  = Wrap(had);
    double difw  = Wrap(spinw - hadw);

    // fill wrapped histograms
    hSinInW  -> Fill(sinw);
    hPhiSpiW -> Fill(spinw);
    hPhiHadW -> Fill(hadw);
    hPhiDifW -> Fill(difw);

    // sine(2x) test ----------------------------------------------------------

    // randomly sample from sine(2x), uniform spin
    //   --> goal is to constrain spin - had to a sine(2x) function
    double sin2  = fSinIn2 -> GetRandom(-1.*TMath::Pi(), 1.*TMath::Pi());
    double spin2 = rando   -> Uniform(-1.*TMath::Pi(), 1.*TMath::Pi());
    double had2  = (spin2 - sin2) / 2.;

    // now calculate unwrapped difference
    double dif2 = spin2 - (2.*had2);

    // fill unwrapped histograms
    hSinIn2  -> Fill(sin2);
    hPhiSpi2 -> Fill(spin2);
    hPhiHad2 -> Fill(had2);
    hPhiDif2 -> Fill(dif2);

    // now wrap angles
    double sin2w  = Wrap(sin2);
    double spin2w = Wrap(spin2);
    double had2w  = Wrap(had2);
    double dif2w  = Wrap(spin2w - (2.*had2w)); 

    // fill wrapped histograms
    hSinIn2W  -> Fill(sin2w);
    hPhiSpi2W -> Fill(spin2w);
    hPhiHad2W -> Fill(had2w);
    hPhiDif2W -> Fill(dif2w);

  }

  // open output file
  TFile *fout = new TFile(out, "recreate");

  // set line / marker colors
  Style(923, 29, hSinIn);
  Style(923, 29, hSinInW);
  Style(923, 29, hSinIn2);
  Style(923, 29, hSinIn2W);

  Style(799, 24, hPhiSpi);
  Style(799, 24, hPhiSpiW);
  Style(799, 24, hPhiSpi2);
  Style(799, 24, hPhiSpi2W);

  Style(899, 25, hPhiHad);
  Style(899, 25, hPhiHadW);
  Style(899, 25, hPhiHad2);
  Style(899, 25, hPhiHad2W);

  Style(879, 26, hPhiDif);
  Style(879, 26, hPhiDifW);
  Style(879, 26, hPhiDif2);
  Style(879, 26, hPhiDif2W);

  // make legends
  TLegend* leg = new TLegend(0.1, 0.1, 0.3, 0.3);
  leg -> SetFillColor(0);
  leg -> SetLineColor(0);
  leg -> SetTextFont(42);
  leg -> SetTextAlign(12);
  leg -> AddEntry(hSinIn, "Input sin(x)");
  leg -> AddEntry(hPhiSpi, "#varphi_{s}");
  leg -> AddEntry(hPhiHad, "#varphi_{h}, sine(x) - #varphi_{s}");
  leg -> AddEntry(hPhiDif, "#Delta#varphi = #varphi_{s} - #varphi_{h}");

  TLegend* leg2 = new TLegend(0.1, 0.1, 0.3, 0.3);
  leg2 -> SetFillColor(0);
  leg2 -> SetLineColor(0);
  leg2 -> SetTextFont(42);
  leg2 -> SetTextAlign(12);
  leg2 -> AddEntry(hSinIn2, "Input sin(2x)");
  leg2 -> AddEntry(hPhiSpi2, "#varphi_{s}");
  leg2 -> AddEntry(hPhiHad2, "#varphi_{h}, (sine(2x) - #varphi_{s})/2");
  leg2 -> AddEntry(hPhiDif2, "#Delta#varphi = #varphi_{s} - 2#varphi_{h}");

  // make plots
  TCanvas* cSine = new TCanvas("cSine", "", 750, 750);
  cSine   -> cd();
  hSinIn  -> Draw();
  hPhiSpi -> Draw("same");
  hPhiHad -> Draw("same");
  hPhiDif -> Draw("same");
  leg     -> Draw();
  fout    -> cd();
  cSine   -> Write();
  cSine   -> Close();

  TCanvas* cSineW = new TCanvas("cSineWrap", "", 750, 750);
  cSineW   -> cd();
  hSinInW  -> Draw();
  hPhiSpiW -> Draw("same");
  hPhiHadW -> Draw("same");
  hPhiDifW -> Draw("same");
  leg      -> Draw();
  fout     -> cd();
  cSineW   -> Write();
  cSineW   -> Close();

  TCanvas* cSine2 = new TCanvas("cSine2", "", 750, 750);
  cSine2   -> cd();
  hSinIn2  -> Draw();
  hPhiSpi2 -> Draw("same");
  hPhiHad2 -> Draw("same");
  hPhiDif2 -> Draw("same");
  leg2     -> Draw();
  fout     -> cd();
  cSine2   -> Write();
  cSine2   -> Close();

  TCanvas* cSine2W = new TCanvas("cSine2Wrap", "", 750, 750);
  cSine2W   -> cd();
  hSinIn2W  -> Draw();
  hPhiSpi2W -> Draw("same");
  hPhiHad2W -> Draw("same");
  hPhiDif2W -> Draw("same");
  leg2      -> Draw();
  fout      -> cd();
  cSine2W   -> Write();
  cSine2W   -> Close();

  // save histograms
  fout      -> cd();
  fSinIn    -> Write();
  fSinIn2   -> Write();
  hSinIn    -> Write();
  hSinIn2   -> Write();
  hPhiSpi   -> Write();
  hPhiSpi2  -> Write();
  hPhiHad   -> Write();
  hPhiHad2  -> Write();
  hPhiDif   -> Write();
  hPhiDif2  -> Write();
  hSinInW   -> Write();
  hSinIn2W  -> Write();
  hPhiSpiW  -> Write();
  hPhiSpi2W -> Write();
  hPhiHadW  -> Write();
  hPhiHad2W -> Write();
  hPhiDifW  -> Write();
  hPhiDif2W -> Write();

  // close output and exit
  fout -> cd();
  fout -> Close();
  return;

}
