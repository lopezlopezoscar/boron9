#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TCutG.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TText.h"
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TAxis.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <TROOT.h>
#include <TRint.h>

using namespace std;
// ----- for run fron Terminal $root -b -l b_spect.C ----- -b supress GUI 
void b_spect()
{
   TStopwatch t;

   TFile *output = new TFile("b_spect_siliconcut.root","recreate");

TString Treename1 = "DATA";
TChain * dchain = new TChain(Treename1);

dchain->Add("sorted10043.root"); // Target 2 (thin)
dchain->Add("sorted10045.root");
dchain->Add("sorted10047.root"); 
dchain->Add("sorted10048.root"); 
dchain->Add("sorted10050.root"); 
dchain->Add("sorted10051.root"); 
dchain->Add("sorted10053.root"); 
dchain->Add("sorted10054.root"); 
dchain->Add("sorted10056.root"); 
dchain->Add("sorted10057.root"); 
dchain->Add("sorted10059.root"); 
dchain->Add("sorted10060.root"); 
dchain->Add("sorted10062.root"); 
dchain->Add("sorted10063.root"); 
dchain->Add("sorted10065.root"); 
dchain->Add("sorted10067.root"); 
dchain->Add("sorted10068.root"); 
dchain->Add("sorted10070.root"); 
dchain->Add("sorted10071.root"); 
dchain->Add("sorted10073.root"); 

int entries = dchain->GetEntries();

for(int j = 0; j < entries; j++)
	{
		dchain->GetEntry(j);
	}

// ---- CUT  ----
  gROOT->ProcessLine(".x CUTpad1tof.C");
  gROOT->ProcessLine(".x CUTpad1X1.C");
  gROOT->ProcessLine(".x CUTptime.C");
  gROOT->ProcessLine(".x CUTatime.C");
  gROOT->ProcessLine(".x CUTprotons.C");
  gROOT->ProcessLine(".x CUTalfas.C"); 
  gROOT->ProcessLine(".x CUTprotons_gs.C");
  gROOT->ProcessLine(".x CUT9B_gs.C");
  gROOT->ProcessLine(".x CUT5Li.C");

// -----------------------    Plot Pad Energy vs TOF for PID
// -----------------------    For make a TCutG of tritons " CUTpad1tof.C "
//TCanvas *cCUTpad1tof = new TCanvas();
 TH2F *hpad1vstof = new TH2F("hpad1vstof","pad1:tof ",600,2800,3400,2500,0,2500);
 dchain->Draw("pad1:tof>>hpad1vstof","!X1flag && !U1flag ","col");
 hpad1vstof->Write();

// -----------------------    Plot Pad Energy vs FP position
// -----------------------    For make a TCutG of pad1 vs X1pos " CUTpad1X1.C "
//TCanvas *cCUTpad1X1 = new TCanvas();
 TH2F *hCUTpad1X1 = new TH2F("hCUTpad1X1","pad1:X1pos ",1200,0,900,2500,0,2500);
 dchain->Draw("pad1:X1pos>>hCUTpad1X1","!X1flag && !U1flag","col");
 hCUTpad1X1->Write();

dchain->SetAlias("newX1pos","X1pos+0.0127483*(tof-3120)-0.000332362*pow(tof-3120,2.)"); // Scattering Angle correction Daniel manual

/*
//--------------TOF v.s. Focal Plane position (bananas)
//tof versus Focal plane position X1pos
TCanvas *c3 = new TCanvas();
   TH2F *htofvsX1pos = new TH2F("htofvsX1pos","TOF vs X1pos",1000,-100,1000,1000,2950,3300);
   dchain->Draw("tof:X1pos>>htofvsX1pos","CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag","col");
   htofvsX1pos->Write();

//tof versus Focal plane position newX1pos
TCanvas *c4 = new TCanvas();
   TH2F *htofvsnewX1pos = new TH2F("htofvsnewX1pos","TOF vs newX1pos",1000,-100,1000,1000,2950,3300);
   dchain->Draw("tof:newX1pos>>htofvsnewX1pos","CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag","col");
   htofvsnewX1pos->Write();
*/

//------- Spectrum + Tritons Cut only
//TCanvas *cFocalPlane1pos = new TCanvas();
 TH1F *hnewX1pos = new TH1F("hnewX1pos","X1pos+Scatter Angle Correction 9B;Focal Plane (arb. units);Counts",1200,-10,800);
 dchain->Draw("newX1pos>>hnewX1pos","CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag","col");
 hnewX1pos->Write();

// ------- Spectrum + Tritons Cut in Kinetic Energy of Tritons
// We take values from Terminal for linear calibration
// p0 = 40.9662, p1 = 0.0111842, p2 = -8.32757e-07 
 //dchain->SetAlias("Et","40.9662+newX1pos*0.0111842-0.000000832757*pow(newX1pos,2.)");
 dchain->SetAlias("Et","40.9662+newX1pos*0.0111842-0.000000584097*pow(newX1pos,2.)"); //(Typo???)
//TCanvas *cFocalPlane2Kin = new TCanvas();
 TH1F *K = new TH1F("K","Kinetic_t + ScatCorr 9B; Kinetic_t (MeV);Counts",1200,40.85,49.38);
 dchain->Draw("Et>>K","CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag","col");
 K->Write();

// ------- --- Spectrum + Tritons Cut + Calibrated in Excitation Energy
// Ex = Eb - Et - Qval - Eloss
  dchain->SetAlias("Exx","50.0-Et-1.08662683205-0.01215"); //
//TCanvas *cFocalPlane3Ex = new TCanvas();
  TH1F *Ex_t = new TH1F("Ex_t","9B Spectrum Calibrated ; Ex (MeV);Counts",1200,-1,8);
  dchain->Draw("Exx>>Ex_t","CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag","col");
  Ex_t->Write();

// ---------* a) Single Spectrum  ---  Spectrum + abs(SiliconTime-tof)<200  
 //TCanvas *cSingles = new TCanvas();
 TH1F *Ex_sinlges = new TH1F("Ex_sinlges","a) 9B Spectrum (Target 2) Singles ; Ex (MeV);Counts",1200,-1,8);
 dchain->Draw("Exx>>Ex_sinlges","CUTpad1tof && CUTpad1X1 && abs(SiliconTime-tof)<200 && !X1flag && !U1flag && SiliconHits==1 &&  TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720","col");
   Ex_sinlges->Write();

// -----------------------    Plot CAKE vs Exx
// -----------------------    For make a TCutG of SiliconEnergy vs Exx " CUTprotons.C & CUTalfas.C"
 //TCanvas *cCakeExx = new TCanvas();
  TH2F *SiliconEnergysvsExx = new TH2F("SiliconEnergysvsExx","SiliconEnergy_:Exx_ 9B;Ex [MeV];Silicon Energy [arb. units] ",1200,-1,9,1200,0,10000);
  dchain->Draw("SiliconEnergy:Exx>>SiliconEnergysvsExx","CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag && abs(SiliconTime-tof)<200 && SiliconHits==1 &&  TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720","col");
  // omitting time gatting
  //dchain->Draw("SiliconEnergy:Exx>>SiliconEnergysvsExx","CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag && SiliconHits==1 &&  TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720","col");
  SiliconEnergysvsExx->Write();


// Plot Exx using CUTS from Plot CAKE vs Exx
//  ------ Spectrum 9B + Protons CUT + low lying region
  //TCanvas *cSpectrum = new TCanvas();
  TH1F *Spectrum_pCut = new TH1F("Spectrum_pCut","9B Spectrum_pCut ; Ex (MeV);Counts",1200,-2,8);
  dchain->Draw("Exx>>Spectrum_pCut","CUTpad1tof && CUTprotons_gs && !X1flag && !U1flag && abs(SiliconTime-tof)<200","col"); //cut9BProtonsThin
  // omitting time gatting
  //dchain->Draw("Exx>>Spectrum_pCut","CUTpad1tof && CUTprotons_gs && !X1flag && !U1flag","col");
  Spectrum_pCut->Write();
  
//  ------ Spectrum 9B + 9B gs CUT + low lying region
  //TCanvas *cSpectrum2 = new TCanvas();
  TH1F *Spectrum_9BgsCut = new TH1F("Spectrum_9BgsCut","9B Spectrum_9BgsCut ; Ex (MeV);Counts",1200,-2,8);
  //dchain->Draw("Exx>>Spectrum_9BgsCut","CUTpad1tof && CUT9B_gs && !X1flag && !U1flag && abs(SiliconTime-tof)<200","col"); //cut9BgsThin
  // omitting time gatting
  dchain->Draw("Exx>>Spectrum_9BgsCut","CUTpad1tof && CUT9B_gs && !X1flag && !U1flag","col");
  Spectrum_9BgsCut->Write();
   
//  ------ Spectrum 9B + 5Li CUT + low lying region
  //TCanvas *cSpectrum3 = new TCanvas();
  TH1F *Spectrum_5LiCut = new TH1F("Spectrum_5LiCut","9B Spectrum_5LiCut ; Ex (MeV);Counts",1200,-2,8);
  //dchain->Draw("Exx>>Spectrum_5LiCut","CUTpad1tof && CUT5Li && !X1flag && !U1flag && abs(SiliconTime-tof)<200","col"); //cut5LiThin
  // omitting time gatting
  dchain->Draw("Exx>>Spectrum_5LiCut","CUTpad1tof && CUT5Li && !X1flag && !U1flag","col");
  Spectrum_5LiCut->Write();

// -------------------  Plot Caketime vs Cake Energy
// -------------------  For make a TCutG of tritons SiliconTimeOffset-tof vs SiliconEnergy" CUTptime.C & CUTatime"
  //TCanvas *cCaketime = new TCanvas();
  TH2F *Caketime = new TH2F("Caketime","All Detectors ;ESilicon [keV];Tsi-Tk600/ns ",1200,0,10000,1200,-3500,-2900);
  dchain->Draw("SiliconTimeOffset-tof:SiliconEnergy>>Caketime"," CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag && SiliconHits==1 &&  TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720 && abs(tof)<1e4 && abs(SiliconTime[0])<1e5","col");   
  Caketime->Write();

// Plot Exx using CUTS from Plot Caketime vs Cake Energy
// ---------* b) Proton Spectrum  ---  Spectrum + abs(SiliconTime-tof)<200  
  //TCanvas *cproton = new TCanvas();
  TH1F *Ex_protons = new TH1F("Ex_protons","b) 9B Spectrum (Target 2) Proton decay gates ; Ex (MeV);Counts",1200,-1,8);
  dchain->Draw("Exx>>Ex_protons","CUTpad1tof && CUTpad1X1 && CUTptime && abs(SiliconTime-tof)<200 && !X1flag && !U1flag && SiliconHits==1 &&  TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720","col");
  Ex_protons->Write();

// ---------* b) Alphas Spectrum  ---  Spectrum + abs(SiliconTime-tof)<200  
  //TCanvas *calfas = new TCanvas();
  TH1F *Ex_alfas = new TH1F("Ex_alfas","c) 9B Spectrum (Target 2) Alpha decay gates ; Ex (MeV);Counts",1200,-1,8);
  dchain->Draw("Exx>>Ex_alfas","CUTpad1tof && CUTpad1X1 && CUTatime && abs(SiliconTime-tof)<200 && !X1flag && !U1flag && SiliconHits==1 &&  TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720","col");
  Ex_alfas->Write();

// ----------- CAKE -----------

/*
//TCanvas *c77711 = new TCanvas();
TH2F *Caketime1 = new TH2F("Caketime1","Detector 1;ESilicon [keV];Tsi-Tk600/ns ",1200,0,10000,1200,-3500,-2900);
dchain->Draw("(SiliconTimeOffset-tof):SiliconEnergy>>Caketime1"," CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag && SiliconHits==1 && TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720 && abs(tof)<1e4 && abs(SiliconTime[0])<1e5 && DetectorHit==1","col");   
Caketime1->Write();

//TCanvas *c77722 = new TCanvas();
TH2F *Caketime2 = new TH2F("Caketime2","Detector 2;ESilicon [keV];Tsi-Tk600/ns ",1200,0,10000,1200,-3500,-2900);
dchain->Draw("(SiliconTimeOffset-tof):SiliconEnergy>>Caketime2"," CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag && SiliconHits==1 && TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720 && abs(tof)<1e4 && abs(SiliconTime[0])<1e5 && DetectorHit==2","col");   
Caketime2->Write();

//TCanvas *c77733 = new TCanvas();
TH2F *Caketime3 = new TH2F("Caketime3","Detector 3;ESilicon [keV];Tsi-Tk600/ns ",1200,0,10000,1200,-3500,-2900);
dchain->Draw("(SiliconTimeOffset-tof):SiliconEnergy>>Caketime3"," CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag && SiliconHits==1 && TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720 && abs(tof)<1e4 && abs(SiliconTime[0])<1e5 && DetectorHit==3","col");   
Caketime3->Write();

//TCanvas *c77744 = new TCanvas();
TH2F *Caketime4 = new TH2F("Caketime4","Detector 4;ESilicon [keV];Tsi-Tk600/ns ",1200,0,10000,1200,-3500,-2900);
dchain->Draw("(SiliconTimeOffset-tof):SiliconEnergy>>Caketime4"," CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag && SiliconHits==1 && TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720 && abs(tof)<1e4 && abs(SiliconTime[0])<1e5 && DetectorHit==4","col");   
Caketime4->Write();

//TCanvas *c77755 = new TCanvas();
TH2F *Caketime5 = new TH2F("Caketime5","Detector 5;ESilicon [keV];Tsi-Tk600/ns ",1200,0,10000,1200,-3500,-2900);
  dchain->Draw("(SiliconTimeOffset-tof):SiliconEnergy>>Caketime5"," CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag && SiliconHits==1 && TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720 && abs(tof)<1e4 && abs(SiliconTime[0])<1e5 && DetectorHit==5","col"); 
  Caketime5->Write();
*/

  dchain->Draw("(5+StripFront)/30.*pi:(StripBack+(DetectorHit-1)*8)*2*pi/40.>>hcake(80,-6.28318,6.28318,60,-3.141592,3.141592)","","polcolz");

t.Print();
}
