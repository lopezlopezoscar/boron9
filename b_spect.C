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

// -----------------------    Focal Plane position
//------------------------    Spectrum + Tritons Cut only
//TCanvas *cFocalPlane1pos = new TCanvas();
 dchain->SetAlias("newX1pos","X1pos+0.0127483*(tof-3120)-0.000332362*pow(tof-3120,2.)"); // Scattering Angle correction Daniel manual
 TH1F *hnewX1pos = new TH1F("hnewX1pos","X1pos+Scatter Angle Correction 9B;Focal Plane (arb. units);Counts",1200,-10,800);
 dchain->Draw("newX1pos>>hnewX1pos","CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag","col");
 hnewX1pos->Write();

// p0 = 4.913, p1 = -0.0106291
//dchain->SetAlias("Exx","7.848 - 0.0106291*newX1pos");//Slope from 26Mg calib - intercept used to fix G.S. = 0
//dchain->SetAlias("Exx","7.938 - newX1pos*0.0111842 + 0.000000584097*pow(newX1pos,2.)");
dchain->SetAlias("Exx","50.0-40.9632-1.08662683205-0.01215 - newX1pos*0.0111842 + 0.000000584097*pow(newX1pos,2.)");
//50.0-40.9662-1.08662683205-0.01215

// --------------------------------------------------
// -----------------------    Singles ---------------
// --------------------------------------------------

// -----------------------    Plot CAKE vs Exx
 //TCanvas *cCakeExx1 = new TCanvas();
  TH2F *SiliconEnergysvsExx1 = new TH2F("SiliconEnergysvsExx1","SiliconEnergy_:Exx_ 9B;Ex [MeV];Silicon Energy [arb. units] ",1200,-1,9,1200,0,10000);
  dchain->Draw("SiliconEnergy:Exx>>SiliconEnergysvsExx1","CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag","col");
  SiliconEnergysvsExx1->Write();

// -----------------------    Excitation Energy
// -----------------------    Singles Spectrum  // Singles = FP ONLY (NO CAKE)
 //TCanvas *cSingles = new TCanvas();
 TH1F *Ex_sinlges = new TH1F("Ex_sinlges","9B Singles Spectrum ; Ex (MeV);Counts",1200,-1,8);
 dchain->Draw("Exx>>Ex_sinlges","CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag","col"); 
   Ex_sinlges->Write();

// -----------------------    Plot Exx using CUTS from Plot CAKE vs Exx
//  ----------------------    Spectrum 9B + Protons CUT + low lying region
  //TCanvas *cpCut1 = new TCanvas();
  TH1F *Spectrum_pCut1 = new TH1F("Spectrum_pCut1","9B Spectrum proton Cut ; Ex (MeV);Counts",1200,-2,8);
  dchain->Draw("Exx>>Spectrum_pCut1","CUTpad1tof && CUTprotons_gs && !X1flag && !U1flag","col"); //cut9BProtonsThin
  Spectrum_pCut1->Write();

//  ------ Spectrum 9B + 9B gs CUT + low lying region
  TCanvas *c9BgsCut1 = new TCanvas();
  TH1F *Spectrum_9BgsCut1 = new TH1F("Spectrum_9BgsCut1","9B Spectrum 9B g.s. Cut ; Ex (MeV);Counts",1200,-2,8);
  dchain->Draw("Exx>>Spectrum_9BgsCut1","CUTpad1tof && CUT9B_gs && !X1flag && !U1flag","col");
  Spectrum_9BgsCut1->Write();
   
//  ------ Spectrum 9B + 5Li CUT + low lying region
  //TCanvas *cAlphaCut1 = new TCanvas();
  TH1F *Spectrum_5LiCut1 = new TH1F("Spectrum_5LiCut1","9B Spectrum alpha Cut ; Ex (MeV);Counts",1200,-2,8);
  dchain->Draw("Exx>>Spectrum_5LiCut1","CUTpad1tof && CUT5Li && !X1flag && !U1flag","col");
  Spectrum_5LiCut1->Write();

// --------------------------------------------------
// -----------------------    Coincidences ----------
// --------------------------------------------------

// -----------------------    Plot CAKE vs Exx
// -----------------------    For make a TCutG of SiliconEnergy vs Exx " CUTprotons.C & CUTalfas.C"
 //TCanvas *cCakeExx2 = new TCanvas();
  TH2F *SiliconEnergysvsExx2 = new TH2F("SiliconEnergysvsExx2","SiliconEnergy_:Exx_ 9B;Ex [MeV];Silicon Energy [arb. units] ",1200,-1,9,1200,0,10000);
  dchain->Draw("SiliconEnergy:Exx>>SiliconEnergysvsExx2","CUTpad1tof && CUTpad1X1 && !X1flag && !U1flag && abs(SiliconTime-tof)<200 && SiliconHits==1 &&  TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720","col");
  SiliconEnergysvsExx2->Write();

// -----------------------    Excitation Energy
// -----------------------    Coincidences  ---  Spectrum + abs(SiliconTime-tof)<200  
 //TCanvas *cCoincidences = new TCanvas();
 // p0 = 4.913, p1 = -0.0106291
 TH1F *Ex_coincidences = new TH1F("Ex_coincidences","a) 9B Coincidences Spectrum ; Ex (MeV);Counts",1200,-1,8);
 dchain->Draw("Exx>>Ex_coincidences","CUTpad1tof && CUTpad1X1 && abs(SiliconTime-tof)<200 && !X1flag && !U1flag && SiliconHits==1 &&  TDCChannelFront!=868 && TDCChannelFront!=870 && TDCChannelFront!=880 && TDCChannelFront!=720","col");
   Ex_coincidences->Write();
  
// -----------------------    Plot Exx using CUTS from Plot CAKE vs Exx
//  ----------------------    Spectrum 9B + Protons CUT + low lying region
  //TCanvas *cpCut_2 = new TCanvas();
  TH1F *Spectrum_pCut2 = new TH1F("Spectrum_pCut2","9B Spectrum proton Cut ; Ex (MeV);Counts",1200,-2,8);
  dchain->Draw("Exx>>Spectrum_pCut2","CUTpad1tof && CUTprotons_gs && !X1flag && !U1flag && abs(SiliconTime-tof)<200","col"); //cut9BProtonsThin
  Spectrum_pCut2->Write();
  
//  ------ Spectrum 9B + 9B gs CUT + low lying region
  //TCanvas *c9BgsCut2 = new TCanvas();
  TH1F *Spectrum_9BgsCut2 = new TH1F("Spectrum_9BgsCut2","9B Spectrum 9B g.s. Cut ; Ex (MeV);Counts",1200,-2,8);
  dchain->Draw("Exx>>Spectrum_9BgsCut2","CUTpad1tof && CUT9B_gs && !X1flag && !U1flag && abs(SiliconTime-tof)<200","col"); //cut9BgsThin
  Spectrum_9BgsCut2->Write();
   
//  ------ Spectrum 9B + 5Li CUT + low lying region
  //TCanvas *cAlphaCut2 = new TCanvas();
  TH1F *Spectrum_5LiCut2 = new TH1F("Spectrum_5LiCut2","9B Spectrum alpha Cut ; Ex (MeV);Counts",1200,-2,8);
  dchain->Draw("Exx>>Spectrum_5LiCut2","CUTpad1tof && CUT5Li && !X1flag && !U1flag && abs(SiliconTime-tof)<200","col"); //cut5LiThin
  Spectrum_5LiCut2->Write();
}
