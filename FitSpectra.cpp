#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TROOT.h>
#include "function_library.h"
#include <TCanvas.h>


bool VerboseFlag = false; //print out problems

void FitSpectra();

double Sp_9B = -0.1859; //MeV

double reduced_width_9B_resonance = 0.54/penetrability(4,1,8,1,1,-Sp_9B,1.25*(pow(8.,1./3.) + 1)) / 1e3; ///1e3 to convert to MeV - 0.54 keV ground-state width

double FitFunction(double *x, double *pars)
{
    double result = 0;
    
    //pars[0] = height of the 9B ground state
    //pars[1] = reduced width of the 9B ground state
    //pars[2] = energy of the 9B ground-state resonance
    //pars[3] =?
    
    if(x[0] - Sp_9B > 0)
    {
        
        double width = pars[1] * 2 * penetrability(4,1,8,1,1,x[0] - Sp_9B,1.25*(pow(8.,1./3.) + 1));
        
        if(VerboseFlag)std::cout << "Energy = " << x[0] - Sp_9B << std::endl;
        if(VerboseFlag)std::cout << "penetrability = " << penetrability(4,1,8,1,1,x[0] - Sp_9B,1.25*(pow(8.,1./3.) + 1)) << std::endl;    
        if(VerboseFlag)std::cout << "width = " << width << std::endl;
        
        result = pars[0] * width / (pow(x[0] - pars[2],2.) + 0.25 * pow(width,2));//all terms in here are in MeV!
        
        if(VerboseFlag)std::cout << "result = " << result << std::endl;
    }
    else
        TF1::RejectPoint();//can't have decays below the threshold
    
    return result;
}

double WidthFunction(double *x, double *pars)
{
    double width = 0;
    if(x[0] - Sp_9B > 0)
    {
        width = pars[0] * 2 * penetrability(4,1,8,1,1,x[0] - Sp_9B,1.25*(pow(8.,1./3.) + 1));
    }
    else
        TF1::RejectPoint();
        
        
    return width;
}
    
    

int main()
{
    FitSpectra();
    
    return 0;
}

void FitSpectra()
{
    TFile *fin = TFile::Open("b_spect_siliconcut.root");
    
    TH1F *hSingles = (TH1F*)fin->Get("Ex_t");
    
    TCanvas *c1 = new TCanvas();
    
    hSingles->Draw();
    
    c1->SetLogy();
    
    TF1 *fitty = new TF1("fitty",FitFunction,Sp_9B+0.01,1,3);
    fitty->SetParameters(1.e5,reduced_width_9B_resonance,0.);
    fitty->Draw("same");
    
//     hSingles->Fit(fitty,"BRMLE");
    
    TCanvas *c2 = new TCanvas();
    
    fitty->Draw("");
    c2->SetLogy();
    
    TCanvas *c3 = new TCanvas();
    c3->SetLogy();
    TF1 *widthF = new TF1("widthF",WidthFunction,Sp_9B+0.01,1,1);
    widthF->SetParameter(0,reduced_width_9B_resonance);
    widthF->Draw();
}
