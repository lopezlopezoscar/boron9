#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TROOT.h>
#include "function_library.h"
#include <TCanvas.h>
#include <TF1Convolution.h>

bool VerboseFlag = false; //print out problems

void FitSpectra();

double Sp_9B = -0.1859; //MeV

double reduced_width_9B_resonance = 0.5*0.54/penetrability(4,1,8,1,1,-Sp_9B,1.25*(pow(8.,1./3.) + 1)) / 1e3; ///1e3 to convert to MeV - 0.54 keV ground-state width

//this is a kludge but I need access to the functions outside the main loop to do calculations soooooo
// TF1 *fitty;
// TF1 *fGaus;

double FitFunction(double *x, double *pars)
{
    double result = 0;
    
    //pars[0] = height of the 9B ground state
    //pars[1] = reduced width of the 9B ground state
    //pars[2] = energy of the 9B ground-state resonance
    //pars[3] =?
    
    if(x[0] - Sp_9B > 0)
    {
        
        double width = pars[1] * 2. * penetrability(4,1,8,1,1,x[0] - Sp_9B,1.25*(pow(8.,1./3.) + 1));
        
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

double GaussianPeak(double *x, double *pars)
{
    //pars[0] = height
    //pars[1] = width
    //NO POSITON - we smear at dawn/around the actual point for each thingy
    
    double result = pars[0] * TMath::Gaus(x[0],0,pars[1]);
    
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
    
// double ConvolvedFunction(double *x, double *pars)
// {
//     TF1Convolution *fConv = new TF1Convolution(fitty,fGaus,false);
//         
//     double result = fConv->EvalNumConv(x[0]);
//     
//     delete fConv;//only you can prevent memory leaks
//     return result;
// }

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
    fitty->SetParameters(1.e3,reduced_width_9B_resonance,0.);
    fitty->Draw("same");
    fitty->SetNpx(1e5);
    
    TF1 *fGaus = new TF1("fGaus",GaussianPeak,-5,5,2);
    fGaus->SetParameters(1.e5,0.02);//guessed these
    fGaus->SetLineColor(3);
    fGaus->Draw("same");
    fGaus->SetNpx(1e5);
    
    //assuming that the spectrum is described by the lineshape combined with a fat old Gaussian for resolution
    TF1Convolution *fConv = new TF1Convolution(fitty,fGaus,false);
    
//     TF1 *fConvolved = new TF1("fConvolved",ConvolvedFunction,Sp_9B+0.01,1,3);
    TF1 *fConvolved = new TF1("fConvolved",*fConv,-5,5,fConv->GetNpar());
    if(VerboseFlag)std::cout << "fConvolved->GetNpar(): " << fConvolved->GetNpar() << std::endl;
    fConvolved->SetParameters(
//                               fitty->GetParameter(0),
                              0.3,//this is an empirical guesstimate :)
                              fitty->GetParameter(1),
                              fitty->GetParameter(2),
                              fGaus->GetParameter(0),
                              fGaus->GetParameter(1)
                );
    
    fConvolved->SetLineColor(4);
    fConvolved->SetNpx(1e5);
    fConvolved->Draw("same");
    
    fConvolved->SetRange(-0.2,1);
    
    hSingles->Fit(fConvolved,"BRLME");//this takes a long time, comment it out if you don't want to fit and are just guessing parameters :)
    
    if(VerboseFlag)std::cout << "fConvolved->Eval(0): " << fConvolved->Eval(0) << std::endl;
    
    TFile *foutput = new TFile("FitResults.root","RECREATE");
    
    TF1 *widthF = new TF1("widthF",WidthFunction,Sp_9B+0.01,1,1);
    widthF->SetParameter(0,reduced_width_9B_resonance);
    
    hSingles->Write();
    fitty->Write();
    fGaus->Write();
    fConv->Write();
    fConvolved->Write();
    widthF->Write();
    
    foutput->Close();
    
//     hSingles->Fit(fitty,"BRMLE");
    
//     TCanvas *c2 = new TCanvas();
    
//     fitty->Draw("");
//     c2->SetLogy();
    
//     TCanvas *c3 = new TCanvas();
//     c3->SetLogy();
    
//     widthF->Draw();
}
