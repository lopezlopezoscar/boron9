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
double Sa_9B = 1.689; //MeV

double reduced_width_9B_resonance_0 = 0.54/penetrability(4,1,8,1,1,-Sp_9B,1.25*(pow(8.,1./3.) + 1)) / 1e3; ///1e3 to convert to MeV - 0.54 keV ground-state width
double reduced_width_9B_resonance_1 = 650/penetrability(4,1,8,1,0,1.86-Sp_9B,1.25*(pow(8.,1./3.) + 1)) / 1e3; //1st excited state 1.85 MeV
double reduced_width_9B_resonance_2 = 81/penetrability(3,2,5,4,2,2.345-Sp_9B,1.25*(pow(5.,1./3.) + pow(4.,1./3.))) / 1e3; //2nd excited state//Using Sp
double reduced_width_9B_resonance_3 = 614/penetrability(4,1,8,1,2,2.75-Sp_9B,1.25*(pow(8.,1./3.) + 1)) / 1e3; //3rd excited state 2.75 MeV
double reduced_width_9B_resonance_4 = 3130/penetrability(4,1,8,1,1,2.79-Sp_9B,1.25*(pow(8.,1./3.) + 1)) / 1e3; //4th excited state 2.79 MeV

//this is a kludge but I need access to the functions outside the main loop to do calculations soooooo
// TF1 *fitty;
// TF1 *fGaus;

double FitFunction(double *x, double *pars)
{
    double result = 0;
    
    //pars[0] = height of the 9B ground state
    //pars[1] = reduced width of the 9B ground state
    //pars[2] = energy of the 9B ground-state resonance
    //pars[3] = height of the excited state
    //pars[4] = reduced width of the excited state
    //pars[5] = energy of the resonance
    
    if(x[0] - Sp_9B > 0)
    {
        // Single states
        
        // Ground State 3/2- l=1
        //double width = pars[1] * 2 * penetrability(4,1,8,1,1,x[0] - Sp_9B,1.25*(pow(8.,1./3.) + 1));
        //result = pars[0] * width / (pow(x[0] - pars[2],2.) + 0.25 * pow(width,2));//all terms in here are in MeV!
        
        // 1st excited state 1/2+ l=0
        //double width1 = pars[1] * 2 * penetrability(4,1,8,1,0,x[0] - Sp_9B,1.25*(pow(8.,1./3.) + 1)); //8Be+p 
        //result += pars[0] * width1 / (pow(x[0] - pars[2],2.) + 0.25 * pow(width1,2));
        
        // 2nd excited state 5/2 l=2
        //double width2 = pars[1] * 2 * penetrability(3,2,5,4,2,x[0] - Sp_9B,1.25*(pow(5.,1./3.) + pow(4.,1./3.))); //5Li + alpha 
        //result += pars[0] * width2 / (pow(x[0] - pars[2],2.) + 0.25 * pow(width2,2));
        
        // 3rd excited state 5/2 l=2
        //double width3 = pars[1] * 2 * penetrability(4,1,8,1,2,x[0] - Sp_9B,1.25*(pow(8.,1./3.) + 1)); //8Be+p 
        //result += pars[0] * width3 / (pow(x[0] - pars[2],2.) + 0.25 * pow(width3,2));
        
        // 4th excited state 1/2 l=1
        double width4 = pars[1] * 2 * penetrability(4,1,8,1,1,x[0] - Sp_9B,1.25*(pow(8.,1./3.) + 1)); //8Be+p 
        result += pars[0] * width4 / (pow(x[0] - pars[2],2.) + 0.25 * pow(width4,2));
    }
    else
        TF1::RejectPoint();//can't have decays below the threshold
    
    return result;
}

int main()
{
    FitSpectra();
    
    return 0;
}

void FitSpectra()
{
    TFile *fin = TFile::Open("b_spect_siliconcut.root");
    
    TH1F *hSingles = (TH1F*)fin->Get("Ex_singles");
    
    TCanvas *c1 = new TCanvas();
    
    hSingles->Draw();
    
    c1->SetLogy();
    
    TF1 *fitty = new TF1("fitty",FitFunction,Sp_9B+0.01,4,3);//1 state 3 parameters -> 2 states 6 parameters
    //fitty->SetParameters(50,reduced_width_9B_resonance_0,0);//GS
    //fitty->SetParameters(300,reduced_width_9B_resonance_1,1.86);//1st
    //fitty->SetParameters(250,reduced_width_9B_resonance_2,2.345);//2nd
    //fitty->SetParameters(650,reduced_width_9B_resonance_3,2.75);//3rd
    fitty->SetParameters(1000,reduced_width_9B_resonance_4,2.79);//4th
    fitty->Draw("same");
    fitty->SetNpx(1e5);

    TFile *foutput = new TFile("FitResults.root","RECREATE");
    
    hSingles->Write();
    fitty->Write();
    foutput->Close();
}
