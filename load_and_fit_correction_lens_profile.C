//root
#include "TVector2.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TString.h"

//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <iomanip>
#include <vector>

#include <time.h>

using namespace std;


string dataFileIn = "./lens_profile/lens_profile_n1.4_long.dat";
string outDataFitParFileIn = "./lens_profile/lens_profile_n1.4_long_pol8_sym_par.dat";

Double_t xGlobalMin =-350.0;
Double_t xGlobalMax = 350.0;
Double_t yGlobalMin =   0.0;
Double_t yGlobalMax =  25.0;

Double_t pol8_symmetric(Double_t *x, Double_t *par);
void load_and_fit_correction_lens_profile_function( string dataFileIn, TGraph *gr);

Int_t load_and_fit_correction_lens_profile(){
  TGraph *gr = new TGraph();
  load_and_fit_correction_lens_profile_function( dataFileIn, gr);
  TF1 *f_fit = (TF1*)gr->GetListOfFunctions()->FindObject("f_pol8_symmetric");
  cout<<f_fit->GetParameter(0)<<endl
      <<f_fit->GetParameter(1)<<endl
      <<f_fit->GetParameter(2)<<endl
      <<f_fit->GetParameter(3)<<endl
      <<f_fit->GetParameter(4)<<endl;
  gr->Draw("APL");
  return 0;
}

Double_t pol8_symmetric(Double_t *x, Double_t *par){
  Double_t p0  = par[0];
  Double_t p2  = par[1];
  Double_t p4  = par[2];
  Double_t p6  = par[3];
  Double_t p8  = par[4];
  return p0 + p2*TMath::Power(x[0],2) + p4*TMath::Power(x[0],4) + p6*TMath::Power(x[0],6) + p8*TMath::Power(x[0],8);
}

void load_and_fit_correction_lens_profile_function( string dataFileIn, TGraph *gr){
  TGraph *grfinal = new TGraph();
  ifstream fileIn (dataFileIn.c_str());
  Double_t x;
  Double_t y;  
  if (fileIn.is_open()){
    while(fileIn>>x>>y){
      gr->SetPoint(gr->GetN(),x,y);
      grfinal->SetPoint(grfinal->GetN(),x,y);
    }
    fileIn.close();
  }
  Double_t xL;
  Double_t yL;
  Double_t xR;
  Double_t yR;
  gr->GetPoint(0,xL,yL);
  gr->GetPoint((gr->GetN()-1),xR,yR);
  gr->Fit("pol8","Q","",xL,xR);
  TF1 *f_tmp_fit = (TF1*)gr->GetListOfFunctions()->FindObject("pol8");
  //fit
  Int_t npar_pol8_symmetric = 5;
  TF1 *f_pol8_symmetric = new TF1("f_pol8_symmetric", pol8_symmetric, xL, xR, npar_pol8_symmetric);
  f_pol8_symmetric->SetParameter(0,f_tmp_fit->GetParameter(0));
  f_pol8_symmetric->SetParameter(1,f_tmp_fit->GetParameter(2));
  f_pol8_symmetric->SetParameter(2,f_tmp_fit->GetParameter(4));
  f_pol8_symmetric->SetParameter(3,f_tmp_fit->GetParameter(6));
  f_pol8_symmetric->SetParameter(4,f_tmp_fit->GetParameter(8));
  //
  gr->Fit("f_pol8_symmetric","Q","",xL,xR);
  TF1 *f_fit = (TF1*)gr->GetListOfFunctions()->FindObject("f_pol8_symmetric");
  FILE *fp01;
  fp01 = fopen(outDataFitParFileIn.c_str(), "w+");
  fprintf(fp01, "%10.5f \n",xL);
  fprintf(fp01, "%10.5f \n",xR);
  for(Int_t ii = 0;ii<npar_pol8_symmetric;ii++)
    fprintf(fp01, "%70.60f \n",f_fit->GetParameter(ii));
  fclose(fp01);
  
  //

}
