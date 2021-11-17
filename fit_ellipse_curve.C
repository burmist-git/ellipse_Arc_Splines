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

////////////////// Circular Arc Splines //////////////////////////
//Approximation of Smooth Planar Curves by Circular Arc Splines //
//http://kaj.uniwersytetradom.pl/docs/Biarcs.pdf                //
////////////////// Circular Arc Splines //////////////////////////

Double_t xMin =-350.0;
Double_t xMax = 350.0;
Double_t yMin =   0.0;
Double_t yMax =  25.0;

/*
Double_t xMin =-500.0;
Double_t xMax = 500.0;
Double_t yMin =-500.0;
Double_t yMax = 500.0;
*/

Bool_t halfEllips(Int_t sign, Double_t x, Double_t a, Double_t b, Double_t xC, Double_t &y);
Bool_t halfEllipsDer(Int_t sign, Double_t x, Double_t a, Double_t b, Double_t xC, Double_t &yDer);

void generate_ellipse(TGraph *gr, Double_t a, Double_t b, Double_t xC, Double_t yC);
void generate_parabola(TGraph *gr, Double_t k, Double_t C, Double_t xMin, Double_t xMax);
void loadData_in_out_all_fit_window(TString dataFileIn, TGraph *gr_window, TGraph *gr_out_window, TGraph *gr_all, Double_t xAbsMin, Double_t xAbsMax);
Bool_t fit(Double_t x0, Double_t y0, Double_t y0Der, Double_t x2, Double_t y2, Double_t y2Der, Double_t &a, Double_t &b, Double_t &xC, Double_t &yC);
Bool_t fit_parabola(Double_t x0, Double_t y0, Double_t y0Der, Double_t x2, Double_t y2, Double_t y2Der, Double_t &k, Double_t &C);

Int_t fit_ellipse_curve( TString dataFileIn = "ellipse.dat"){
  ///////////////////////// LOAD DATA ////////////////////////////
  TGraph *gr_window = new TGraph();
  TGraph *gr_out_window = new TGraph();
  TGraph *gr_all = new TGraph();
  TF1 *f_fit = new TF1();
  //Double_t xAbsMin = 270;
  //Double_t xAbsMax = 298;
  //Double_t xAbsMin = 240;
  //Double_t xAbsMax = 270;
  //Double_t xAbsMin = 210;
  //Double_t xAbsMax = 240;
  //Double_t xAbsMin = 180;
  //Double_t xAbsMax = 210;
  //Double_t xAbsMin = 150;
  //Double_t xAbsMax = 180;
  //Double_t xAbsMin = 150;
  //Double_t xAbsMax = 160;
  //Double_t xAbsMin = 140;
  //Double_t xAbsMax = 150;
  //Double_t xAbsMin = 130;
  //Double_t xAbsMax = 140;
  Double_t xAbsMin = 50.0;
  Double_t xAbsMax = 100.0;
  //Double_t xAbsMin = 50.0;
  //Double_t xAbsMax = 70.0;
  //Double_t xAbsMin = 10.0;
  //Double_t xAbsMax = 50;
  TString outGifFile = "tmp.gif";
  loadData_in_out_all_fit_window("./lens_profile/lens_profile_n1.4_long.dat", gr_window, gr_out_window, gr_all, xAbsMin, xAbsMax);
  //
  Double_t xL;
  Double_t yL;
  Double_t xR;
  Double_t yR;
  gr_all->GetPoint(0,xL,yL);
  gr_all->GetPoint((gr_all->GetN()-1),xR,yR);
  gr_all->Fit("pol8","Q","",xL,xR);
  f_fit = (TF1*)gr_all->GetListOfFunctions()->FindObject("pol8");
  ///////////////////////// LOAD DATA ////////////////////////////
  //
  ///////////////////////// Define fit regions ///////////////////
  Double_t x0 = xAbsMin;
  Double_t x2 = xAbsMax;
  Double_t y0 = f_fit->Eval(x0);
  Double_t y2 = f_fit->Eval(x2);
  Double_t y0Der = f_fit->Derivative(x0);
  Double_t y2Der = f_fit->Derivative(x2);
  ///////////////////////// Define fit regions ///////////////////
  //
  ///////////////////////// fit  /////////////////////////////////
  Double_t a;
  Double_t b;
  Double_t xC;
  Double_t yC;
  Double_t k;
  Double_t C;
  Bool_t fit_result = fit(x0, y0, y0Der, x2, y2, y2Der, a, b, xC, yC);
  Bool_t fit_result_parabola = fit_parabola(x0, y0, y0Der, x2, y2, y2Der, k, C);
  ///////////////////////// fit  /////////////////////////////////
  //
  //////////////////// support graphs ////////////////////////////
  //
  TGraph *gr_support = new TGraph();
  gr_support->SetMarkerStyle(20);
  gr_support->SetMarkerColor(kMagenta);
  gr_support->SetPoint(0,x0,y0);
  gr_support->SetPoint(1,x2,y2);
  //
  TGraph *gr_elips = new TGraph();
  if(fit_result){
    generate_ellipse(gr_elips, a, b, xC, yC);
  }
  TGraph *gr_parabola = new TGraph();
  if(fit_result_parabola){
    generate_parabola(gr_parabola, k, C, -300, 300);
  }
  //////////////////// support graphs ////////////////////////////
  //
  //////////////////// plotting //////////////////////////////////
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,1000);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  gPad->SetGridx();
  gPad->SetGridy();
  //
  c1->SetRightMargin(0.01);
  c1->SetLeftMargin(0.05);
  c1->SetTopMargin(0.01);
  c1->SetBottomMargin(0.05);
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr_window);
  mg->Add(gr_out_window);
  mg->Add(gr_support);
  //mg->Add(gr_all);
  if(fit_result)
    mg->Add(gr_elips);
  if(fit_result_parabola){
    gr_parabola->SetMarkerColor(kBlue+2);
    mg->Add(gr_parabola);
  }
  mg->GetXaxis()->SetLimits(xMin,xMax);
  mg->Draw("AP");
  mg->SetMinimum(yMin);
  mg->SetMaximum(yMax);
  c1->SaveAs(outGifFile.Data());
  //////////////////// plotting //////////////////////////////////
  return 0;
}

Bool_t halfEllips(Int_t sign, Double_t x, Double_t a, Double_t b, Double_t xC, Double_t yC, Double_t &y){
  Double_t det = a*a - (x - xC)*(x - xC);
  if(det>=0.0){
    y = sign*b/a*TMath::Sqrt(det) + yC;
    return true;
  }
  return false;
}

Bool_t halfEllipsDer(Int_t sign, Double_t x, Double_t a, Double_t b, Double_t xC, Double_t &yDer){
  Double_t det = a*a - (x - xC)*(x - xC);
  if(det>0.0){
    yDer = -sign*b/a/TMath::Sqrt(det)*(x-xC);
    return true;
  }
  return false;
}

void loadData_in_out_all_fit_window(TString dataFileIn, TGraph *gr_window, TGraph *gr_out_window, TGraph *gr_all, Double_t xAbsMin, Double_t xAbsMax){
  gr_window->SetTitle("");
  gr_window->SetMarkerColor(kGreen+2);
  gr_out_window->SetTitle("");
  gr_out_window->SetMarkerColor(kRed);
  gr_all->SetTitle("");
  gr_all->SetMarkerColor(kBlack);
  ifstream fileIn (dataFileIn.Data());
  Double_t x;
  Double_t y;  
  if (fileIn.is_open()) {
    while(fileIn>>x>>y){
      if( TMath::Abs(x)>xAbsMin && TMath::Abs(x)<xAbsMax){
	gr_window->SetPoint(gr_window->GetN(),x,y);
      }
      else{
	gr_out_window->SetPoint(gr_out_window->GetN(),x,y);
      }
      gr_all->SetPoint(gr_all->GetN(),x,y);
    }
    fileIn.close();
  }
  else cout<<"Unable to open file"<<endl; 
}

Bool_t fit(Double_t x0, Double_t y0, Double_t y0Der, Double_t x2, Double_t y2, Double_t y2Der, Double_t &a, Double_t &b, Double_t &xC, Double_t &yC){
  Double_t denominator = y0Der*x2 - y2Der*x0;
  if(denominator == 0.0){
    cout<<"denominator = "<<denominator<<endl;
    return false;
  }
  yC = (y0*y0Der*x2 - y2*y2Der*x0)/(denominator);
  if(x2 == 0.0){
    cout<<"x2 = "<<x2<<endl;
    return false;
  }
  Double_t b_div_a_2 = (yC - y2)*y2Der/x2;
  if(b_div_a_2 <= 0.0){
    return false;
  }
  Double_t a2 = (yC-y0)*(yC-y0)/b_div_a_2 + x0*x0;
  if(a2 <= 0.0){
    cout<<"yC          = "<<yC<<endl
	<<"x0          = "<<x0<<endl
	<<"y0          = "<<y0<<endl
	<<"x2          = "<<x2<<endl
	<<"y2          = "<<y2<<endl
	<<"yDer0       = "<<y0Der<<endl
	<<"yDer2       = "<<y2Der<<endl;
    cout<<"denominator = "<<denominator<<endl
	<<"b_div_a_2   = "<<b_div_a_2<<endl
	<<"a2          = "<<a2<<endl;
    return false;
  }
  a = TMath::Sqrt(a2);
  b = a*TMath::Sqrt(b_div_a_2);
  xC = 0;
  cout<<"a  = "<<a<<endl
      <<"b  = "<<b<<endl
      <<"xC = "<<xC<<endl
      <<"yC = "<<yC<<endl;
  return true;
}

void generate_parabola(TGraph *gr, Double_t k, Double_t C, Double_t xMin, Double_t xMax){
  Double_t x;
  Double_t y;
  Int_t nn = 10000;
  for(Int_t i = 0; i<nn;i++){
    x = xMin + (xMax - xMin)/(nn-1)*i;
    y = k*x*x + C;
    gr->SetPoint(gr->GetN(),x,y);
  }
}

void generate_ellipse(TGraph *gr, Double_t a, Double_t b, Double_t xC, Double_t yC){
  Double_t phi_min = 0.0;
  Double_t phi_max = 2.0*TMath::Pi();
  Double_t phi;
  Double_t r;
  Double_t x;
  Double_t y;
  Int_t nn = 1000000;
  for(Int_t i = 0; i<nn;i++){
    phi = phi_min + (phi_max - phi_min)/(nn-1)*i;
    r = a*b/TMath::Sqrt((b*TMath::Cos(phi))*(b*TMath::Cos(phi)) + (a*TMath::Sin(phi))*(a*TMath::Sin(phi)));
    TVector2 v_tmp;
    v_tmp.SetMagPhi(r,phi);
    x = v_tmp.X() + xC;
    y = v_tmp.Y() + yC;
    gr->SetPoint(gr->GetN(),x,y);
  }
}

Bool_t fit_parabola(Double_t x0, Double_t y0, Double_t y0Der, Double_t x2, Double_t y2, Double_t y2Der, Double_t &k, Double_t &C){
  if(x0 == 0.0)
    return false;
  k = y0Der/2/x0;
  C = y2 - k*x2*x2;
  return true;
}
