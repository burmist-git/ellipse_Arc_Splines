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

#include <time.h>

using namespace std;

Bool_t halfEllips(Int_t sign, Double_t x, Double_t a, Double_t b, Double_t xC, Double_t &y);
Bool_t halfEllipsDer(Int_t sign, Double_t x, Double_t a, Double_t b, Double_t xC, Double_t &yDer);

Int_t generate_ellipse_curve(){
    
  Double_t xMin = -7.0;
  Double_t xMax =  7.0;
  Double_t yMin = -7.0;
  Double_t yMax =  7.0;

  Double_t a = 5.0;
  Double_t b = 2.0;
  Double_t yC = -2.0;

  Double_t x0 = -2.8;
  Double_t y0;
  Double_t y0Der;
  Double_t x2 = -2.6;
  Double_t y2;
  Double_t y2Der;
  halfEllips(1, x0, a, b, xC, y0);
  halfEllipsDer(1, x0, a, b, xC, y0Der);
  halfEllips(1, x2, a, b, xC, y2);
  halfEllipsDer(1, x2, a, b, xC, y2Der);
  
  const Int_t nn = 10000;

  TString outFileName = "ellipse.dat";
  
  ////////////////// INITIAL DATA /////////////////////////
  Double_t xarr[nn];
  Double_t yarr[nn];

  Double_t x = 0.0;
  Double_t y = -999.0;
  Int_t sign = 1;
  Int_t ncount = 0;
  
  for(Int_t i = 0;i<nn;i++){
    x = xMin + (xMax - xMin)/(nn-1)*i;
    if(halfEllips( sign, x, a, b, xC, y)){
      xarr[ncount] = x;
      yarr[ncount] = y;
      ncount++;
    }
  }
  TGraph *gr = new TGraph(ncount,xarr,yarr);
  gr->SetTitle("");
  gr->SetLineWidth(3.0);  

  ///////////////// plotting //////////////////////////////
  TGraph *gr_support = new TGraph();
  gr_support->SetPoint(0,0,0);
  gr_support->SetPoint(1,x0,y0);
  gr_support->SetPoint(2,x2,y2);
  //
  gr_support->SetMarkerStyle(20);  
  gr_support->SetMarkerColor(kMagenta);  
  //
  const Int_t nn_tl_P0 = 10000;
  const Int_t nn_tl_P2 = 10000;
  Double_t x_tl_P0[nn_tl_P0];
  Double_t y_tl_P0[nn_tl_P0]; 
  Double_t x_tl_P2[nn_tl_P2];
  Double_t y_tl_P2[nn_tl_P2]; 
  //
  for(Int_t j = 0;j<nn_tl_P0;j++){
    x_tl_P0[j] = xMin + (xMax - xMin)/(nn_tl_P0-1)*j;
    y_tl_P0[j] = y0Der*(x_tl_P0[j] - x0) + y0; 
  }
  for(Int_t j = 0;j<nn_tl_P2;j++){
    x_tl_P2[j] = xMin + (xMax - xMin)/(nn_tl_P2-1)*j;
    y_tl_P2[j] = y2Der*(x_tl_P2[j] - x2) + y2; 
  }
  TGraph *gr_tl_P0 = new TGraph(nn_tl_P0,x_tl_P0,y_tl_P0);
  gr_tl_P0->SetMarkerColor(kBlue);
  TGraph *gr_tl_P2 = new TGraph(nn_tl_P2,x_tl_P2,y_tl_P2);
  gr_tl_P2->SetMarkerColor(kBlue);

  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,1000);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  gPad->SetGridx();
  gPad->SetGridy();

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr_support);
  mg->Add(gr_tl_P0);
  mg->Add(gr_tl_P2);
  mg->GetXaxis()->SetLimits(xMin,xMax);
  mg->Draw("AP");
  mg->SetMinimum(yMin);
  mg->SetMaximum(yMax);

  ///////////////// plotting //////////////////////////////


  /////////////////////////////////////////////////////////
  FILE * fp;
  fp = fopen (outFileName.Data(), "w+");
  fprintf(fp, "%20.15f %20.15f %20.15f \n",x0,y0,y0Der);
  fprintf(fp, "%20.15f %20.15f %20.15f \n",x2,y2,y2Der);
  fclose(fp);
  /////////////////////////////////////////////////////////

  
  return 0;
}

Bool_t halfEllips(Int_t sign, Double_t x, Double_t a, Double_t b, Double_t xC, Double_t &y){
  Double_t det = a*a - (x - xC)*(x - xC);
  if(det>=0.0){
    y = sign*b/a*TMath::Sqrt(det);
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
