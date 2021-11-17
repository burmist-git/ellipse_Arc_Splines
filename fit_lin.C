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

Double_t xGlobalMin =-350.0;
Double_t xGlobalMax = 350.0;
Double_t yGlobalMin =   0.0;
Double_t yGlobalMax =  25.0;

struct fit_state_str {
  Double_t x0;
  Double_t x2;
  Double_t y0;
  Double_t y2;
  Double_t y0Der;
  Double_t y2Der;
  //
  TString fitCurveTypeName;
  //
  Double_t a;
  Double_t b;
  Double_t xC;
  Double_t yC;
  Double_t k;
  Double_t C;
  Double_t R;
  //
  Double_t fit_x0;
  Double_t fit_x2;
  Double_t fit_y0;
  Double_t fit_y2;
  Double_t fit_y0Der;
  Double_t fit_y2Der;
  //
  Bool_t fit_result;
  //
  Double_t pol0_val;
  //
  vector<Double_t> vec_x_fit;
  vector<Double_t> vec_y_fit;
  //
  fit_state_str(){
    x0    = -999.0;
    x2    = -999.0;
    y0    = -999.0;
    y2    = -999.0;
    y0Der = -999.0;
    y2Der = -999.0;
    fitCurveTypeName = "";
    a  = -999.0;
    b  = -999.0;
    xC = -999.0;
    yC = -999.0;
    k  = -999.0;
    C  = -999.0;
    R  = -999.0;
    //
    fit_x0    = -999.0;
    fit_x2    = -999.0;
    fit_y0    = -999.0;
    fit_y2    = -999.0;
    fit_y0Der = -999.0;
    fit_y2Der = -999.0;
    //
    fit_result = false;
    //
    pol0_val = -999.0;
  }
  void print_state(){
    cout<<"x0               "<<x0<<endl
	<<"x2               "<<x2<<endl
	<<"y0               "<<y0<<endl;
    cout<<"y2               "<<y2<<endl
	<<"y0Der            "<<y0Der<<endl
	<<"y2Der            "<<y2Der<<endl
	<<"fitCurveTypeName "<<fitCurveTypeName<<endl
	<<"a                "<<a<<endl;
    cout<<"b                "<<b<<endl
	<<"xC               "<<xC<<endl
	<<"yC               "<<yC<<endl;
    cout<<"k                "<<k<<endl
	<<"C                "<<C<<endl
	<<"R                "<<R<<endl
	<<"fit_result       "<<fit_result<<endl;
    cout<<"fit_x0           "<<fit_x0<<endl
	<<"fit_x2           "<<fit_x2<<endl
	<<"fit_y0           "<<fit_y0<<endl
	<<"fit_y2           "<<fit_y2<<endl
	<<"fit_y0Der        "<<fit_y0Der<<endl
	<<"fit_y2Der        "<<fit_y2Der<<endl;
    cout<<"pol0_val         "<<pol0_val<<endl;
    cout<<"vec_x_fit.size() "<<vec_x_fit.size()<<endl
	<<"vec_y_fit.size() "<<vec_y_fit.size()<<endl;
  }
};

Bool_t halfEllips(Int_t sign, Double_t x, Double_t a, Double_t b, Double_t xC, Double_t yC, Double_t &y);
Bool_t halfEllipsDer(Int_t sign, Double_t x, Double_t a, Double_t b, Double_t xC, Double_t &yDer);
//
Bool_t halfCircle(Int_t sign, Double_t x, Double_t R, Double_t xC, Double_t yC, Double_t &y);
Bool_t halfCircleDer(Int_t sign, Double_t x, Double_t R, Double_t xC, Double_t &yDer);
//
Bool_t parabola(Double_t x, Double_t k, Double_t C, Double_t &y);
Bool_t parabolaDer(Double_t x, Double_t k, Double_t &yDer);
//
void generate_ellipse(vector<Double_t> &vec_x_fit, vector<Double_t> &vec_y_fit, Double_t a, Double_t b, Double_t xC, Double_t yC, Double_t xMin, Double_t xMax);
void generate_circle(vector<Double_t> &vec_x_fit, vector<Double_t> &vec_y_fit, Double_t R, Double_t xC, Double_t yC, Double_t xMin, Double_t xMax);
void generate_parabola(vector<Double_t> &vec_x_fit, vector<Double_t> &vec_y_fit, Double_t k, Double_t C, Double_t xMin, Double_t xMax);
void generate_pol0(vector<Double_t> &vec_x_fit, vector<Double_t> &vec_y_fit, Double_t pol0_val, Double_t xMin, Double_t xMax);
//
void loadData(TString dataFileIn, TGraph *gr);
//
Bool_t fit_ellipse(Double_t x0, Double_t y0, Double_t y0Der, Double_t x2, Double_t y2, Double_t y2Der, Double_t &a, Double_t &b, Double_t &xC, Double_t &yC);
Bool_t fit_parabola(Double_t x0, Double_t y0, Double_t y0Der, Double_t x2, Double_t y2, Double_t y2Der, Double_t &k, Double_t &C);
Bool_t fit_circle(Double_t x0, Double_t y0, Double_t y0Der, Double_t x2, Double_t y2, Double_t y2Der, Double_t &R, Double_t &xC, Double_t &yC);
Bool_t fit_pol0(Double_t x0, Double_t y0, Double_t x2, Double_t &pol0_val);
Bool_t fit( Double_t x0, Double_t y0, Double_t y0Der,
	    Double_t x2, Double_t y2, Double_t y2Der,
	    TString fitCurveTypeName, Double_t xMin, Double_t xMax, fit_state_str *fit_state);
//
void plot(TGraph *gr, const vector<fit_state_str*> &v_fit_state_str);
//
Int_t fit_ellipse_curve_parabola_circle(){
  TString dataFileIn = "./lens_profile/lens_profile_n1.4_long.dat";
  //
  ///////////////////////// LOAD DATA ////////////////////////////
  TGraph *gr = new TGraph();
  loadData(dataFileIn, gr);
  Double_t xL;
  Double_t yL;
  Double_t xR;
  Double_t yR;
  gr->GetPoint(0,xL,yL);
  gr->GetPoint((gr->GetN()-1),xR,yR);
  gr->Fit("pol8","Q","",xL,xR);
  TF1 *f_fit = (TF1*)gr->GetListOfFunctions()->FindObject("pol8");
  ///////////////////////// FIT //////////////////////////////////
  TString fitCurveTypeName;
  Double_t x0;
  Double_t y0;
  Double_t y0Der;
  Double_t x2;
  Double_t y2;
  Double_t y2Der;
  vector<fit_state_str*> v_fit_state_str;
  //
  Double_t dx = 10;
  //01
  fitCurveTypeName = "ellipse";
  x0 = 290;
  x2 = 299;
  y0 = f_fit->Eval(x0);
  y0Der = f_fit->Derivative(x0);
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state))
    v_fit_state_str.push_back(fit_state);
  //for(Int_t kk = 1;kk<18;kk++){
  for(Int_t kk = 1;kk<18;kk++){
    fitCurveTypeName = "ellipse";
    x2 = v_fit_state_str.at(kk-1)->x0;
    x0 = x2 - dx;
    y0 = f_fit->Eval(x0);
    y0Der = f_fit->Derivative(x0);
    y2 = v_fit_state_str.at(kk-1)->fit_y0;
    y2Der = v_fit_state_str.at(kk-1)->fit_y0Der;
    fit_state_str *fit_state_tmp = new fit_state_str();
    if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state_tmp))
      v_fit_state_str.push_back(fit_state_tmp);
  }
  //
  fitCurveTypeName = "parabola";
  x0 = 90;
  x2 = 100;
  y0 = f_fit->Eval(x0);
  y0Der = f_fit->Derivative(x0);
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state_parabola_first = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state_parabola_first)){
    x0 = 90;
    x2 = TMath::Sqrt((v_fit_state_str.at(v_fit_state_str.size()-1)->fit_y0 - fit_state_parabola_first->C)/fit_state_parabola_first->k);
    y0 = f_fit->Eval(x0);
    y0Der = f_fit->Derivative(x0);
    y2 = f_fit->Eval(x2);
    y2Der = f_fit->Derivative(x2);
    fit_state_str *fit_state_parabola = new fit_state_str();
    if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state_parabola))
      v_fit_state_str.push_back(fit_state_parabola);
  }
  //
  fitCurveTypeName = "pol0";
  x0 = TMath::Sqrt((v_fit_state_str.at(v_fit_state_str.size()-2)->fit_y0 - fit_state_parabola_first->C)/fit_state_parabola_first->k);
  x2 = v_fit_state_str.at(v_fit_state_str.size()-2)->x0;
  y0 = v_fit_state_str.at(v_fit_state_str.size()-2)->fit_y0;
  y0Der = 0.0;
  y2 = v_fit_state_str.at(v_fit_state_str.size()-2)->fit_y0;
  y2Der = 0.0;
  fit_state_str *fit_state_pol0 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state_pol0))
      v_fit_state_str.push_back(fit_state_pol0);

  /*
  //01
  fitCurveTypeName = "ellipse";
  x0 = 280;
  x2 = 299;
  y0 = f_fit->Eval(x0);
  y0Der = f_fit->Derivative(x0);
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state01 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state01))
    v_fit_state_str.push_back(fit_state01);
  //02
  fitCurveTypeName = "ellipse";
  x0 = 260;
  x2 = 280;
  y0 = f_fit->Eval(x0);
  y0Der = f_fit->Derivative(x0);
  y2 = fit_state01->fit_y0;
  y2Der = fit_state01->fit_y0Der;
  fit_state_str *fit_state02 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state02))
    v_fit_state_str.push_back(fit_state02);
  //03
  fitCurveTypeName = "ellipse";
  x0 = 240;
  x2 = 260;
  y0 = f_fit->Eval(x0);
  y0Der = f_fit->Derivative(x0);
  y2 = fit_state02->fit_y0;
  y2Der = fit_state02->fit_y0Der;
  fit_state_str *fit_state03 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state03))
    v_fit_state_str.push_back(fit_state03);
  */
  /*
  //01
  fitCurveTypeName = "ellipse";
  x0 = 160;
  x2 = 180;
  y0 = f_fit->Eval(x0);
  y0Der = f_fit->Derivative(x0);
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state01 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state01))
    v_fit_state_str.push_back(fit_state01);
  //02
  fitCurveTypeName = "ellipse";
  x0 = 180;
  x2 = 200;
  y0 = fit_state01->fit_y2;
  y0Der = fit_state01->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state02 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state02))
    v_fit_state_str.push_back(fit_state02);
  //03
  fitCurveTypeName = "ellipse";
  x0 = 200;
  x2 = 220;
  y0 = fit_state02->fit_y2;
  y0Der = fit_state02->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state03 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state03))
    v_fit_state_str.push_back(fit_state03);
  //04
  fitCurveTypeName = "ellipse";
  x0 = 220;
  x2 = 240;
  y0 = fit_state03->fit_y2;
  y0Der = fit_state03->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state04 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state04))
    v_fit_state_str.push_back(fit_state04);
  //05
  fitCurveTypeName = "ellipse";
  x0 = 240;
  x2 = 260;
  y0 = fit_state04->fit_y2;
  y0Der = fit_state04->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state05 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state05))
    v_fit_state_str.push_back(fit_state05);
  //06
  fitCurveTypeName = "ellipse";
  x0 = 260;
  x2 = 280;
  y0 = fit_state05->fit_y2;
  y0Der = fit_state05->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state06 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state06))
    v_fit_state_str.push_back(fit_state06);
  //07
  fitCurveTypeName = "ellipse";
  x0 = 280;
  x2 = 299;
  y0 = fit_state06->fit_y2;
  y0Der = fit_state06->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state07 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state07))
    v_fit_state_str.push_back(fit_state07);
  */
  
  /*
  //01
  fitCurveTypeName = "parabola";
  x0 = 1;
  x2 = 30;
  y0 = f_fit->Eval(x0);
  y0Der = f_fit->Derivative(x0);
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state01 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state01))
    v_fit_state_str.push_back(fit_state01);
  //02
  fitCurveTypeName = "parabola";
  x0 = 30;
  x2 = 50;
  y0 = fit_state01->fit_y2;
  y0Der = fit_state01->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state02 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state02))
    v_fit_state_str.push_back(fit_state02);
  //03
  fitCurveTypeName = "circle";
  x0 = 50;
  x2 = 130;
  y0 = fit_state02->fit_y2;
  y0Der = fit_state02->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state03 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state03))
    v_fit_state_str.push_back(fit_state03);
  //04
  fitCurveTypeName = "ellipse";
  x0 = 130;
  x2 = 140;
  y0 = fit_state03->fit_y2;
  y0Der = fit_state03->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state04 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x0, fit_state04))
    v_fit_state_str.push_back(fit_state04);
  //05
  fitCurveTypeName = "ellipse";
  x0 = 140;
  x2 = 150;
  y0 = fit_state04->fit_y2;
  y0Der = fit_state04->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state05 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x0, fit_state05))
    v_fit_state_str.push_back(fit_state05);
  //06
  fitCurveTypeName = "ellipse";
  x0 = 150;
  x2 = 160;
  y0 = fit_state05->fit_y2;
  y0Der = fit_state05->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state06 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x0, fit_state06))
    v_fit_state_str.push_back(fit_state06);
  //07
  fitCurveTypeName = "ellipse";
  x0 = 160;
  x2 = 170;
  y0 = fit_state06->fit_y2;
  y0Der = fit_state06->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state07 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x0, fit_state07))
    v_fit_state_str.push_back(fit_state07);
  //08
  fitCurveTypeName = "ellipse";
  x0 = 170;
  x2 = 190;
  y0 = fit_state07->fit_y2;
  y0Der = fit_state07->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state08 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x0, fit_state08))
    v_fit_state_str.push_back(fit_state08);
  */
  
  /*
  //03
  fitCurveTypeName = "ellipse";
  x0 = 270;
  x2 = 299;
  y0 = f_fit->Eval(x0);
  y0Der = f_fit->Derivative(x0);
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Derivative(x2);
  fit_state_str *fit_state03 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x0, fit_state03))
    v_fit_state_str.push_back(fit_state03);
  */
  
  /*
  //05
  fitCurveTypeName = "ellipse";
  xMin = 240;
  xMax = 270;
  fit_state_str *fit_state05 = new fit_state_str();
  if(fit( f_fit, fitCurveTypeName, xMin, xMax, fit_state05))
    v_fit_state_str.push_back(fit_state05);
  //06
  fitCurveTypeName = "ellipse";
  xMin = 210;
  xMax = 240;
  fit_state_str *fit_state06 = new fit_state_str();
  if(fit( f_fit, fitCurveTypeName, xMin, xMax, fit_state06))
    v_fit_state_str.push_back(fit_state06);
  //07
  fitCurveTypeName = "ellipse";
  xMin = 180;
  xMax = 210;
  fit_state_str *fit_state07 = new fit_state_str();
  if(fit( f_fit, fitCurveTypeName, xMin, xMax, fit_state07))
    v_fit_state_str.push_back(fit_state07);
  //08
  fitCurveTypeName = "ellipse";
  xMin = 150;
  xMax = 180;
  fit_state_str *fit_state08 = new fit_state_str();
  if(fit( f_fit, fitCurveTypeName, xMin, xMax, fit_state08))
    v_fit_state_str.push_back(fit_state08);
  */



  /*
  //02
  fitCurveTypeName = "parabola";
  x0 = 20;
  x2 = 80;
  y0 = fit_state01->fit_y2;
  y0Der = fit_state01->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Eval(x2);
  fit_state_str *fit_state02 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state02))
    v_fit_state_str.push_back(fit_state02);
  //03
  fitCurveTypeName = "circle";
  x0 = 80;
  x2 = 200;
  y0 = fit_state02->fit_y2;
  y0Der = fit_state02->fit_y2Der;
  y2 = f_fit->Eval(x2);
  y2Der = f_fit->Eval(x2);
  fit_state_str *fit_state03 = new fit_state_str();
  if(fit( x0, y0, y0Der, x2, y2, y2Der, fitCurveTypeName, x0, x2, fit_state03))
    v_fit_state_str.push_back(fit_state03);
  */
  /*
  //04
  fitCurveTypeName = "parabola";
  xMin = 60;
  xMax = 80;
  fit_state_str *fit_state04 = new fit_state_str();
  if(fit( f_fit, fitCurveTypeName, xMin, xMax, fit_state04))
    v_fit_state_str.push_back(fit_state04);
  //05
  fitCurveTypeName = "parabola";
  xMin = 80;
  xMax = 100;
  fit_state_str *fit_state05 = new fit_state_str();
  if(fit( f_fit, fitCurveTypeName, xMin, xMax, fit_state05))
    v_fit_state_str.push_back(fit_state05);
  //06
  fitCurveTypeName = "parabola";
  xMin = 100;
  xMax = 110;
  fit_state_str *fit_state06 = new fit_state_str();
  if(fit( f_fit, fitCurveTypeName, xMin, xMax, fit_state06))
    v_fit_state_str.push_back(fit_state06);
  //07
  fitCurveTypeName = "circle";
  xMin = 110;
  xMax = 130;
  fit_state_str *fit_state07 = new fit_state_str();
  if(fit( f_fit, fitCurveTypeName, xMin, xMax, fit_state07))
    v_fit_state_str.push_back(fit_state07);
  */
  //
  ///////////////////////// PLOT /////////////////////////////////
  plot(gr,v_fit_state_str);
  return 0;
}

Bool_t fit( Double_t x0, Double_t y0, Double_t y0Der,
	    Double_t x2, Double_t y2, Double_t y2Der,
	    TString fitCurveTypeName, Double_t xMin, Double_t xMax, fit_state_str *fit_state){
  ///////////////////////// fit //////////////////////////////////
  Double_t a = -999.0;
  Double_t b = -999.0;
  Double_t xC = -999.0;
  Double_t yC = -999.0;
  Double_t k = -999.0;
  Double_t C = -999.0;
  Double_t R = -999.0;
  Double_t pol0_val = -999.0;
  Bool_t fit_result = false;
  if(fitCurveTypeName == "ellipse"){
    fit_result = fit_ellipse(x0, y0, y0Der, x2, y2, y2Der, a, b, xC, yC);
    if(fit_result){
      generate_ellipse(fit_state->vec_x_fit, fit_state->vec_y_fit, a, b, xC, yC, x0, x2);
      Double_t ytmp;
      Double_t yDertmp;
      fit_state->fit_x0 = x0;
      fit_state->fit_x2 = x2;
      if(halfEllips(-1, x0, a, b, xC, yC, ytmp))
	fit_state->fit_y0 = ytmp;
      if(halfEllips(-1, x2, a, b, xC, yC, ytmp))
	fit_state->fit_y2 = ytmp;
      if(halfEllipsDer(-1, x0, a, b, xC, yDertmp))
	fit_state->fit_y0Der = yDertmp;
      if(halfEllipsDer(-1, x2, a, b, xC, yDertmp))
	fit_state->fit_y2Der = yDertmp;
    }
    else
      return false;
  }
  else if(fitCurveTypeName == "parabola"){
    fit_result = fit_parabola(x0, y0, y0Der, x2, y2, y2Der, k, C);
    if(fit_result){
      //generate_parabola(fit_state->vec_x_fit, fit_state->vec_y_fit, k, C, x0, x2);
      generate_parabola(fit_state->vec_x_fit, fit_state->vec_y_fit, k, C, 0.0, x2);
      Double_t ytmp;
      Double_t yDertmp;
      fit_state->fit_x0 = x0;
      fit_state->fit_x2 = x2;
      if(parabola(x0, k, C, ytmp))
	fit_state->fit_y0 = ytmp;
      if(parabola(x2, k, C, ytmp))
	fit_state->fit_y2 = ytmp;
      if(parabolaDer(x0, k, yDertmp))
	fit_state->fit_y0Der = yDertmp;
      if(parabolaDer(x2, k, yDertmp))
	fit_state->fit_y2Der = yDertmp;
    }
    else
      return false;
  }
  else if(fitCurveTypeName == "circle"){
    fit_result = fit_circle(x0, y0, y0Der, x2, y2, y2Der, R, xC, yC);
    if(fit_result){
      generate_circle(fit_state->vec_x_fit, fit_state->vec_y_fit, R, xC, yC, x0, x2);
      Double_t ytmp;
      Double_t yDertmp;
      fit_state->fit_x0 = x0;
      fit_state->fit_x2 = x2;
      if(halfCircle(-1, x0, R, xC, yC, ytmp))
	fit_state->fit_y0 = ytmp;
      if(halfCircle(-1, x2, R, xC, yC, ytmp))
	fit_state->fit_y2 = ytmp;
      if(halfCircleDer(-1, x0, R, xC, yDertmp))
	fit_state->fit_y0Der = yDertmp;
      if(halfCircleDer(-1, x2, R, xC, yDertmp))
	fit_state->fit_y2Der = yDertmp;
    }
    else
      return false;
  }
  else if(fitCurveTypeName == "pol0"){
    fit_result = fit_pol0(x0, y0, x2, pol0_val);
    if(fit_result){
      generate_pol0(fit_state->vec_x_fit, fit_state->vec_y_fit, pol0_val, x0, x2);
      Double_t ytmp;
      Double_t yDertmp;
      fit_state->fit_x0 = x0;
      fit_state->fit_x2 = x2;
      fit_state->pol0_val = pol0_val;
    }
    else
      return false;
  }
  else
    return false;
  //
  fit_state->x0 = x0;
  fit_state->x2 = x2;
  fit_state->y0 = y0;
  fit_state->y2 = y2;
  fit_state->y0Der = y0Der;
  fit_state->y2Der = y2Der;
  fit_state->fitCurveTypeName = fitCurveTypeName;
  fit_state->a = a;
  fit_state->b = b;
  fit_state->xC = xC;
  fit_state->yC = yC;
  fit_state->k = k;
  fit_state->C = C;
  fit_state->R = R;
  //
  fit_state->fit_result = fit_result;
  ///////////////////////// fit  /////////////////////////////////
  return true;
}

void plot(TGraph *gr, const vector<fit_state_str*> &v_fit_state_str){
  /*
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
  */
  //////////////////// support graphs ////////////////////////////
  vector<TGraph*> v_gr_support;
  for (unsigned int i=0; i<v_fit_state_str.size(); i++){
    v_fit_state_str.at(i)->print_state();
    TGraph *grtmp = new TGraph();
    for(unsigned int ii = 0; ii< v_fit_state_str.at(i)->vec_x_fit.size();ii++)
      grtmp->SetPoint(grtmp->GetN(),v_fit_state_str.at(i)->vec_x_fit.at(ii),v_fit_state_str.at(i)->vec_y_fit.at(ii));
    v_gr_support.push_back(grtmp);
  }
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
  mg->Add(gr);
  for(unsigned int ii = 0; ii< v_gr_support.size();ii++)
    mg->Add(v_gr_support.at(ii));
  mg->GetXaxis()->SetLimits(xGlobalMin,xGlobalMax);
  mg->Draw("AP");
  mg->SetMinimum(yGlobalMin);
  mg->SetMaximum(yGlobalMax);
  //c1->SaveAs(outGifFile.Data());
  //////////////////// plotting //////////////////////////////////
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

Bool_t halfCircle(Int_t sign, Double_t x, Double_t R, Double_t xC, Double_t yC, Double_t &y){
  return halfEllips(sign, x, R, R, xC, yC, y);
}

Bool_t halfCircleDer(Int_t sign, Double_t x, Double_t R, Double_t xC, Double_t &yDer){
  return halfEllipsDer(sign, x, R, R, xC, yDer);
}

Bool_t parabola(Double_t x, Double_t k, Double_t C, Double_t &y){
  y = k*x*x+C;
  return true;
}

Bool_t parabolaDer(Double_t x, Double_t k, Double_t &yDer){
  yDer = 2*k*x;
  return true;
}

void loadData(TString dataFileIn, TGraph *gr){
  gr->SetTitle("");
  gr->SetMarkerColor(kBlack);
  ifstream fileIn (dataFileIn.Data());
  Double_t x;
  Double_t y;  
  if (fileIn.is_open()){
    while(fileIn>>x>>y)
      gr->SetPoint(gr->GetN(),x,y);
    fileIn.close();
  }
  else cout<<"Unable to open file"<<endl; 
}

Bool_t fit_pol0(Double_t x0, Double_t y0, Double_t x2, Double_t &pol0_val){
  pol0_val = y0;
  return true;
}

Bool_t fit_ellipse(Double_t x0, Double_t y0, Double_t y0Der, Double_t x2, Double_t y2, Double_t y2Der, Double_t &a, Double_t &b, Double_t &xC, Double_t &yC){
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
  //Double_t b_div_a_2 = (yC - y2)*y2Der/x2;
  Double_t b_div_a_2 = (yC - y0)*y0Der/x0;
  if(b_div_a_2 <= 0.0){
    cout<<"yC          = "<<yC<<endl
	<<"x0          = "<<x0<<endl
	<<"y0          = "<<y0<<endl
	<<"x2          = "<<x2<<endl
	<<"y2          = "<<y2<<endl
	<<"yDer0       = "<<y0Der<<endl
	<<"yDer2       = "<<y2Der<<endl;
    cout<<"b_div_a_2 = "<<b_div_a_2<<endl;
    return false;
  }
  //Double_t a2 = (yC-y0)*(yC-y0)/b_div_a_2 + x0*x0;
  Double_t a2 = (yC-y2)*(yC-y2)/b_div_a_2 + x2*x2;
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
  //cout<<"a  = "<<a<<endl
  //    <<"b  = "<<b<<endl
  //    <<"xC = "<<xC<<endl
  //    <<"yC = "<<yC<<endl;
  return true;
}

void generate_parabola(vector<Double_t> &vec_x_fit, vector<Double_t> &vec_y_fit, Double_t k, Double_t C, Double_t xMin, Double_t xMax){
  Double_t x;
  Double_t y;
  Int_t nn = 100000;
  for(Int_t i = 0; i<nn;i++){
    x = xMin + (xMax - xMin)/(nn-1)*i;
    y = k*x*x + C;
    vec_x_fit.push_back(x);
    vec_y_fit.push_back(y);
  }
}

void generate_pol0(vector<Double_t> &vec_x_fit, vector<Double_t> &vec_y_fit, Double_t pol0_val, Double_t xMin, Double_t xMax){
  Double_t x;
  Double_t y;
  Int_t nn = 10000;
  for(Int_t i = 0; i<nn;i++){
    x = xMin + (xMax - xMin)/(nn-1)*i;
    y = pol0_val;
    vec_x_fit.push_back(x);
    vec_y_fit.push_back(y);
  }
}

void generate_ellipse(vector<Double_t> &vec_x_fit, vector<Double_t> &vec_y_fit, Double_t a, Double_t b, Double_t xC, Double_t yC, Double_t xMin, Double_t xMax){
  Double_t phi_min = 0.0;
  Double_t phi_max = 2.0*TMath::Pi();
  Double_t phi;
  Double_t r;
  Double_t x;
  Double_t y;
  Int_t nn = 10000000;
  for(Int_t i = 0; i<nn;i++){
    phi = phi_min + (phi_max - phi_min)/(nn-1)*i;
    r = a*b/TMath::Sqrt((b*TMath::Cos(phi))*(b*TMath::Cos(phi)) + (a*TMath::Sin(phi))*(a*TMath::Sin(phi)));
    TVector2 v_tmp;
    v_tmp.SetMagPhi(r,phi);
    x = v_tmp.X() + xC;
    y = v_tmp.Y() + yC;
    if(phi>=TMath::Pi() && phi<=2*TMath::Pi()){
      //if(phi>=0 && phi<=TMath::Pi()){
      if(x>=xMin && x<=xMax){
	vec_x_fit.push_back(x);
	vec_y_fit.push_back(y);
      }
    }
  }
}

void generate_circle(vector<Double_t> &vec_x_fit, vector<Double_t> &vec_y_fit, Double_t R, Double_t xC, Double_t yC, Double_t xMin, Double_t xMax){
  generate_ellipse(vec_x_fit, vec_y_fit, R, R, xC, yC, xMin, xMax);
}

Bool_t fit_parabola(Double_t x0, Double_t y0, Double_t y0Der, Double_t x2, Double_t y2, Double_t y2Der, Double_t &k, Double_t &C){
  if(x0 == 0.0)
    return false;
  k = y0Der/2/x0;
  C = y0 - k*x0*x0;
  return true;
}

Bool_t fit_circle(Double_t x0, Double_t y0, Double_t y0Der, Double_t x2, Double_t y2, Double_t y2Der, Double_t &R, Double_t &xC, Double_t &yC){
  Double_t denominator = y2Der - y0Der;
  if(denominator == 0.0)
    return false;
  yC = (x2 - x0 - y0*y0Der + y2*y2Der)/(denominator);
  xC = (y0 - yC)*y0Der + x0;
  R = TMath::Sqrt((y0-yC)*(y0-yC) + (x0-xC)*(x0-xC));
  return true;
}
