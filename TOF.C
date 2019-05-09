#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TRandom3.h>
#include <TMinuit.h>
#include <TFile.h>
#include <TLegend.h>
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TCanvas.h"

#define a1 2.00e-2 //Linear conversion factors from time to ADC counts. Produced by timecal.C
#define da1 2.00e-4
#define of1 -3.2e-1
#define dof1 2e-2
#define a2 2.28e-2
#define da2 3.00e-4
#define of2 -3.5e-1
#define dof2 2e-2
#define a3 -1.44e-1 //Linear conversion factors from position to time. Produced by positions.C
#define da3 2e-3
#define of3 4.89e1
#define dof3 2e-1
#define c_b 14
#define sp_res 0.85
#define Zu 175
#define Dcx 140


TString ToString(int num){
	ostringstream start;
	start<<num;
	TString start1=start.str();
	return start1;	
}


int TOF(){

	double trash;
	double temp1;
	double temp2;
	double tof;
	double dist;
	
	TString filename = "data/distribution/1V6_barcal3.dat";

	TH1D * hist_tof = new TH1D("hist_tof", "hist_tof", 1400, -20, 50);
	TH1D * hist_tof_cen = new TH1D("hist_tof_cen", "hist_tof_cen", 1400, -20, 50);
	TH1D * hist_dist = new TH1D("hist_dist", "hist_dis", 100, -100, 400);
	TH1D * hist_beta = new TH1D("hist_beta", "hist_beta", 1000, 0, 10);
	TH2D * hist_dist_tof = new TH2D("hist_dist_tof", "hist_dist_tof", 120, -30, 330, 50, -5, 20);
	
	std::ifstream myfile;
	myfile.open(filename);

   while ( true ) {
     myfile >> trash;
     myfile >> temp1;
     myfile >> temp2;
     if(temp2 != 0){
		//tof = temp2/a2 - temp1/2/a1 - of2/a2 + of1/2/a1 + of3/2 - 15;
		tof = temp2/a2 - temp1/2/a1 - 2;
		//dist = temp1/a1/a3 - of1/a1/a3 -of3/a3 -30/a3;
		dist = -temp1/a1/a3-46;
		hist_dist -> Fill(dist);
		hist_dist_tof -> Fill(dist, tof);
		if ( dist > 250 && dist < 500) hist_tof_cen -> Fill(tof);
		dist = TMath::Sqrt((dist-140)*(dist-140)+Zu*Zu);
		hist_tof -> Fill(tof); 
		hist_beta -> Fill(dist/tof/30);
	 }
     if( myfile.eof() ) break;
   };
   
   TCanvas * c1 = new TCanvas("c1", "c1", 1);
   c1 -> cd();
   hist_tof -> GetXaxis() -> SetTitle("Time of flight [ns]");
   hist_tof -> GetYaxis() -> SetTitle("Number of events");
   hist_tof -> Draw();
   TCanvas * c2 = new TCanvas("c2", "c2", 1);
   c2 -> cd();
   hist_dist -> GetXaxis() -> SetTitle("Position on the scintillator [cm]");
   hist_dist -> GetYaxis() -> SetTitle("Number of events");
   hist_dist -> Draw();
   TCanvas * c3 = new TCanvas("c3", "c3", 1);
   c3 -> cd();
   hist_beta -> GetXaxis() -> SetTitle("#beta_{#mu}");
   hist_beta -> GetYaxis() -> SetTitle("Number of events");
   hist_beta -> Draw();
   TCanvas * c4 = new TCanvas("c4", "c4", 1);
   c4 -> cd();
   hist_tof_cen -> Draw();
   TCanvas * c5 = new TCanvas("c5", "c5", 1);
   c5 -> cd();
   hist_dist_tof -> GetXaxis() -> SetTitle("Position on the scintillator [cm]");
   hist_dist_tof -> GetYaxis() -> SetTitle("Time of flight [ns]");
   hist_dist_tof -> SetStats(0);
   hist_dist_tof -> Draw("COLZ");
   
return 0;
}
