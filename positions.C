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
#define a2 2.28e-2
#define da2 3.00e-4

TString ToString(int num){
	ostringstream start;
	start<<num;
	TString start1=start.str();
	return start1;	
}

void time_plot(TString filename, TH1D * hist, double * mean, double * stddev){
	//Plots to an histogram the time difference between the 2 signals of the upper photomultipliers. Returns the mean and RMS of the histogram. 
	double trash;
	double temp1;
	double temp2;
	
	std::ifstream myfile;
	myfile.open(filename);

   while ( true ) {
     myfile >> trash;
     myfile >> temp1;
     myfile >> temp2;
     hist -> Fill(temp2/a2);//Convert from ADC to time.
     if( myfile.eof() ) break;
   };
      
   myfile.close();	
   
   *mean = hist -> GetMean();
   *stddev = hist -> GetRMS()/TMath::Sqrt(hist -> GetEntries()-1);
}

int positions(){
	//Plots the time differences of the 2 signals from the upper photomultipliers versus the position of the small scintillator.
	double pos[8];
	double dpos[8];
	double mean[8];
	double stddev[8];
	double dmean[8];
	double dstddev[8];
	TCanvas * canv[8];
	TH1D * histos[8];
	TString filename;
	
	//Cycle on the positions of the small scintillator.
	for(int i = 0; i < 8; i++){
		filename = "data/positions/1G" + ToString(i+16) + "_XCAL_" + ToString(50+i*25) + ".dat";
		pos[i] = 50+i*25;//Positions
		dpos[i] = 0.5;//Position uncertainty
		histos[i] = new TH1D("hist_" + ToString(i), "i", 100, 0, 0.75/a2);//Initialize the histogram
		time_plot(filename, histos[i], &mean[i], &stddev[i]);
		cout << "mean " << mean[i] <<" stddev " << stddev[i] << endl;
		canv[i] = new TCanvas ("c_" + ToString(i), "c_" + ToString(i), 1);//Draw the histogram
		canv[i] -> cd();
		histos[i] -> SetTitle(" ");
		histos[i] -> GetXaxis() -> SetTitle("Time difference [ns]");
		histos[i] -> GetYaxis() -> SetTitle("Number of events");
		histos[i] -> Draw();
	}
	
	TCanvas * c_pos = new TCanvas("c_pos", "c_pos", 1);
	TGraphErrors * gr_pos = new TGraphErrors(8, pos, mean, dpos, stddev);//Fit linearly times differences versus positions
	gr_pos -> SetTitle(" ");
	gr_pos -> GetXaxis() -> SetTitle("Position [cm]");
	gr_pos -> GetYaxis() -> SetTitle("Time difference [ns]");
	TF1 * fit_lin = new TF1("fit_lin", "[0]+[1]*x", 0, 300);
	fit_lin -> SetParameter(0, 25);
	fit_lin -> SetParameter(1, -0.07);
	gr_pos -> Fit("fit_lin", "", "R", 45, 230);
	cout << "Speed of light in the scintillator " << -1/(gr_pos -> GetFunction("fit_lin") -> GetParameter(1)) * 1e7 << 
	" +- " << -(gr_pos -> GetFunction("fit_lin") -> GetParError(1))/(gr_pos -> GetFunction("fit_lin") -> GetParameter(1))/(gr_pos -> GetFunction("fit_lin") -> GetParameter(1)) * 1e7 << endl;
	cout << gr_pos -> GetFunction("fit_lin") -> GetChisquare() << endl;
	gr_pos -> Draw();
	
	
	
	
	
	return 0;
	}
