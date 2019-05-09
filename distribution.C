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
#define c_b 14

TString ToString(int num){
	ostringstream start;
	start<<num;
	TString start1=start.str();
	return start1;	
}

int distribution(){
	double trash;
	double temp1;
	double temp2;
	//double diff;
	//double bin_size=2;
	int count = 0;

	TH1D * hist_tot = new TH1D("hist_tot", "hist_tot", 100, 0, 2000*0.040283*c_b/2);
	TH1D * hist_coin = new TH1D("hist_coin", "hist_coin", 100, 0, 2000*0.040283*c_b/2);
	
	TCanvas * c_1 = new TCanvas("c_1", "c_1", 1);
	TCanvas * c_2 = new TCanvas("c_2", "c_2", 1);
	
	TString filename = "data/distribution/1V6_barcal3.dat";
	//TString filename = "data/distribution/2L9_barcal4.dat";
	
	std::ifstream myfile;
	myfile.open(filename);

   while ( true ) {
     myfile >> trash;
     myfile >> temp1;
     myfile >> temp2;
     //diff = (temp1-0.62197265625)/a1;
     //cout << diff << endl;
     //if (diff > 0. && diff < bin_size) bin_size = diff;
     hist_tot -> Fill(temp1/a1*c_b/2);
     if(temp2/a2!=0) hist_coin -> Fill(temp1/a1*c_b/2);//Convert from ADC to time.
     //count++;
     //if(count == 100000) break;
     if( myfile.eof() ) break;
   };
      
   myfile.close();	
	
	//cout << bin_size;
	
	c_1 -> cd();
	hist_tot -> Draw();
	c_2 -> cd();
	hist_coin -> Draw();
	
	return 0;
}
