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

TString ToString(int num){
	ostringstream start;
	start<<num;
	TString start1=start.str();
	return start1;	
}

void averager(TString filename ,double * t1, double * t2){
	//Takes a file with three columns and calculates the mean of the last two returning them in variables *t1 and *t2
	double trash;
	double count = 0;
	double temp1;
	double temp2;
	
	std::ifstream myfile;
	myfile.open(filename);

   while ( true ) {
     myfile >> trash;
     myfile >> temp1;
     myfile >> temp2;
     *t1 += temp1;
     *t2 += temp2;
     count ++;
     if( myfile.eof() ) break;
   };
      
   myfile.close();
	
	*t1 = *t1 / count;
	*t2 = *t2 / count;
	count = 0;
	}

void stddev(TString filename , double t1, double t2, double * dt1, double * dt2){
	//Takes a file with three columns and calculates the standard deviation of the last two returning them in variables *dt1 and *dt2.
	//The already calculated means t1 and t2 are used as an input.
	double trash;
	double count = 0;
	double temp1;
	double temp2;
	
	std::ifstream myfile;
	myfile.open(filename);

   while ( true ) {
     myfile >> trash;
     myfile >> temp1;
     myfile >> temp2;
     *dt1 += (temp1 - t1)*(temp1 - t1);
     *dt2 += (temp2 - t2)*(temp2 - t2);
     count ++;
     if( myfile.eof() ) break;
   };
      
   myfile.close();
	
	*dt1 = *dt1 / count / (count-1);
	*dt2 = *dt2 / count / (count-1);
	*dt1 = TMath::Sqrt(*dt1);
	*dt2 = TMath::Sqrt(*dt2);
	count = 0;
	}

int timecal(){
	// This program is used to calibrate the output of two TAC-ADC systems 
	double times[13];
	double dtimes[13];
	double t1[13] = {0};
	double dt1[13] = {0};
	double t2[13] = {0};
	double dt2[13] = {0};
	TString filename;
	
	for (int i = 0; i<13; i++){
		//cycle on the times used to calibrate from 10ns to 130ns 
		filename = "data/timecal/1G"+ToString(i+3)+"_timecal_"+ToString(10*(i+1))+"ns.dat"; //file to analyze 
		times[i] = 10*(i+1); //input time
		dtimes[i] = 1; //uncertainty on input time intervals
		averager(filename, &t1[i], &t2[i]);//get mean of the ADC distributions
		stddev(filename, t1[i], t2[i], &dt1[i], &dt2[i]);//get RMS of the ADC distributions
		cout << times[i] << " " << t1[i] << " " << t2[i] << endl;
	}
	
	TGraphErrors * cal1 = new TGraphErrors(13, times, t1, dtimes, dt1);//Create the graphs for the calibration
	TGraphErrors * cal2 = new TGraphErrors(13, times, t2, dtimes, dt2);
	
	cal1 -> SetTitle(" ");
	cal1 -> GetXaxis() -> SetTitle("time [ns]");
	cal1 -> GetYaxis() -> SetTitle("ADC value [au]");
	cal2 -> SetLineColor(2);
	cal1 -> Draw("AP");
	cal2 -> Draw("Psame");
	
	cal1 -> Fit("pol1", "", "R", 15, 115);//Fit to calibrate
	cal2 -> Fit("pol1", "", "R", 15, 115);
	
	return 0;
	}
