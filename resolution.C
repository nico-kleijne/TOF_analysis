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
#define sp_res 0.85


TString ToString(int num){
	ostringstream start;
	start<<num;
	TString start1=start.str();
	return start1;	
}

void time_plot(TString filename, TH1D * hist, double * mean, double * stddev, int column){
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
     if (column == 0) if(temp1/a1!=0) hist -> Fill(temp1/a1);//Convert from ADC to time.
     if (column == 1) if(temp2/a2!=0) hist -> Fill(temp2/a2);
     if( myfile.eof() ) break;
   };
      
   myfile.close();	
   
   *mean = hist -> GetMean();
   *stddev = hist -> GetRMS();
   //*stddev = hist -> GetRMS()/TMath::Sqrt(hist -> GetEntries()-1);
}


int resolution(){
	
	double mean[4][3][3]; // first index position (50, 100, 150, 200) , second index threshold (50, 100, 150) , third index scintillators (12, 13, 23);
	double stddev[4][3][3];
	double scin_res[4][3][3];// first index position (50, 100, 150, 200) , second index threshold (50, 100, 150) , third index scintillators (1, 2, 3);
	
	TString filename;
		
	TH1D * histos[4][3][3];
	TCanvas * canv[4][3][3];
	
	for (int i = 0; i < 4; i++){
		histos[i][1][0] = new TH1D("histo_"+ToString(i)+"_1_0", "histo_"+ToString(i)+"_1_0", 100, 0, 100);
		canv[i][1][0] = new TCanvas("c_"+ToString(i)+"_1_0", "c_"+ToString(i)+"_1_0", 1);
		filename = "data/resolutions/2L" + ToString(i+1) + "_fine_tcal12_" + ToString(50*(i+1)) + "_100.dat";
		time_plot(filename, histos[i][1][0], &mean[i][1][0], &stddev[i][1][0], 0);
		cout << "Pos: " << ToString(50*(i+1)) << " Thres: 100 Diff: 12 " << stddev[i][1][0] << endl;
		canv[i][1][0] -> cd();
		histos[i][1][0] -> Draw();
		}
		
	for (int i = 0; i < 4; i++){
		histos[i][1][1] = new TH1D("histo_"+ToString(i)+"_1_1", "histo_"+ToString(i)+"_1_1", 100, 0, 100);
		canv[i][1][1] = new TCanvas("c_"+ToString(i)+"_1_1", "c_"+ToString(i)+"_1_1", 1);
		filename = "data/resolutions/2L" + ToString(i+1) + "_fine_tcal12_" + ToString(50*(i+1)) + "_100.dat";
		time_plot(filename, histos[i][1][1], &mean[i][1][1], &stddev[i][1][1], 1);
		cout << "Pos: " << ToString(50*(i+1)) << " Thres: 100 Diff: 13 " << stddev[i][1][1] << endl;
		canv[i][1][1] -> cd();
		histos[i][1][1] -> Draw();
		}
		
	for (int i = 0; i < 4; i++){
		histos[i][0][0] = new TH1D("histo_"+ToString(i)+"_0_0", "histo_"+ToString(i)+"_0_0", 100, 0, 100);
		canv[i][0][0] = new TCanvas("c_"+ToString(i)+"_0_0", "c_"+ToString(i)+"_0_0", 1);
		filename = "data/resolutions/2L" + ToString(i+5) + "_fine_tcal12_" + ToString(50*(i+1)) + "_50.dat";
		time_plot(filename, histos[i][0][0], &mean[i][0][0], &stddev[i][0][0], 0);
		cout << "Pos: " << ToString(50*(i+1)) << " Thres: 50 Diff: 12 " << stddev[i][0][0] << endl;
		canv[i][0][0] -> cd();
		histos[i][0][0] -> Draw();
		}
		
	for (int i = 0; i < 4; i++){
		histos[i][0][1] = new TH1D("histo_"+ToString(i)+"_0_1", "histo_"+ToString(i)+"_0_1", 100, 0, 100);
		canv[i][0][1] = new TCanvas("c_"+ToString(i)+"_0_1", "c_"+ToString(i)+"_0_1", 1);
		filename = "data/resolutions/2L" + ToString(i+5) + "_fine_tcal12_" + ToString(50*(i+1)) + "_50.dat";
		time_plot(filename, histos[i][0][1], &mean[i][0][1], &stddev[i][0][1], 1);
		cout << "Pos: " << ToString(50*(i+1)) << " Thres: 50 Diff: 13 " << stddev[i][0][1] << endl;
		canv[i][0][1] -> cd();
		histos[i][0][1] -> Draw();
		}
		
	for (int i = 0; i < 4; i++){
		histos[i][2][0] = new TH1D("histo_"+ToString(i)+"_2_0", "histo_"+ToString(i)+"_2_0", 100, 0, 100);
		canv[i][2][0] = new TCanvas("c_"+ToString(i)+"_2_0", "c_"+ToString(i)+"_2_0", 1);
		filename = "data/resolutions/2M" + ToString(i+1) + "_fine_tcal12_" + ToString(50*(i+1)) + "_150.dat";
		time_plot(filename, histos[i][2][0], &mean[i][2][0], &stddev[i][2][0], 0);
		cout << "Pos: " << ToString(50*(i+1)) << " Thres: 150 Diff: 12 " << stddev[i][2][0] << endl;
		canv[i][2][0] -> cd();
		histos[i][2][0] -> Draw();
		}
		
	for (int i = 0; i < 4; i++){
		histos[i][2][1] = new TH1D("histo_"+ToString(i)+"_2_1", "histo_"+ToString(i)+"_2_1", 100, 0, 100);
		canv[i][2][1] = new TCanvas("c_"+ToString(i)+"_2_1", "c_"+ToString(i)+"_2_1", 1);
		filename = "data/resolutions/2M" + ToString(i+1) + "_fine_tcal12_" + ToString(50*(i+1)) + "_150.dat";
		time_plot(filename, histos[i][2][1], &mean[i][2][1], &stddev[i][2][1], 1);
		cout << "Pos: " << ToString(50*(i+1)) << " Thres: 150 Diff: 13 " << stddev[i][2][1] << endl;
		canv[i][2][1] -> cd();
		histos[i][2][1] -> Draw();
		}
		
	for (int i = 0; i < 4; i++){
		histos[i][2][2] = new TH1D("histo_"+ToString(i)+"_2_2", "histo_"+ToString(i)+"_2_2", 100, 0, 100);
		canv[i][2][2] = new TCanvas("c_"+ToString(i)+"_2_2", "c_"+ToString(i)+"_2_2", 1);
		filename = "data/resolutions/2M" + ToString(i+5) + "_fine_tcal23_" + ToString(50*(i+1)) + "_150.dat";
		time_plot(filename, histos[i][2][2], &mean[i][2][2], &stddev[i][2][2], 0);
		cout << "Pos: " << ToString(50*(i+1)) << " Thres: 150 Diff: 23 " << stddev[i][2][2] << endl;
		canv[i][2][2] -> cd();
		histos[i][2][2] -> Draw();
		}
		
	for (int i = 0; i < 2; i++){
		histos[i][1][2] = new TH1D("histo_"+ToString(i)+"_1_2", "histo_"+ToString(i)+"_1_2", 100, 0, 100);
		canv[i][1][2] = new TCanvas("c_"+ToString(i)+"_1_2", "c_"+ToString(i)+"_1_2", 1);
		filename = "data/resolutions/2M" + ToString(i+9) + "_fine_tcal23_" + ToString(50*(i+1)) + "_100.dat";
		time_plot(filename, histos[i][1][2], &mean[i][1][2], &stddev[i][1][2], 0);
		cout << "Pos: " << ToString(50*(i+1)) << " Thres: 100 Diff: 23 " << stddev[i][2][2] << endl;
		canv[i][1][2] -> cd();
		histos[i][1][2] -> Draw();
		}
	
	for (int i = 2; i < 4; i++){
		histos[i][1][2] = new TH1D("histo_"+ToString(i)+"_1_2", "histo_"+ToString(i)+"_1_2", 100, 0, 100);
		canv[i][1][2] = new TCanvas("c_"+ToString(i)+"_1_2", "c_"+ToString(i)+"_1_2", 1);
		filename = "data/resolutions/2G" + ToString(i-1) + "_fine_tcal23_" + ToString(50*(i+1)) + "_100.dat";
		time_plot(filename, histos[i][1][2], &mean[i][1][2], &stddev[i][1][2], 0);
		cout << "Pos: " << ToString(50*(i+1)) << " Thres: 100 Diff: 23 " << stddev[i][1][2] << endl;
		canv[i][1][2] -> cd();
		histos[i][1][2] -> Draw();
		}
		
	for (int i = 0; i < 4; i++){
		histos[i][0][2] = new TH1D("histo_"+ToString(i)+"_0_2", "histo_"+ToString(i)+"_0_2", 100, 0, 100);
		canv[i][0][2] = new TCanvas("c_"+ToString(i)+"_0_2", "c_"+ToString(i)+"_0_2", 1);
		filename = "data/resolutions/2G" + ToString(i+3) + "_fine_tcal23_" + ToString(50*(i+1)) + "_50.dat";
		time_plot(filename, histos[i][0][2], &mean[i][0][2], &stddev[i][0][2], 0);
		cout << "Pos: " << ToString(50*(i+1)) << " Thres: 50 Diff: 23 " << stddev[i][0][2] << endl;
		canv[i][0][2] -> cd();
		histos[i][0][2] -> Draw();
		}		
		
		cout << "____________________________________________" << endl;
		
		
		
	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 3; j++){
			scin_res[i][j][0] = TMath::Sqrt((stddev[i][j][0]*stddev[i][j][0]+stddev[i][j][1]*stddev[i][j][1]-stddev[i][j][2]*stddev[i][j][2]-sp_res*sp_res)/2);
			scin_res[i][j][1] = TMath::Sqrt((stddev[i][j][2]*stddev[i][j][2]+stddev[i][j][0]*stddev[i][j][0]-stddev[i][j][1]*stddev[i][j][1]-sp_res*sp_res)/2);
			scin_res[i][j][2] = TMath::Sqrt((stddev[i][j][1]*stddev[i][j][1]+stddev[i][j][2]*stddev[i][j][2]-stddev[i][j][0]*stddev[i][j][0]-sp_res*sp_res)/2);
			cout << "Pos: " << ToString(50*(i+1)) << " Thres: " << ToString(50*(j+1)) << " PMT1: " << scin_res[i][j][0] << endl;
			cout << "Pos: " << ToString(50*(i+1)) << " Thres: " << ToString(50*(j+1)) << " PMT2: " << scin_res[i][j][1] << endl;
			cout << "Pos: " << ToString(50*(i+1)) << " Thres: " << ToString(50*(j+1)) << " PMT3: " << scin_res[i][j][2] << endl;
		}
	}	
	return 0;
}
