#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <string>
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
#include <TRandom2.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TFitResultPtr.h>
#include <RConfig.h>
#include <TStopwatch.h>


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
#define Dcy 2
#define res_bin 100

TRandom3* Rx = new TRandom3(23);
//Velocità nella luce nel vuoto
Double_t C = 29.9792458 ;
//Dimensioni fisichedel sistema
//Velocità di propagazione della luce nella sbarra
Double_t beta_s = 0.4697 * C;

//Lunghezze in metri
//Lunghezze lastra sopra
Double_t Ux = 279;
Double_t Uy = 4;

//Lunghezze lastra sotto
Double_t Dx = 14;
Double_t Dy = 23;


//Qulache variabile per la simulazione
//Massa del muone in GeV
Double_t M_mu = 0.105;

//Intervallo di energia della simulazione in GeV
Double_t Emin = M_mu;
Double_t Emax = 1 ; 


//Variabile per simulare la lettura del primo TAC
Double_t c2=-0.35, s2=0.06, delay2= 30.5;

//Variabile per simulare la lettura del secondo TAC
Double_t c1=-0.32, s1=0.07, delay1= 30.5;

TString ToString(int num){
	ostringstream start;
	start<<num;
	TString start1=start.str();
	return start1;	
}

void res_adder(TString filename, TH1D * resolution13, TH1D * resolution23){
	
	double trash;
	double temp1;
	double temp2;
	double mean1 = 0;
	double mean2 = 0;
	int counter1 = 0;
	int counter2 = 0;
	
	std::ifstream myfile;
	myfile.open(filename);

   while ( true ) {
     myfile >> trash;
     myfile >> temp1;
     myfile >> temp2;
     if(temp1 != 0){
		mean1 += temp1;
		counter1 ++;		
	 }
	 if(temp2 != 0){
		mean2 += temp2;
		counter2 ++;		
	 }
     if( myfile.eof() ) break;
   };
   
   myfile.close();
	
	mean1 = mean1 / counter1;
	mean2 = mean2 / counter2;
	
	
	myfile.open(filename); 
	
	while ( true ) {
     myfile >> trash;
     myfile >> temp1;
     myfile >> temp2;
     if(temp1 != 0){
		resolution23 -> Fill((temp1 - mean1)/a1);		
	 }
	 if(temp2 != 0){
		resolution13 -> Fill((temp2 - mean2)/a2);
	 }
     if( myfile.eof() ) break;
   };
	
	myfile.close();
	
	}

double random_extracter(TH1D * resolution, double distribution[500]){
	
	int partial_integral;
	int partial_sum = 0;
	double result = 0;
	partial_integral = Rx -> Integer(resolution -> GetEntries());
	for (int i = 0; i<500; i++){
		partial_sum += distribution[i];
		if (partial_sum >= partial_integral + 1){
			result = (i-250)*0.2-0.1;
			//cout << "partial_sum " << partial_sum << endl;
			//cout << partial_integral << endl;
			break;
			}			
	//cout << i << endl;
	}
	//cout << result << endl;
	return result;
}

int TOF(){

	double trash;
	double temp1;
	double temp2;
	double tof;
	double dist;
	double distribution13[500];
	double distribution23[500];
	double Cusu;
	int count_test = 0;
	
	unsigned i=0, N = 4.35e6;
	//Counter vari
	unsigned Missed =0;
	

	//Variabili dell'evento
	Double_t Xu=0, Yu=0, Xd=0, Yd=0, C_Theta=0, Phi=0;
	//Variabile di un evento andato a segno
	Double_t E=0, beta_mu=0, Ts=0, Td=0, Td2=0, Vs=0, Vd=0;
	
	//Variabile di ricostruzione
	Double_t Tsr=0, Tdr =0, Tdr2, Xur=0, Xur2=0, beta_mur=0;
	
	//TString filename = "data/distribution/1V6_barcal3.dat";
	//TString filename = "data/resolutions/2L1_fine_tcal12_50_100.dat";
	TString filename = "data/2G7_beta1_mis1.dat";

	TH1D * hist_tof = new TH1D("hist_tof", "hist_tof", 70, -20, 50);
	TH1D * hist_tof_cen = new TH1D("hist_tof_cen", "hist_tof_cen", 70, -20, 50);
	TH1D * hist_dist = new TH1D("hist_dist", "hist_dis", 100, -100, 400);
	TH1D * hist_beta = new TH1D("hist_beta", "hist_beta", 300, 0, 10);
	TH2D * hist_dist_tof = new TH2D("hist_dist_tof", "hist_dist_tof", 120, -30, 330, 50, -5, 20);	
	TH1F * histo_beta = new TH1F("histo_beta","histo_beta", 300, 0, 10);
	TH1F * histo_TOF = new TH1F("histoTOF","histoTOF", 70, -20, 50);
	TH1F * histo_x = new TH1F("histox","histox", 100, -100, 400);
	TH1D * resolution13 = new TH1D("resolution13","resolution13", 500, -50, 50);
	TH1D * resolution23 = new TH1D("resolution23","resolution23", 500, -50, 50);
	TH1D * rand_dist = new TH1D("rand_dist","rand_dist", 500, -50, 50);
	//TH1D * his_tdr = new TH1F("hist_tdr","hist_tdr", 100, -50, 50);

	
	std::ifstream myfile;
	myfile.open(filename);

   while ( true ) {
     myfile >> trash;
     myfile >> temp1;
     myfile >> temp2;
     if(temp2 != 0){
		//tof = temp2/a2 - temp1/2/a1 - of2/a2 + of1/2/a1 + of3/2 - 15;
		//tof = temp2/a2 - temp1/2/a1 - 2.2665;
		tof = temp2/2/a2 - temp1/2/a1 + 6.99;
		//dist = temp1/a1/a3 - of1/a1/a3 -of3/a3 -30/a3;
		dist = -(temp1/a1+temp2/a2)/a3-203;
		hist_dist -> Fill(dist);
		hist_dist_tof -> Fill(dist, tof);
		if ( dist > 250 && dist < 500) hist_tof_cen -> Fill(tof);
		dist = TMath::Sqrt((dist-140)*(dist-140)+Zu*Zu);
		hist_tof -> Fill(tof); 
		hist_beta -> Fill(dist/tof/30);
	 }
     if( myfile.eof() ) break;
   };
   
	res_adder("data/resolutions/2M9_fine_tcal23_50_100.dat", resolution13, resolution23);
	res_adder("data/resolutions/2M10_fine_tcal23_100_100.dat", resolution13, resolution23);
	res_adder("data/resolutions/2G1_fine_tcal23_150_100.dat", resolution13, resolution23);
	res_adder("data/resolutions/2G2_fine_tcal23_200_100.dat", resolution13, resolution23);
	
	TCanvas * c1_prov = new TCanvas("c1_prov", "c1_prov", 1);
	c1_prov -> cd();
	resolution13 -> Draw();
	
	TCanvas * c2_prov = new TCanvas("c2_prov", "c2_prov", 1);
	c2_prov -> cd();
	resolution23 -> Draw();
	
	for(int i = 0; i < 500; i++){
		distribution13[i] = resolution13 -> GetBinContent(i);
		distribution23[i] = resolution23 -> GetBinContent(i);
		//count_test += distribution23[i];
		//cout << distribution23[i] << endl;
	}
	
	//cout << count_test << endl;
   
	for (i=0; i<N; i++){
			
			//Estraggo posizione e direzione iniziale
			Xu = (Rx -> Rndm()) * Ux;
			Yu = (Rx -> Rndm()) * Uy;
			C_Theta = TMath::Power(Rx -> Rndm(),1.0/3.0);
			Phi = (Rx -> Rndm()) * 2 * (TMath::Pi());
			
			//Calcolo posizione finale
			Xd = Xu - Zu * TMath::Cos(Phi) * sqrt(1 - pow(C_Theta, 2.0))/C_Theta;
			Yd = Yu - Zu * TMath::Sin(Phi) * sqrt(1 - pow(C_Theta, 2.0))/C_Theta;
			if(Dcx - Dx/2 <= Xd && Xd<= Dcx + Dx/2 && Dcy - Dy/2 <= Yd && Yd <= Dcy + Dy/2){
				//Genero l'energia e la velocità
				E = (Rx -> Rndm()) *( Emax - Emin) + Emin;
				//beta_mu = sqrt (1 - pow(M_mu/E, 2.0)) * C;
				beta_mu = C;
				
				//Genero i tempi in lettura nella barra
				Ts = (Ux - 2 * Xu)/beta_s + delay1;
				//Vs = Rx -> Gaus(a1 * Ts + c1,s1);
				
				//Genero i TOF
				Td = sqrt(pow(Zu,2.0) + pow( Xu-Xd, 2.0) + pow( Yu-Yd, 2.0))/beta_mu + delay2 - Xu/beta_s;
				Td2 = sqrt(pow(Zu,2.0) + pow( Xu-Xd, 2.0) + pow( Yu-Yd, 2.0))/beta_mu + delay1 + Xu/beta_s;
				//Vd = Rx -> Gaus(a2 * Td + c2,s2);
				
				//Ricostruisco
				//Tdr = Td;
				//Tdr = Td + Rx -> Gaus(0, 4.5);
				Tdr = Td + random_extracter(resolution13, distribution13);
				//Tdr2 = Td2 + Rx -> Gaus(0, 4.5);
				Tdr2 = Td2 + random_extracter(resolution23, distribution23);
				//Tsr = Ts;
				//Tsr = Ts + random_extracter(resolution13, distribution13);
				Tsr = Ts;
				Xur = (Ux-(Tsr - delay1)*beta_s)/2;
				Xur2 = (Tdr2 - Tdr - delay1 +delay2)*beta_s/2;
				beta_mur = sqrt(pow(Zu,2.0) + pow( Xur-Dcx, 2.0))/(Tdr -delay2 + Xur2/beta_s);
				
				histo_beta -> Fill(beta_mur/C);
				histo_TOF -> Fill((Tdr+Tdr2)/2-delay1);
				histo_x -> Fill(Xur2);
				// Stampo le variabili generate e la velcità reale
				//fout  << Vs << '\t'  << Vd << '\t'  << beta_mu << endl;
			}
			else {Missed ++;}
	}

	for (int j = 0; j < 10000; j++){
		Cusu = random_extracter(resolution23, distribution23);
		rand_dist -> Fill(Cusu);
	}


	TCanvas * c_rand_dist = new TCanvas("c_rand_dist", "c_rand_dist", 1);
	c_rand_dist -> cd();
	rand_dist -> Draw();
	
	//TCanvas * ctdr = new TCanvas("ctdr", "ctdr", 1);
	//ctdr -> cd();
	//hist_tdr -> Draw();
   
   TCanvas * c1 = new TCanvas("c1", "c1", 1);
   c1 -> cd();
   hist_tof -> GetXaxis() -> SetTitle("Time of flight [ns]");
   hist_tof -> GetYaxis() -> SetTitle("Number of events");
   histo_TOF -> SetLineColor(2);
   histo_TOF -> Scale(hist_tof -> GetEntries() / histo_TOF -> GetEntries());
   histo_TOF -> Draw("hist");
   hist_tof -> Draw("same");
   cout << hist_tof -> GetMean() << endl;
   cout << histo_TOF -> GetMean() << endl;
   TCanvas * c2 = new TCanvas("c2", "c2", 1);
   c2 -> cd();
   hist_dist -> GetXaxis() -> SetTitle("Position on the scintillator [cm]");
   hist_dist -> GetYaxis() -> SetTitle("Number of events");
   hist_dist -> Draw();
   histo_x -> SetLineColor(2);
   histo_x -> Scale(hist_dist -> GetEntries() / histo_x -> GetEntries());
   histo_x -> Draw("samehist");
   TCanvas * c3 = new TCanvas("c3", "c3", 1);
   c3 -> cd();
   hist_beta -> GetXaxis() -> SetTitle("#beta_{#mu}");
   hist_beta -> GetYaxis() -> SetTitle("Number of events");
   histo_beta -> SetLineColor(2);
   histo_beta -> Scale(hist_beta -> GetEntries() / histo_beta -> GetEntries());
   histo_beta -> Draw("hist");
   hist_beta -> Draw("same");
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
