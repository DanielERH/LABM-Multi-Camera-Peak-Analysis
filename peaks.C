//
// Daniel R. V3.0
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
 #include "TSpectrum2.h"
 #include "TCanvas.h"
 #include "TRandom.h"
 #include "TH2.h"
 #include "TF2.h"
 #include "TMath.h"
 #include "TROOT.h"
#include "TSystem.h"
#include <chrono>

#define maxpeaks 400
#define maxeventnum 30000

string pdfname;
string imagenameanddir;
string pngnames[maxeventnum][8];

int EventCount=0, CameraCount=0, camMaxEvent[8][5]={0}, frameindex=0;
int camminnp[8]={0}, cammaxnp[8]={0}; 
double camMaxRates[8][5]={0}, EventData[maxeventnum][8][650][3]={0};
float MinTimeStamp=0, MaxTimeStamp=0, MinHCurrent=0, MaxHCurrent=0, MinLCurrent=0, MaxLCurrent=0, MinHLumi=0, MaxHLumi=0, MinLLumi=0, MaxLLumi=0, MinECLLumi=0, MaxECLLumi=0, MinHMQC2LEX=0, MaxHMQC2LEX=0, MinHMQC1LEX=0, MaxHMQC1LEX=0, MinHMQC1REX=0, MaxHMQC1REX=0, MinHMQC2REX=0, MaxHMQC2REX=0, MinHMQC2LEY=0, MaxHMQC2LEY=0, MinHMQC1LEY=0, MaxHMQC1LEY=0, MinHMQC1REY=0, MaxHMQC1REY=0, MinHMQC2REY=0, MaxHMQC2REY=0, MinLMQC2RPX=0, MaxLMQC2RPX=0, MinLMQC1RPX=0, MaxLMQC1RPX=0, MinLMQC1LPX=0, MaxLMQC1LPX=0, MinLMQC2LPX=0, MaxLMQC2LPX=0, MinLMQC2RPY=0, MaxLMQC2RPY=0, MinLMQC1RPY=0, MaxLMQC1RPY=0, MinLMQC1LPY=0, MaxLMQC1LPY=0, MinLMQC2LPY=0, MaxLMQC2LPY=0, dummy=0;
float Allpvs[maxeventnum][22]={0};

char imagename[50], n[200];

TCanvas *c1 = new TCanvas("c1","Image Heat Map");

TSpectrum2 *s;
s = new TSpectrum2(2*maxpeaks);

//Histograms
#define NHISTS 65
int ihist=1;
TH1F *camPVhists[NHISTS];

//Functions
void readevent(int num, int div);
void analyzeimage();
void FindMaxPeaks();
void MaxRGBHM(int evnt);
void produceHM(int eventnm, int camnum);
void hists(int eventsnum);
void createNtuple(int evntsnm);

//Directories
//Main directory
string maindir="runs/run6/";

//Camera png images directory
string pngdir= maindir + string("pngPics/");

//Camera png images heatmap generated save directory
string heatmdir= maindir + string("images/");  

//Histograms graphs save directory
string histdir= maindir + string("graphs/"); 

//Data File directory
sprintf(n, "runs/run6/data-run6-all-20240502-Example");
fstream labmdata (n);

//NTuple File directory
string rootfilenamedummy = maindir + string("ntuple/run6.root");
string treenamedummy = string("run6");

const char *rootfilename = rootfilenamedummy.c_str();
const char *treename = treenamedummy.c_str();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void peaks(){

   do{
      EventCount++;
      readevent(EventCount, 1); //Read each event (if number is 1, otherwhise every "x" event) to save number of peaks and rates per peak for each camera
      if (labmdata.peek()==EOF) break; //we peek the next character and if is the end of file, then we break the loop
   }while (true);
   labmdata.close();

   FindMaxPeaks();

   for(int i=1; i<=EventCount; i++){
      MaxRGBHM(i);
   }
   
   hists(EventCount);
   
   createNtuple(EventCount);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void readevent(int num, int div){ 
   const char* frames[] = {"\\", "|", "/", "-"};
   int pvscontup=0;
   
   labmdata >> dummy;
   if(num==1){
      MaxTimeStamp=dummy;
      MinTimeStamp=dummy;
   }else{
      if (dummy>MaxTimeStamp) MaxTimeStamp=dummy;
      if (dummy<MinTimeStamp) MinTimeStamp=dummy;
   }
   Allpvs[num][pvscontup]=dummy; //TimeStamp; pvscontup=0
   pvscontup++;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   if(num==1){
      MaxHCurrent=dummy;
      MinHCurrent=dummy;
   }else{
      if (dummy>MaxHCurrent) MaxHCurrent=dummy;
      if (dummy<MinHCurrent) MinHCurrent=dummy;
   }
   Allpvs[num][pvscontup]=dummy; //HCurrent; pvscontup=1
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxLCurrent=dummy;
      MinLCurrent=dummy;
   }else{
      if (dummy>MaxLCurrent) MaxLCurrent=dummy;
      if (dummy<MinLCurrent) MinLCurrent=dummy;
   }
   Allpvs[num][pvscontup]=dummy; //LCurrent; pvscontup=2
   pvscontup++;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   if(num==1){
      MaxHLumi=dummy;
      MinHLumi=dummy;
   }else{
      if (dummy>MaxHLumi) MaxHLumi=dummy;
      if (dummy<MinHLumi) MinHLumi=dummy;
   }
   Allpvs[num][pvscontup]=dummy; //HLumi; pvscontup=3
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxLLumi=dummy;
      MinLLumi=dummy;
   }else{
      if (dummy>MaxLLumi) MaxLLumi=dummy;
      if (dummy<MinLLumi) MinLLumi=dummy;
   }	
   Allpvs[num][pvscontup]=dummy; //LLumi; pvscontup=4
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxECLLumi=dummy;
      MinECLLumi=dummy;
   }else{
      if (dummy>MaxECLLumi) MaxECLLumi=dummy;
      if (dummy<MinECLLumi) MinECLLumi=dummy;
   }	
   Allpvs[num][pvscontup]=dummy; //ECLLumi; pvscontup=5
   pvscontup++;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   if(num==1){
      MaxHMQC2LEX=dummy;
      MinHMQC2LEX=dummy;
   }else{
      if (dummy>MaxHMQC2LEX) MaxHMQC2LEX=dummy;
      if (dummy<MinHMQC2LEX) MinHMQC2LEX=dummy;
   }
   Allpvs[num][pvscontup]=dummy; //HMQC2LEX; pvscontup=6
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxHMQC1LEX=dummy;
      MinHMQC1LEX=dummy;
   }else{
      if (dummy>MaxHMQC1LEX) MaxHMQC1LEX=dummy;
      if (dummy<MinHMQC1LEX) MinHMQC1LEX=dummy;
   }	
   Allpvs[num][pvscontup]=dummy; //HMQC1LEX; pvscontup=7
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxHMQC1REX=dummy;
      MinHMQC1REX=dummy;
   }else{
      if (dummy>MaxHMQC1REX) MaxHMQC1REX=dummy;
      if (dummy<MinHMQC1REX) MinHMQC1REX=dummy;
   }		
   Allpvs[num][pvscontup]=dummy; //HMQC1REX; pvscontup=8
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxHMQC2REX=dummy;
      MinHMQC2REX=dummy;
   }else{
      if (dummy>MaxHMQC2REX) MaxHMQC2REX=dummy;
      if (dummy<MinHMQC2REX) MinHMQC2REX=dummy;
   }	
   Allpvs[num][pvscontup]=dummy; //HMQC2REX; pvscontup=9
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxHMQC2LEY=dummy;
      MinHMQC2LEY=dummy;
   }else{
      if (dummy>MaxHMQC2LEY) MaxHMQC2LEY=dummy;
      if (dummy<MinHMQC2LEY) MinHMQC2LEY=dummy;
   }	
   Allpvs[num][pvscontup]=dummy; //HMQC2LEY; pvscontup=10
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxHMQC1LEY=dummy;
      MinHMQC1LEY=dummy;
   }else{
      if (dummy>MaxHMQC1LEY) MaxHMQC1LEY=dummy;
      if (dummy<MinHMQC1LEY) MinHMQC1LEY=dummy;
   }	
   Allpvs[num][pvscontup]=dummy; //HMQC1LEY; pvscontup=11
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxHMQC1REY=dummy;
      MinHMQC1REY=dummy;
   }else{
      if (dummy>MaxHMQC1REY) MaxHMQC1REY=dummy;
      if (dummy<MinHMQC1REY) MinHMQC1REY=dummy;
   }
   Allpvs[num][pvscontup]=dummy; //HMQC1REY; pvscontup=12
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxHMQC2REY=dummy;
      MinHMQC2REY=dummy;
   }else{
      if (dummy>MaxHMQC2REY) MaxHMQC2REY=dummy;
      if (dummy<MinHMQC2REY) MinHMQC2REY=dummy;
   }	
   Allpvs[num][pvscontup]=dummy; //HMQC2REY; pvscontup=13
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxLMQC2RPX=dummy;
      MinLMQC2RPX=dummy;
   }else{
      if (dummy>MaxLMQC2RPX) MaxLMQC2RPX=dummy;
      if (dummy<MinLMQC2RPX) MinLMQC2RPX=dummy;
   }		
   Allpvs[num][pvscontup]=dummy; //LMQC2RPX; pvscontup=14
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxLMQC1RPX=dummy;
      MinLMQC1RPX=dummy;
   }else{
      if (dummy>MaxLMQC1RPX) MaxLMQC1RPX=dummy;
      if (dummy<MinLMQC1RPX) MinLMQC1RPX=dummy;
   }	
   Allpvs[num][pvscontup]=dummy; //LMQC1RPX; pvscontup=15
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxLMQC1LPX=dummy;
      MinLMQC1LPX=dummy;
   }else{
      if (dummy>MaxLMQC1LPX) MaxLMQC1LPX=dummy;
      if (dummy<MinLMQC1LPX) MinLMQC1LPX=dummy;
   }		
   Allpvs[num][pvscontup]=dummy; //LMQC1LPX; pvscontup=16
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxLMQC2LPX=dummy;
      MinLMQC2LPX=dummy;
   }else{
      if (dummy>MaxLMQC2LPX) MaxLMQC2LPX=dummy;
      if (dummy<MinLMQC2LPX) MinLMQC2LPX=dummy;
   }	
   Allpvs[num][pvscontup]=dummy; //LMQC2LPX; pvscontup=17
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxLMQC2RPY=dummy;
      MinLMQC2RPY=dummy;
   }else{
      if (dummy>MaxLMQC2RPY) MaxLMQC2RPY=dummy;
      if (dummy<MinLMQC2RPY) MinLMQC2RPY=dummy;
   }	
   Allpvs[num][pvscontup]=dummy; //LMQC2RPY; pvscontup=18
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxLMQC1RPY=dummy;
      MinLMQC1RPY=dummy;
   }else{
      if (dummy>MaxLMQC1RPY) MaxLMQC1RPY=dummy;
      if (dummy<MinLMQC1RPY) MinLMQC1RPY=dummy;
   }	
   Allpvs[num][pvscontup]=dummy; //LMQC1RPY; pvscontup=19
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxLMQC1LPY=dummy;
      MinLMQC1LPY=dummy;
   }else{
      if (dummy>MaxLMQC1LPY) MaxLMQC1LPY=dummy;
      if (dummy<MinLMQC1LPY) MinLMQC1LPY=dummy;
   }
   Allpvs[num][pvscontup]=dummy; //LMQC1LPY; pvscontup=20
   pvscontup++;
   labmdata >> dummy;
   if(num==1){
      MaxLMQC2LPY=dummy;
      MinLMQC2LPY=dummy;
   }else{
      if (dummy>MaxLMQC2LPY) MaxLMQC2LPY=dummy;
      if (dummy<MinLMQC2LPY) MinLMQC2LPY=dummy;
   }
   Allpvs[num][pvscontup]=dummy; //LMQC2LPY; pvscontup=21
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata >> dummy;
   labmdata.ignore(1, '\n');

   if (isdigit(labmdata.peek())==0){ //if the next character is NOT a digit then we proceed
      do{		
         labmdata >> imagename; //save the png image name
         
         if(imagename[6]==string("1"))CameraCount=0; //Check what camera we are working on right now
         if(imagename[6]==string("2"))CameraCount=1;
         if(imagename[6]==string("3"))CameraCount=2;
         if(imagename[6]==string("4"))CameraCount=3;
         if(imagename[6]==string("5"))CameraCount=4;
         if(imagename[6]==string("6"))CameraCount=5;
         if(imagename[6]==string("7"))CameraCount=6;
         if(imagename[6]==string("8"))CameraCount=7;
         
         pngnames[num][CameraCount]=string(imagename);
     
         std::cout << "\rLoading " << frames[frameindex] << std::flush;
         if(frameindex==3){
            frameindex=0;
         }else{
            frameindex++;
         }     

         if(num%div==0){
            imagenameanddir = pngdir + imagename;
            analyzeimage();
         }
         labmdata.ignore(1, '\n');
         if (labmdata.peek()==EOF) break; //after the last png name saved, we skip to the next character and if is the end of file, then we break from the camera loop.
      }while(isdigit(labmdata.peek())==0); //after the last png name saved, if the next character is NOT a digit, then we continue in the camera loop.
   } 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void analyzeimage(){
   // Process the camera png image to save number of peaks and rates of each peak per camera, per event
   
   int xbins=0, ybins=0, binpks=0;
   
   extern std::string imagenameanddir;

   TASImage* image = new TASImage(imagenameanddir.c_str()); //open the png image using TASImage

   //Get the png image x and y pixels, and the argb array
   UInt_t yPixels = image->GetHeight();
   UInt_t xPixels = image->GetWidth();
   UInt_t *argb   = image->GetArgbArray();
   
   TH2D *ImageHeatMap = new TH2D("ImageHeatMap", " ", xPixels, 0, xPixels, yPixels, 0, yPixels); //create the png image heat map histogram
   
   //Fill out the png heat map histogram
   for(int row=0; row<xPixels; ++row) {
      for(int col=0; col<yPixels; ++col) {
         UInt_t raw = argb[col*xPixels+row];
         float red = float((raw >> 16) & 0xff)/256;
         float green = float((raw >> 8) & 0xff)/256;
         float blue = float((raw >> 0) & 0xff)/256;

         ImageHeatMap->SetBinContent(row+1, yPixels-col, red+green+blue);
      }
   }
			
   ImageHeatMap->Rebin2D(60,60);
			
   //Use TSpectrum2 to search for peaks on the heatmap

   Int_t nfound = s->Search(ImageHeatMap,1,"colz,noMarkov",0.05);
   Double_t *xpeaks = s->GetPositionX(); //array of x position of each peak found
   Double_t *ypeaks = s->GetPositionY(); //array of y position of each peak found
   
   EventData[EventCount][CameraCount][0][0]=nfound; //Save number of peaks
   for (int i=1; i<=nfound; i++){
      EventData[EventCount][CameraCount][i][0]=ImageHeatMap->GetBinContent(ImageHeatMap->GetXaxis()->FindBin(xpeaks[i]),ImageHeatMap->GetXaxis()->FindBin(ypeaks[i])); //Save rate for each peak
      EventData[EventCount][CameraCount][i][1]=xpeaks[i]; //Save rate for each peak
      EventData[EventCount][CameraCount][i][2]=ypeaks[i]; //Save rate for each peak
   }
 
   delete ImageHeatMap;
   delete image;			
 }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FindMaxPeaks(){
   //Save the 5 highest rates and their events per camera
   for (int k=0; k<8; k++){
      for (int i=1; i<=EventCount; i++){
         if (EventData[i][k][0][0]>=2){
            for (int j=1; j<=EventData[i][k][0][0]; j++){
            
               if (camMaxRates[k][0]<=EventData[i][k][j][0]) {
                  camMaxRates[k][4]=camMaxRates[k][3];
                  camMaxRates[k][3]=camMaxRates[k][2];
                  camMaxRates[k][2]=camMaxRates[k][1];
                  camMaxRates[k][1]=camMaxRates[k][0];
                  camMaxRates[k][0]=EventData[i][k][j][0];
                  
                  camMaxEvent[k][4]=camMaxEvent[k][3];
                  camMaxEvent[k][3]=camMaxEvent[k][2];
                  camMaxEvent[k][2]=camMaxEvent[k][1];
                  camMaxEvent[k][1]=camMaxEvent[k][0];
                  camMaxEvent[k][0]=i;
               } else if (camMaxRates[k][1]<=EventData[i][k][j][0]) {
                  camMaxRates[k][4]=camMaxRates[k][3];
                  camMaxRates[k][3]=camMaxRates[k][2];
                  camMaxRates[k][2]=camMaxRates[k][1];
                  camMaxRates[k][1]=EventData[i][k][j][0];

                  camMaxEvent[k][4]=camMaxEvent[k][3];
                  camMaxEvent[k][3]=camMaxEvent[k][2];
                  camMaxEvent[k][2]=camMaxEvent[k][1];
                  camMaxEvent[k][1]=i;
               } else if (camMaxRates[k][2]<=EventData[i][k][j][0]) {
                  camMaxRates[k][4]=camMaxRates[k][3];
                  camMaxRates[k][3]=camMaxRates[k][2];
                  camMaxRates[k][2]=EventData[i][k][j][0];

                  camMaxEvent[k][4]=camMaxEvent[k][3];
                  camMaxEvent[k][3]=camMaxEvent[k][2];
                  camMaxEvent[k][2]=i;
               } else if (camMaxRates[k][3]<=EventData[i][k][j][0]) {
                  camMaxRates[k][4]=camMaxRates[k][3];
                  camMaxRates[k][3]=EventData[i][k][j][0];

                  camMaxEvent[k][4]=camMaxEvent[k][3];
                  camMaxEvent[k][3]=i;
               } else if (camMaxRates[k][4]<EventData[i][k][j][0]) {
                  camMaxRates[k][4]=EventData[i][k][j][0];
                  
                  camMaxEvent[k][4]=i;
               }
               
            }
         }
      }
   }		
 }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MaxRGBHM(int evnt){ 
   
   int cameracounter=0;
      
   for(int k=0; k<8; k++){ //for each camera
      for(int j=0; j<5; j++){ //Check all the Max rate events for each camera	
         if((evnt==camMaxEvent[k][j])and(cameracounter==0)){ //If the camera Max rate event and the event we are working matches, then print heat map; 
            imagenameanddir = pngdir + pngnames[evnt][k].c_str();
            produceHM(evnt, k);
            cameracounter=1; 
         }
      }
      cameracounter=0;
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void produceHM(int eventnm, int camnum){
   // Process the camera png image to create heatmaps with peaks highlighted

   extern std::string imagenameanddir;

   TASImage* image = new TASImage(imagenameanddir.c_str()); //open the png image using TASImage
			
   //Get the png image x and y pixels, and the argb array
   UInt_t yPixels = image->GetHeight();
   UInt_t xPixels = image->GetWidth();
   UInt_t *argb   = image->GetArgbArray();
   
   TH2D *ImageHeatMap = new TH2D("ImageHeatMap", " ", xPixels, 0, xPixels, yPixels, 0, yPixels); //create the png image heat map histogram
   
   //Fill out the png heat map histogram
   for(int row=0; row<xPixels; ++row) {
      for(int col=0; col<yPixels; ++col) {
         UInt_t raw = argb[col*xPixels+row];
         float red = float((raw >> 16) & 0xff)/256;
         float green = float((raw >> 8) & 0xff)/256;
         float blue = float((raw >> 0) & 0xff)/256;

         ImageHeatMap->SetBinContent(row+1, yPixels-col, red+green+blue);
      }
   }
			
   ImageHeatMap->Rebin2D(60,60);
			
   //Use TSpectrum2 to search for peaks on the heatmap

   Int_t nfound = s->Search(ImageHeatMap,1,"colz,noMarkov",0.05);
   Double_t *xpeaks = s->GetPositionX(); //array of x position of each peak found
   Double_t *ypeaks = s->GetPositionY(); //array of y position of each peak found
			
   auto m = new TMarker();
   m->SetMarkerStyle(23);
   m->SetMarkerColor(kRed);
   gPad->SetLogz();
			
   c1->Clear();
   pdfname = heatmdir + pngnames[eventnm][camnum].c_str() + std::string(".pdf(");
   ImageHeatMap->Draw("CONTZ");
   c1->Update();
			
   for (int i=0;i<nfound;i++) {
      m->DrawMarker(xpeaks[i],ypeaks[i]); //Set red triangular markers in each peak found
   }
			
   ImageHeatMap->SetTitle("Image Heat Map");
   c1->Update();
   c1->SaveAs(pdfname.c_str());
   c1->Clear();   
		
   pdfname = heatmdir + pngnames[eventnm][camnum].c_str() + std::string(".pdf)");
   ImageHeatMap->Draw("SURFZ");
   c1->Update();
   ImageHeatMap->SetTitle("Image Surf");
   c1->SaveAs(pdfname.c_str());
   c1->Clear();
   
   delete ImageHeatMap;
   delete image;			
 }
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hists(int eventsnum) {
   //Fill and draw histograms for the PVs
   
   int histcount=0, histcheck1=0, histcheck2=0;
   double maxrgb=0, max2rgb=0, minrgb=0, maxX=0, max2X=0, maxY=0, max2Y=0, minX=0, minY=0, calc=0;
   
   for(int i=1; i<=eventsnum; i++){
      if(i==1){
         camPVhists[ihist] = new TH1F("TimeStamp","Time Stamp", 100, MinTimeStamp, MaxTimeStamp); // ihsit = 1;
            camPVhists[ihist]->SetTitle("Time Stamp of each entry");
            camPVhists[ihist]->SetXTitle("Time Stamp");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
   	
         camPVhists[ihist] = new TH1F("HCurrent","HER Current", 100, MinHCurrent, MaxHCurrent); // ihsit = 2;
            camPVhists[ihist]->SetTitle("HER Current");
            camPVhists[ihist]->SetXTitle("HER Current (A)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
    	
         camPVhists[ihist] = new TH1F("LCurrent","LER Current", 100, MinLCurrent, MaxLCurrent); // ihsit = 3;
            camPVhists[ihist]->SetTitle("LER Current");
            camPVhists[ihist]->SetXTitle("LER Current (A)");   
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
    	
         camPVhists[ihist] = new TH1F("HLumi","HER Luminosity", 100, MinHLumi, MaxHLumi); // ihsit = 4;
            camPVhists[ihist]->SetTitle("HER Luminosity");
            camPVhists[ihist]->SetXTitle("HER Luminosity (cm^{-2}s^{-1})");
            camPVhists[ihist]->SetYTitle("dN"); 
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
    	
         camPVhists[ihist] = new TH1F("LLumi","LER Luminosity", 100, MinLLumi, MaxLLumi); // ihsit = 5;
            camPVhists[ihist]->SetTitle("LER Luminosity");
            camPVhists[ihist]->SetXTitle("LER Luminosity (cm^{-2}s^{-1})");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("ECLLumi","ECL Luminosity", 100, MinECLLumi, MaxECLLumi); // ihsit = 6;
            camPVhists[ihist]->SetTitle("ECL Luminosity");
            camPVhists[ihist]->SetXTitle("ECL Luminosity (cm^{-2}s^{-1})");
            camPVhists[ihist]->SetYTitle("dN");   
            camPVhists[ihist]->SetOption("HIST"); 
            camPVhists[ihist]->SetStats(0);  
            ++ihist;
    	
         camPVhists[ihist] = new TH1F("HMQC2LEX","HER MQC2LE x Position", 100, MinHMQC2LEX, MaxHMQC2LEX); // ihsit = 7;
            camPVhists[ihist]->SetTitle("HER MQC2LE x Position");
            camPVhists[ihist]->SetXTitle("HER MQC2LE x (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("HMQC1LEX","HER MQC1LE x Position", 100, MinHMQC1LEX, MaxHMQC1LEX); // ihsit = 8;
            camPVhists[ihist]->SetTitle("HER MQC1LE x Position");
            camPVhists[ihist]->SetXTitle("HER MQC1LE x (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("HMQC1REX","HER MQC1RE x Position", 100, MinHMQC1REX, MaxHMQC1REX); // ihsit = 9;
            camPVhists[ihist]->SetTitle("HER MQC1RE x Position");
            camPVhists[ihist]->SetXTitle("HER MQC1RE x (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("HMQC2REX","HER MQC2RE x Position", 100, MinHMQC2REX, MaxHMQC2REX); // ihsit = 10;
            camPVhists[ihist]->SetTitle("HER MQC2RE x Position");
            camPVhists[ihist]->SetXTitle("HER MQC2RE x (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("HMQC2LEY","HER MQC2LE y Position", 100, MinHMQC2LEY, MaxHMQC2LEY); // ihsit = 11;
            camPVhists[ihist]->SetTitle("HER MQC2LE y Position");
            camPVhists[ihist]->SetXTitle("HER MQC2LE y (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("HMQC1LEY","HER MQC1LE y Position", 100, MinHMQC1LEY, MaxHMQC1LEY); // ihsit = 12;
            camPVhists[ihist]->SetTitle("HER MQC1LE y Position");
            camPVhists[ihist]->SetXTitle("HER MQC1LE y (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("HMQC1REY","HER MQC1RE y Position", 100, MinHMQC1REY, MaxHMQC1REY); // ihsit = 13;
            camPVhists[ihist]->SetTitle("HER MQC1RE y Position");
            camPVhists[ihist]->SetXTitle("HER MQC1RE y (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("HMQC2REY","HER MQC2RE y Position", 100, MinHMQC2REY, MaxHMQC2REY); // ihsit = 14;
            camPVhists[ihist]->SetTitle("HER MQC2RE y Position");
            camPVhists[ihist]->SetXTitle("HER MQC2RE y (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
    	
         camPVhists[ihist] = new TH1F("LMQC2LEX","LER MQC2RP x Position", 100, MinLMQC2RPX, MaxLMQC2RPX); // ihsit = 15;
            camPVhists[ihist]->SetTitle("LER MQC2RP x Position");
            camPVhists[ihist]->SetXTitle("LER MQC2RP x (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("LMQC1LEX","LER MQC1RP x Position", 100, MinLMQC1RPX, MaxLMQC1RPX); // ihsit = 16;
            camPVhists[ihist]->SetTitle("LER MQC1RP x Position");
            camPVhists[ihist]->SetXTitle("LER MQC1RP x (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST"); 
            camPVhists[ihist]->SetStats(0);  
            ++ihist;

         camPVhists[ihist] = new TH1F("LMQC1REX","LER MQC1LP x Position", 100, MinLMQC1LPX, MaxLMQC1LPX); // ihsit = 17;
            camPVhists[ihist]->SetTitle("LER MQC1LP x Position");
            camPVhists[ihist]->SetXTitle("LER MQC1LP x (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("LMQC2REX","LER MQC2LP x Position", 100, MinLMQC2LPX, MaxLMQC2LPX); // ihsit = 18;
            camPVhists[ihist]->SetTitle("LER MQC2LP x Position");
            camPVhists[ihist]->SetXTitle("LER MQC2LP x (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("LMQC2LEY","LER MQC2RP y Position", 100, MinLMQC2RPY, MaxLMQC2RPY); // ihsit = 19;
            camPVhists[ihist]->SetTitle("LER MQC2RP y Position");
            camPVhists[ihist]->SetXTitle("LER MQC2RP y (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST"); 
            camPVhists[ihist]->SetStats(0);  
            ++ihist;

         camPVhists[ihist] = new TH1F("LMQC1LEY","LER MQC1RP y Position", 100, MinLMQC1RPY, MaxLMQC1RPY); // ihsit = 20;
            camPVhists[ihist]->SetTitle("LER MQC1RP y Position");
            camPVhists[ihist]->SetXTitle("LER MQC1RP y (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("LMQC1REY","LER MQC1LP y Position", 100, MinLMQC1LPY, MaxLMQC1LPY); // ihsit = 21;
            camPVhists[ihist]->SetTitle("LER MQC1LP y Position");
            camPVhists[ihist]->SetXTitle("LER MQC1LP y (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("LMQC2REY","LER MQC2LP y Position", 100, MinLMQC2LPY, MaxLMQC2LPY); // ihsit = 22;
            camPVhists[ihist]->SetTitle("LER MQC2LP y Position");
            camPVhists[ihist]->SetXTitle("LER MQC2LP y (cm)");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         //Camera 1
         camPVhists[ihist] = new TH1F("Cam1NPeaks","Camera 1 Number of Peaks distribution", 100, camminnp[0], cammaxnp[0]); // ihsit = 23;
            camPVhists[ihist]->SetTitle("Camera 1 #Peaks distribution");
            camPVhists[ihist]->SetXTitle("#Peaks");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("Cam1NPeaksPE","Camera 1 Number of Peaks per Event", EventCount, 1, EventCount); // ihsit = 24;
            camPVhists[ihist]->SetTitle("Camera 1 #Peaks per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("#Peaks");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 
            
         camPVhists[ihist] = new TH1F("Cam1RGBE","Camera 1 Highest RGB per Event", EventCount, 1, EventCount); // ihsit = 25;
            camPVhists[ihist]->SetTitle("Camera 1 Highest RGB per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("RGB Values");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 

         camPVhists[ihist] = new TH1F("Cam1HHdist","Camera 1 Highest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 26;
            camPVhists[ihist]->SetTitle("Camera 1 Highest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
            
         camPVhists[ihist] = new TH1F("Cam1HLdist","Camera 1 Highest and Lowest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 27;
            camPVhists[ihist]->SetTitle("Camera 1 Highest and Lowest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
            
         //Camera 2
         camPVhists[ihist] = new TH1F("Cam2NPeaks","Camera 2 Number of Peaks distribution", 100, camminnp[1], cammaxnp[1]); // ihsit = 28;
            camPVhists[ihist]->SetTitle("Camera 2 #Peaks distribution");
            camPVhists[ihist]->SetXTitle("#Peaks");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("Cam2NPeaksPE","Camera 2 Number of Peaks per Event", EventCount, 1, EventCount); // ihsit = 29;
            camPVhists[ihist]->SetTitle("Camera 2 #Peaks per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("#Peaks");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 
            
         camPVhists[ihist] = new TH1F("Cam2RGBE","Camera 2 Highest RGB per Event", EventCount, 1, EventCount); // ihsit = 30;
            camPVhists[ihist]->SetTitle("Camera 2 Highest RGB per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("RGB Values");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 

         camPVhists[ihist] = new TH1F("Cam2HHdist","Camera 2 Highest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 31;
            camPVhists[ihist]->SetTitle("Camera 2 Highest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
            
         camPVhists[ihist] = new TH1F("Cam2HLdist","Camera 2 Highest and Lowest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 32;
            camPVhists[ihist]->SetTitle("Camera 2 Highest and Lowest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         //Camera 3
         camPVhists[ihist] = new TH1F("Cam3NPeaks","Camera 3 Number of Peaks distribution", 100, camminnp[2], cammaxnp[2]); // ihsit = 33;
            camPVhists[ihist]->SetTitle("Camera 3 #Peaks distribution");
            camPVhists[ihist]->SetXTitle("#Peaks");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("Cam3NPeaksPE","Camera 3 Number of Peaks per Event", EventCount, 1, EventCount); // ihsit = 34;
            camPVhists[ihist]->SetTitle("Camera 3 #Peaks per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("#Peaks");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 
            
         camPVhists[ihist] = new TH1F("Cam3RGBE","Camera 3 Highest RGB per Event", EventCount, 1, EventCount); // ihsit = 35;
            camPVhists[ihist]->SetTitle("Camera 3 Highest RGB per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("RGB Values");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 

         camPVhists[ihist] = new TH1F("Cam3HHdist","Camera 3 Highest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 36;
            camPVhists[ihist]->SetTitle("Camera 3 Highest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
            
         camPVhists[ihist] = new TH1F("Cam3HLdist","Camera 3 Highest and Lowest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 37;
            camPVhists[ihist]->SetTitle("Camera 3 Highest and Lowest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         //Camera 4
         camPVhists[ihist] = new TH1F("Cam4NPeaks","Camera 4 Number of Peaks distribution", 100, camminnp[3], cammaxnp[3]); // ihsit = 38;
            camPVhists[ihist]->SetTitle("Camera 4 #Peaks distribution");
            camPVhists[ihist]->SetXTitle("#Peaks");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("Cam4NPeaksPE","Camera 4 Number of Peaks per Event", EventCount, 1, EventCount); // ihsit = 39;
            camPVhists[ihist]->SetTitle("Camera 4 #Peaks per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("#Peaks");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 
            
         camPVhists[ihist] = new TH1F("Cam4RGBE","Camera 4 Highest RGB per Event", EventCount, 1, EventCount); // ihsit = 40;
            camPVhists[ihist]->SetTitle("Camera 4 Highest RGB per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("RGB Values");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 

         camPVhists[ihist] = new TH1F("Cam4HHdist","Camera 4 Highest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 41;
            camPVhists[ihist]->SetTitle("Camera 4 Highest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
            
         camPVhists[ihist] = new TH1F("Cam4HLdist","Camera 4 Highest and Lowest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 42;
            camPVhists[ihist]->SetTitle("Camera 4 Highest and Lowest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         //Camera 5
         camPVhists[ihist] = new TH1F("Cam5NPeaks","Camera 5 Number of Peaks distribution", 100, camminnp[4], cammaxnp[4]); // ihsit = 43;
            camPVhists[ihist]->SetTitle("Camera 5 #Peaks distribution");
            camPVhists[ihist]->SetXTitle("#Peaks");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("Cam5NPeaksPE","Camera 5 Number of Peaks per Event", EventCount, 1, EventCount); // ihsit = 44;
            camPVhists[ihist]->SetTitle("Camera 5 #Peaks per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("#Peaks");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 
            
         camPVhists[ihist] = new TH1F("Cam5RGBE","Camera 5 Highest RGB per Event", EventCount, 1, EventCount); // ihsit = 45;
            camPVhists[ihist]->SetTitle("Camera 5 Highest RGB per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("RGB Values");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 

         camPVhists[ihist] = new TH1F("Cam5HHdist","Camera 5 Highest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 46;
            camPVhists[ihist]->SetTitle("Camera 5 Highest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
            
         camPVhists[ihist] = new TH1F("Cam5HLdist","Camera 5 Highest and Lowest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 47;
            camPVhists[ihist]->SetTitle("Camera 5 Highest and Lowest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         //Camera 6
         camPVhists[ihist] = new TH1F("Cam6NPeaks","Camera 6 Number of Peaks distribution", 100, camminnp[5], cammaxnp[5]); // ihsit = 48;
            camPVhists[ihist]->SetTitle("Camera 6 #Peaks distribution");
            camPVhists[ihist]->SetXTitle("#Peaks");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("Cam6NPeaksPE","Camera 6 Number of Peaks per Event", EventCount, 1, EventCount); // ihsit = 49;
            camPVhists[ihist]->SetTitle("Camera 6 #Peaks per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("#Peaks");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 
            
         camPVhists[ihist] = new TH1F("Cam6RGBE","Camera 6 Highest RGB per Event", EventCount, 1, EventCount); // ihsit = 50;
            camPVhists[ihist]->SetTitle("Camera 6 Highest RGB per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("RGB Values");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 

         camPVhists[ihist] = new TH1F("Cam6HHdist","Camera 6 Highest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 51;
            camPVhists[ihist]->SetTitle("Camera 6 Highest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
            
         camPVhists[ihist] = new TH1F("Cam6HLdist","Camera 6 Highest and Lowest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 52;
            camPVhists[ihist]->SetTitle("Camera 6 Highest and Lowest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
            
         //Camera 7
         camPVhists[ihist] = new TH1F("Cam7NPeaks","Camera 7 Number of Peaks distribution", 100, camminnp[6], cammaxnp[6]); // ihsit = 53;
            camPVhists[ihist]->SetTitle("Camera 7 #Peaks distribution");
            camPVhists[ihist]->SetXTitle("#Peaks");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("Cam7NPeaksPE","Camera 7 Number of Peaks per Event", EventCount, 1, EventCount); // ihsit = 54;
            camPVhists[ihist]->SetTitle("Camera 7 #Peaks per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("#Peaks");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 
            
         camPVhists[ihist] = new TH1F("Cam7RGBE","Camera 7 Highest RGB per Event", EventCount, 1, EventCount); // ihsit = 55;
            camPVhists[ihist]->SetTitle("Camera 7 Highest RGB per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("RGB Values");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 

         camPVhists[ihist] = new TH1F("Cam7HHdist","Camera 7 Highest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 56;
            camPVhists[ihist]->SetTitle("Camera 7 Highest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
            
         camPVhists[ihist] = new TH1F("Cam7HLdist","Camera 7 Highest and Lowest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 57;
            camPVhists[ihist]->SetTitle("Camera 7 Highest and Lowest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
            
         //Camera 8
         camPVhists[ihist] = new TH1F("Cam8NPeaks","Camera 8 Number of Peaks distribution", 100, camminnp[7], cammaxnp[7]); // ihsit = 58;
            camPVhists[ihist]->SetTitle("Camera 8 #Peaks distribution");
            camPVhists[ihist]->SetXTitle("#Peaks");
            camPVhists[ihist]->SetYTitle("dN");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;

         camPVhists[ihist] = new TH1F("Cam8NPeaksPE","Camera 8 Number of Peaks per Event", EventCount, 1, EventCount); // ihsit = 59;
            camPVhists[ihist]->SetTitle("Camera 8 #Peaks per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("#Peaks");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 
            
         camPVhists[ihist] = new TH1F("Cam8RGBE","Camera 8 Highest RGB per Event", EventCount, 1, EventCount); // ihsit = 60;
            camPVhists[ihist]->SetTitle("Camera 8 Highest RGB per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("RGB Values");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist; 

         camPVhists[ihist] = new TH1F("Cam8HHdist","Camera 8 Highest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 61;
            camPVhists[ihist]->SetTitle("Camera 8 Highest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
            
         camPVhists[ihist] = new TH1F("Cam8HLdist","Camera 8 Highest and Lowest RGB peaks distance per event", EventCount, 1, EventCount); // ihsit = 62;
            camPVhists[ihist]->SetTitle("Camera 8 Highest and Lowest RGB peaks distance per Event");
            camPVhists[ihist]->SetXTitle("Event");
            camPVhists[ihist]->SetYTitle("distance");
            camPVhists[ihist]->SetOption("HIST");
            camPVhists[ihist]->SetStats(0);
            ++ihist;
      }

      for (int j=0; j<22; j++){
         camPVhists[j+1]->Fill(Allpvs[i][j]);      
      }
      
      histcount=23;
      for (int k=0; k<8; k++){
         if(camminnp[k]>EventData[i][k][0][0]) camminnp[k]=EventData[i][k][0][0];
         if(cammaxnp[k]<EventData[i][k][0][0]) cammaxnp[k]=EventData[i][k][0][0];
         
         if (EventData[i][k][0][0]!=0) camPVhists[histcount]->Fill(EventData[i][k][0][0]);
         histcount++;

         camPVhists[histcount]->Fill(i, EventData[i][k][0][0]);
         histcount++;

         for (int l=1; l<=EventData[i][k][0][0]; l++){
            if(EventData[i][k][0][0]>=2){
               if(maxrgb<EventData[i][k][l][0]){
                  max2rgb=maxrgb;
                  max2X=maxX;
                  max2Y=maxY;
               }else if(max2rgb<EventData[i][k][l][0]){
                  max2rgb=EventData[i][k][l][0];
                  max2X=EventData[i][k][l][1];
                  max2Y=EventData[i][k][l][2];     
               }
               histcheck1=0;
            }else{
               histcheck1=1;
            }
         
            if(maxrgb<EventData[i][k][l][0]){            
               maxrgb=EventData[i][k][l][0];
               maxX=EventData[i][k][l][1];
               maxY=EventData[i][k][l][2];
            }
            if(EventData[i][k][0][0]>=3){
               if(minrgb>EventData[i][k][l][0]){ 
                  minrgb=EventData[i][k][l][0];
                  minX=EventData[i][k][l][1];
                  minY=EventData[i][k][l][2];
               }
               histcheck2=0;
            }else{
               histcheck2=1;
            }
         }   
               
         camPVhists[histcount]->Fill(i, maxrgb);
         maxrgb=0; minrgb=0;
         histcount++;
               
         if (histcheck1==0){
            calc=((maxX - max2X)*(maxX - max2X))+((maxY - max2Y)*(maxY - max2Y));
            calc=TMath::Sqrt(calc);
            camPVhists[histcount]->Fill(i, calc);
         }
         calc=0; max2X=0; max2Y=0;
         histcount++;
               
         if (histcheck2==0){
            calc=((maxX - minX)*(maxX - minX))+((maxY - minY)*(maxY - minY));
            calc=TMath::Sqrt(calc);
            camPVhists[histcount]->Fill(i, calc);
         }
         calc=0; maxX=0; maxY=0; minX=0; minY=0;
         histcount++;
      } 

      if (i==eventsnum) { //Only draw histograms once, when its the last event
         //Print the histograms
	
         camPVhists[1]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("TimeStamp.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[2]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("HCurrent.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[3]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("LCurrent.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[4]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("HLumi.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[5]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("LLumi.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[6]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("ECLLumi.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[7]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("HMQC2LEX.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[8]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("HMQC1LEX.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[9]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("HMQC1REX.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[10]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("HMQC2REX.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[11]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("HMQC2LEY.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[12]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("HMQC1LEY.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[13]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("HMQC1REY.pdf");      
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[14]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("HMQC2REY.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[15]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("LMQC2LEX.pdf");      
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[16]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("LMQC1LEX.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[17]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("LMQC1REX.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[18]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("LMQC2REX.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[19]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("LMQC2LEY.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[20]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("LMQC1LEY.pdf");
         c1->SaveAs(pdfname.c_str());      
         c1->Clear();

         camPVhists[21]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("LMQC1REY.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[22]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("LMQC2REY.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();
         
         camPVhists[23]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam1NPeakDis.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[24]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam1NPeakEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[25]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam1MaxRGBEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[26]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam1MaxPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[27]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam1MaxMinPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[28]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam2NPeakDis.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[29]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam2NPeakEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[30]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam2MaxRGBEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[31]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam2MaxPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[32]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam2MaxMinPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[33]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam3NPeakDis.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[34]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam3NPeakEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[35]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam3MaxRGBEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[36]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam3MaxPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[37]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam3MaxMinPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[38]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam4NPeakDis.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[39]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam4NPeakEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[40]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam4MaxRGBEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[41]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam4MaxPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[42]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam4MaxMinPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[43]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam5NPeakDis.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[44]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam5NPeakEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[45]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam5MaxRGBEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[46]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam5MaxPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[47]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam5MaxMinPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[48]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam6NPeakDis.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[49]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam6NPeakEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[50]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam6MaxRGBEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[51]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam6MaxPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[52]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam6MaxMinPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[53]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam7NPeakDis.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[54]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam7NPeakEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[55]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam7MaxRGBEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[56]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam7MaxPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[57]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam7MaxMinPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[58]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam8NPeakDis.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[59]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam8NPeakEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[60]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam8MaxRGBEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[61]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam8MaxPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         camPVhists[62]->Draw("HIST");
         c1->Update();
         pdfname = histdir + string("Cam8MaxMinPeaksrEvnt.pdf");
         c1->SaveAs(pdfname.c_str());
         c1->Clear();

         for (int m=1; m<=62; m++){
            delete camPVhists[m];
         }
      }
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void createNtuple(int evntsnm) {
   
   TFile *labmfile = new TFile(rootfilename, "RECREATE");
   TTree *labmtree = new TTree(treename, "LABM Cameras branches");

   int TimeStampT;
   float HERCurrentT;
   float LERCurrentT; 
   float HERLumiT;
   float LERLumiT;
   float ECLLumiT;
   float HERMQC2LEX;
   float HERMQC1LEX;
   float HERMQC1REX;
   float HERMQC2REX;
   float HERMQC2LEY;
   float HERMQC1LEY;
   float HERMQC1REY;
   float HERMQC2REY;
   float LERMQC2RPX;
   float LERMQC1RPX;
   float LERMQC1LPX;
   float LERMQC2LPX;
   float LERMQC2RPY;
   float LERMQC1RPY;
   float LERMQC1LPY;
   float LERMQC2LPY;  
   
   string camera1png, camera2png, camera3png, camera4png, camera5png, camera6png, camera7png, camera8png;
   
   int npeakcam1, npeakcam2, npeakcam3, npeakcam4, npeakcam5, npeakcam6, npeakcam7, npeakcam8;
   
   labmtree->Branch("TimeStampT", &TimeStampT, "TimeStampT/I");
   labmtree->Branch("HERCurrentT", &HERCurrentT, "HERCurrentT/F");
   labmtree->Branch("LERCurrentT", &LERCurrentT, "LERCurrentT/F");
   labmtree->Branch("HERLumiT", &HERLumiT, "HERLumiT/F");
   labmtree->Branch("LERLumiT", &LERLumiT, "LERLumiT/F");
   labmtree->Branch("ECLLumiT", &ECLLumiT, "ECLLumiT/F");
   labmtree->Branch("HERMQC2LEX", &HERMQC2LEX, "HERMQC2LEX/F");
   labmtree->Branch("HERMQC1LEX", &HERMQC1LEX, "HERMQC1LEX/F");
   labmtree->Branch("HERMQC1REX", &HERMQC1REX, "HERMQC1REX/F");
   labmtree->Branch("HERMQC2REX", &HERMQC2REX, "HERMQC2REX/F");
   labmtree->Branch("HERMQC2LEY", &HERMQC2LEY, "HERMQC2LEY/F");
   labmtree->Branch("HERMQC1LEY", &HERMQC1LEY, "HERMQC1LEY/F");
   labmtree->Branch("HERMQC1REY", &HERMQC1REY, "HERMQC1REY/F");
   labmtree->Branch("HERMQC2REY", &HERMQC2REY, "HERMQC2REY/F");
   labmtree->Branch("LERMQC2RPX", &LERMQC2RPX, "LERMQC2RPX/F");
   labmtree->Branch("LERMQC1RPX", &LERMQC1RPX, "LERMQC1RPX/F");
   labmtree->Branch("LERMQC1LPX", &LERMQC1LPX, "LERMQC1LPX/F");
   labmtree->Branch("LERMQC2LPX", &LERMQC2LPX, "LERMQC2LPX/F");
   labmtree->Branch("LERMQC2RPY", &LERMQC2RPY, "LERMQC2RPY/F");
   labmtree->Branch("LERMQC1RPY", &LERMQC1RPY, "LERMQC1RPY/F");
   labmtree->Branch("LERMQC1LPY", &LERMQC1LPY, "LERMQC1LPY/F");
   labmtree->Branch("LERMQC2LPY", &LERMQC2LPY, "LERMQC2LPY/F");
   labmtree->Branch("camera1png", &camera1png);
   labmtree->Branch("camera2png", &camera2png);
   labmtree->Branch("camera3png", &camera3png);
   labmtree->Branch("camera4png", &camera4png);
   labmtree->Branch("camera5png", &camera5png);
   labmtree->Branch("camera6png", &camera6png);
   labmtree->Branch("camera7png", &camera7png);
   labmtree->Branch("camera8png", &camera8png);
   labmtree->Branch("npeakcam1", &npeakcam1, "npeakcam1/I");
   labmtree->Branch("npeakcam2", &npeakcam2, "npeakcam2/I");
   labmtree->Branch("npeakcam3", &npeakcam3, "npeakcam3/I");
   labmtree->Branch("npeakcam4", &npeakcam4, "npeakcam4/I");
   labmtree->Branch("npeakcam5", &npeakcam5, "npeakcam5/I");
   labmtree->Branch("npeakcam6", &npeakcam6, "npeakcam6/I");
   labmtree->Branch("npeakcam7", &npeakcam7, "npeakcam7/I");
   labmtree->Branch("npeakcam8", &npeakcam8, "npeakcam8/I");
   
   for(int i=1; i<=evntsnm; i++){

      TimeStampT=Allpvs[i][0];
      HERCurrentT=Allpvs[i][1];
      LERCurrentT=Allpvs[i][2]; 
      HERLumiT=Allpvs[i][3];
      LERLumiT=Allpvs[i][4];
      ECLLumiT=Allpvs[i][5];
      HERMQC2LEX=Allpvs[i][6];
      HERMQC1LEX=Allpvs[i][7];
      HERMQC1REX=Allpvs[i][8];
      HERMQC2REX=Allpvs[i][9];
      HERMQC2LEY=Allpvs[i][10];
      HERMQC1LEY=Allpvs[i][11];
      HERMQC1REY=Allpvs[i][12];
      HERMQC2REY=Allpvs[i][13];
      LERMQC2RPX=Allpvs[i][14];
      LERMQC1RPX=Allpvs[i][15];
      LERMQC1LPX=Allpvs[i][16];
      LERMQC2LPX=Allpvs[i][17];
      LERMQC2RPY=Allpvs[i][18];
      LERMQC1RPY=Allpvs[i][19];
      LERMQC1LPY=Allpvs[i][20];
      LERMQC2LPY=Allpvs[i][21];  
      
      camera1png=pngnames[i][0];
      camera2png=pngnames[i][1];
      camera3png=pngnames[i][2];
      camera4png=pngnames[i][3];
      camera5png=pngnames[i][4];
      camera6png=pngnames[i][5];
      camera7png=pngnames[i][6];
      camera8png=pngnames[i][7];
   
      npeakcam1=int(EventData[i][0][0][0]);
      npeakcam2=int(EventData[i][1][0][0]);
      npeakcam3=int(EventData[i][2][0][0]);
      npeakcam4=int(EventData[i][3][0][0]);
      npeakcam5=int(EventData[i][4][0][0]);
      npeakcam6=int(EventData[i][5][0][0]);
      npeakcam7=int(EventData[i][6][0][0]);
      npeakcam8=int(EventData[i][7][0][0]);
      
      labmtree->Fill();
      
   }

   labmtree->Write();
   labmfile->Close();

}
 
