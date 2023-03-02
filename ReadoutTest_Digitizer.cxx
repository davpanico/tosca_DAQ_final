#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <thread>
#include <memory>
#include <chrono>
#include <ctime>
#include <sstream>
#include <algorithm>
#include "cmath"



#include "CAENDigitizer.h"
#include "keyb.h"


//ROOT Headers
#include "TPaletteAxis.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TObject.h"
#include "TList.h"
#include "TNamed.h"
#include "THttpCallArg.h"
#include "THttpEngine.h"
#include "THttpWSHandler.h"
#include <mutex>
#include <TMemFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TFrame.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <THttpServer.h>
#include "TTree.h"
#include "TString.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <TH1D.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TList.h>
#include <TStyle.h>
#include <TH2.h>
#include <TLegend.h>



#define CAEN_USE_DIGITIZERS
#define IGNORE_DPP_DEPRECATED
using std::string;
using std::stringstream;

//#define WFRMPLOT
#define MATRICES
#define WRITEFILE



#define MAXNB 1         /* Number of connected boards */

int checkCommand() {
	int c = 0;
	if(!kbhit())
			return 0;

	c = getch();
	switch (c) {
		case 's': 
			return 9;
			break;
		case 'k':
			return 1;
			break;
		case 'q':
			return 2;
			break;
	}
	return 0;
}









//++++++++++++++++++++++++++++++++++++++++++++++++++++ PULSE PROCESSING +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define beg_baseline_fit_length  20
#define beg_baseline_fit_length_baf2  40
#define end_baseline_fit_length 100
#define the_noise_level 30
float threshold	= 2600;

struct baseline_data{
  Bool_t good_beg,good_end,stable,ok;
  float value;
};


float abss(float xx){if(xx<0){return -xx;}return xx;}


			      

//NEGSIG: GET THRESHOLD X IN SAMPLES
Short_t get_thr_x(float wf[]){    
  for(Int_t ii=0;ii<400;ii++){                             //Max sample to look for the threshold
    if(wf[ii]<threshold&&wf[ii+1]<threshold){return ii;}
  }
  return -10; 
}


////NEGSIG: LED IN SAMPLES 
Float_t get_led(float wf[]){
    int x_f = 0;
    int	x_s = 0;
    float b,m;
    for(UShort_t indx=1;indx<1023;indx++){
	if(wf[indx]>threshold && wf[indx+1]<threshold){x_f = indx; x_s = indx+1;
	break;}
	}
    m = (wf[x_f]-wf[x_s])/(x_f-x_s);
    b = wf[x_f]-(m*x_f);
    return (threshold-b)/m; 
}



////NEGSIG: LED_OLS IN SAMPLES 
Float_t get_led2(float wf[]){
    int x0 = 0;	
    int x1 = 0;
    int x2 = 0;
    int x3 = 0;
    float b,m;
    for(UShort_t tind=1;tind<1023;tind++){
	if(wf[tind]>threshold && wf[tind+1]<threshold && wf[tind+2]<threshold && wf[tind-1]>threshold){x0 = tind; x1 = tind+1; x2 = tind-1; x3 = tind+2;
	break;}
	}
	float xavg = 0.25 * (x0+x1+x2+x3);
	float yavg = 0.25 * (wf[x0]+wf[x1]+wf[x2]+wf[x3]);
	float sumup = (x0-xavg)*(wf[x0]-yavg)+(x1-xavg)*(wf[x1]-yavg)+(x2-xavg)*(wf[x2]-yavg)+(x3-xavg)*(wf[x3]-yavg);
	float sumdwn = pow(2,x0-xavg)+pow(2,x1-xavg)+ pow(2,x2-xavg) + pow(2,x3-xavg);
        m  = sumup/sumdwn ; 
        b = yavg - m*xavg ;
    return (threshold-b)/m; 
}




////NEGSIG: CFD IN SAMPLES 
Float_t get_cft(float wf2[],UShort_t risetime,Short_t thr_x,Float_t fraction, float baseline){
  UShort_t points=600;
  float wf[1023];
  Float_t delayed[points],inv[points],sum[points];
  UShort_t delay=(UShort_t)((1.-fraction)*(1.*risetime)+.5);
  for(UShort_t ip=0;ip<1023;ip++){
    wf[ip] = -wf2[ip]+baseline;
  }
  for(UShort_t ii=0;ii<points;ii++){
    inv[ii]=-fraction*wf[thr_x-1+ii];
    delayed[ii]=wf[thr_x-1-delay+ii];
    sum[ii]=inv[ii]+delayed[ii];
  }
  UShort_t jb0=0,ja0=0;
  for(UShort_t ii=1;ii<points;ii++){if(sum[ii]>0){ja0=ii;jb0=ii-1;break;}}
//  printf("sum A: %f 	sum B: %f \n",sum[ja0],sum[jb0]);
  if((int)(sum[ja0]-sum[jb0]*100)==0){return -10;}
  Float_t zc=-sum[jb0]/(sum[ja0]-sum[jb0])+jb0;
  return (1.*thr_x-1+zc);
}

////NEGSIG: CFD IN SAMPLES TEST DAVIDE
//Float_t get_cft(float wf2[],UShort_t risetime,Short_t thr_x,Float_t fraction, float baseline){
//  UShort_t points=1000;
//  float wf[1023];
//  Float_t delayed[points],inv[points],sum[points];
//  UShort_t delay=(UShort_t)(0.5*(1.*risetime));
//  for(UShort_t ip=0;ip<1023;ip++){
//    wf[ip] = -wf2[ip]+baseline;
//  }
//  for(UShort_t ii=0;ii<points;ii++){
//    inv[ii]=-fraction*wf[ii];
//    delayed[ii]=wf[ii-delay];
//    sum[ii]=inv[ii]+delayed[ii];
//  }
//  UShort_t jb0=0,ja0=0;
//  for(UShort_t ii=1;ii<points;ii++){if(sum[ii]>0){ja0=ii;jb0=ii-1;break;}}
//  printf("sum A: %f 	sum B: %f \n",sum[ja0],sum[jb0]);
//  if((int)(sum[ja0]-sum[jb0]*100)==0){return -10;}
//  Float_t zc=sum[jb0]/(sum[ja0]-sum[jb0])-jb0;
//  return (zc);
//}



//NEGSIG: BASELINE CALC
baseline_data get_baseline(float wf[]){
  UShort_t samples = 1023;
  TGraph *gr = new TGraph(1023); 
  UShort_t min_beg=4096,min_end=4096,max_beg=0,max_end=0;
  for(UShort_t ii=0;ii<samples;ii++){
    gr->SetPoint(ii,ii,wf[ii]);
    if(ii<beg_baseline_fit_length){
      if(min_beg>wf[ii]){min_beg=wf[ii];}
      if(max_beg<wf[ii]){max_beg=wf[ii];}}
    else if(ii>samples-1-end_baseline_fit_length){
      if(min_end>wf[ii]){min_end=wf[ii];}
      if(max_end<wf[ii]){max_end=wf[ii];}}
  }
  baseline_data b;
  b.value=0.;
  b.stable=false;
  if(max_beg-min_beg>the_noise_level){b.good_beg=false;}else{b.good_beg=true;}
  if(max_end-min_end>the_noise_level){b.good_end=false;}else{b.good_end=true;}
  Float_t base_line_beg=0.,base_line_end=0.;
  if(b.good_beg){TF1 *fbase0 = new TF1("fbase0","pol0",0,beg_baseline_fit_length);gr->Fit("fbase0","QRNC");base_line_beg=fbase0->GetParameter(0);delete fbase0;}
  if(b.good_end){TF1 *fbase1 = new TF1("fbase1","pol0",samples-1-end_baseline_fit_length,samples-1);gr->Fit("fbase1","QRNC");
  base_line_end=fbase1->GetParameter(0);delete fbase1;}
  delete gr;
  if(b.good_beg&&b.good_end){if(abss(base_line_beg-base_line_end)>the_noise_level){b.stable=false;}else{b.stable=true;}}
  if(b.stable||b.good_beg){b.value=base_line_beg;}
  if(!b.good_beg&&b.good_end){b.value=base_line_end;}
  if(b.good_beg||b.good_end){b.ok=true;}else{b.ok=false;}
  return b;
}





   



int main(int argc, char* argv[])
{	CAEN_DGTZ_ErrorCode ret;
	int	handle[MAXNB];
	CAEN_DGTZ_BoardInfo_t BoardInfo;
	CAEN_DGTZ_EventInfo_t eventInfo;
	//CAEN_DGTZ_UINT16_EVENT_t *Evt = NULL;
	CAEN_DGTZ_X742_EVENT_t *Evt = NULL;
	char *buffer = NULL;
	int MajorNumber;
	int i,b;
	int c = 0;
	int Nevent = 0;
	char * evtptr = NULL;
	uint32_t size,bsize;
	uint32_t numEvents;
	Short_t threshold_sampl[32];
	Float_t CFD_sampl[32];
	float bsln[32];
	i = sizeof(CAEN_DGTZ_TriggerMode_t);
	
 	FILE *fptr;
	fptr = fopen("/run/media/root/JYU23/TEST/test_plot.txt","w");   // TO EDIT BASED ON THE MACHINE !!!!!
	//fptr = fopen("/run/media/root/JYU23/DATA/test_plot.txt","w");   // TO EDIT BASED ON THE MACHINE !!!!!
	
	
	
	/*******ROOT INIT******/
	
	TApplication app("App", &argc, argv);
	float timearr[1024];

	for (int ipt = 0; ipt<1023; ipt ++) {
	timearr[ipt] = ipt*.2  ;        // sample in ns
	}



        #ifdef WFRMPLOT

	TCanvas *canvas = new TCanvas("Live Waveforms", "Live Waveforms",200,10,800,600); 
	canvas->SetGrid();
	TH1F *f = canvas->DrawFrame(0,0,220,5000,"Waveform Monitor; Time [ns]; ADC Channel");
	#endif


	#ifdef MATRICES
	//TApplication appm("App", &argc, argv);
	
	Double_t w = 800;
    	Double_t h = 600;
	TCanvas *canvas2 = new TCanvas("Matrices", "Matrices",w,h);
    	canvas2->SetWindowSize(w + 4, h + 28);
 	auto toftof = new TH2F("toftof","",200,28,34,200,28,34);
	toftof->GetYaxis()->SetTitle("ARM1 TOF [ns]");
	toftof->GetXaxis()->SetTitle("ARM2 TOF [ns]");
	toftof->SetStats(0);
	gStyle->SetPalette(1);
	gPad->Update();
	canvas2->Update();
	toftof->Draw("COLZ");
	canvas2->Update();
	#endif

        //THttpServer* server = new THttpServer("http:8080");
	
	stringstream st;
    	string nomefile;
	time_t now = time(NULL);


    	st <<"file:/run/media/root/JYU23/TEST/"<< "TOSCA_"<< now <<".root";
        //st <<"file:/run/media/root/JYU23/DATA/"<< "TOSCA_"<< now <<".root";
            
    	st >> nomefile;
    	const char *nome = nomefile.c_str();

	TFile *rf = new TFile(nome, "recreate");
   	TTree *T = new TTree("T", "Treedist");
// 	float DataTr0[1024];
//	float Ch0[1024];
//	float Ch1[1024];
//	T->Branch("Tr0",&DataTr0,"Tr0[1024]/F");
//	T->Branch("Ch0",&Ch0,"Ch0[1024]/F");
//	T->Branch("Ch1",&Ch1,"Ch1[1024]/F");
//	T->Branch("time",&timearr,"time[1024]/F");
	T->Branch("TCh1",&CFD_sampl[1],"TCh1/F");
	T->Branch("TCh0",&CFD_sampl[0],"TCh0/F");
	T->Branch("TCh16",&CFD_sampl[16],"TCh16/F");
	T->Branch("TCh17",&CFD_sampl[17],"TCh17/F");
	float TOF_CFD1;   //ARM1 sp1-st1  ch1-ch0
	float TOF_CFD2;     //ARM2 sp2-st2  ch9-ch8
	T->Branch("TOF_CFD1",&TOF_CFD1,"TOF_CFD1/F");
	T->Branch("TOF_CFD2",&TOF_CFD2,"TOF_CFD2/F");
	
						     

   



/* **************************************************************** Digitizer INIT Routine ********************************************************** */

    for(b=0; b<MAXNB; b++){
        

        /* ************** Print SW Info ********************* */
        

	printf("\n ********************************************************************* \n");
          printf(" ************    Readout Software For TOSCA Detectors    ************* \n");
          printf(" ************               ----------                   ************* \n");
          printf(" ************                       for CAEN V1742 board ************* \n");
	  printf(" ************                              Version: 2.23 ************* \n");
          printf(" ************     D. Panico, A. Di Nitto, S. Di Costanzo ************* \n");
	  printf(" ********************************************************************* \n");
                                

	

	/* *************** A4818 INIT *********************** */

	int32_t PID = 23365;
        
	
        ret = CAEN_DGTZ_OpenDigitizer2(CAEN_DGTZ_USB_A4818,&PID,0,0,&handle[b]);     //V1742
	//ret = CAEN_DGTZ_OpenDigitizer2(CAEN_DGTZ_USB_A4818,&PID,0,0,&handle[1]);     //V1725
	
        if(ret != CAEN_DGTZ_Success) {
            printf("\n!!!!!!!!!!!!!         Error: CANNOT OPEN DIGITIZER!     !!!!!!!!!!!!!!!!! \n");
            goto QuitProgram;
        }


        

        /* Once we have the handler to the digitizer, we use it to call the other functions */

        ret = CAEN_DGTZ_GetInfo(handle[b], &BoardInfo);
        printf("\n|-----> Connected to CAEN Digitizer Model %s, recognized as board %d\n", BoardInfo.ModelName, b);
        printf("\t|-----> ROC FPGA Release is %s\n", BoardInfo.ROC_FirmwareRel);
        printf("\t|-----> AMC FPGA Release is %s\n", BoardInfo.AMC_FirmwareRel);
	    
		// Check firmware revision (Quit if DPP is ON) */
		sscanf(BoardInfo.AMC_FirmwareRel, "%d", &MajorNumber);
		if (MajorNumber >= 128) {
			printf("This digitizer has a DPP firmware!\n");
			goto QuitProgram;
		}


			

        /* *************** SAMPLING PARAMETERS ******************** */

        ret = CAEN_DGTZ_Reset(handle[b]);                                               /* Reset Digitizer */
	ret = CAEN_DGTZ_SetDRS4SamplingFrequency(handle[b],CAEN_DGTZ_DRS4_5GHz);        /* Set the sampling frequency (Default = 5GHz) */
	ret = CAEN_DGTZ_GetInfo(handle[b], &BoardInfo);                                 /* Get Board Info */
	ret = CAEN_DGTZ_SetRecordLength(handle[b],1024);                                /* Set the lenght of each waveform (in samples) */
	ret = CAEN_DGTZ_SetGroupEnableMask(handle[b],15);                                /* Enable all groups  */
	//ret = CAEN_DGTZ_WriteRegister(handle[b],
	ret = CAEN_DGTZ_SetPostTriggerSize(handle[b],15);

	ret = CAEN_DGTZ_WriteRegister(handle[b],0x1098,0xF6C00);                        /* DC offset GR0 -> 0x6C00 */
	ret = CAEN_DGTZ_WriteRegister(handle[b],0x1198,0xF6C00);       			/* DC offset GR1 -> 0x6C00 */
	ret = CAEN_DGTZ_WriteRegister(handle[b],0x1298,0xF6C00);			/* DC offset GR2 -> 0x6C00 */
	ret = CAEN_DGTZ_WriteRegister(handle[b],0x1398,0xF6C00);			/* DC offset GR3 -> 0x6C00 */
	

        /* ************** TRIGGER PARAMETERS ********************** */

      //ret = CAEN_DGTZ_SetGroupTriggerThreshold(handle[b],0,32768);                    /* Set selfTrigger threshold */
      //ret = CAEN_DGTZ_SetGroupSelfTrigger(handle[b],CAEN_DGTZ_TRGMODE_ACQ_ONLY,15);   /* Set trigger on channel 0 to be ACQ_ONLY */
      	ret = CAEN_DGTZ_SetSWTriggerMode(handle[b],CAEN_DGTZ_TRGMODE_DISABLED);         /* Set the behaviour when a SW tirgger arrives */

	ret = CAEN_DGTZ_SetFastTriggerMode(handle[b],CAEN_DGTZ_TRGMODE_ACQ_ONLY);       /* Enable the TRn as the local trigger  */
	ret = CAEN_DGTZ_SetFastTriggerDigitizing(handle[b], CAEN_DGTZ_ENABLE);          /* Enable the TRn input signal digitization */
	ret = CAEN_DGTZ_SetGroupFastTriggerDCOffset(handle[b],15, 32768);          	/* Set the TRn input signal DC Offset */
	ret = CAEN_DGTZ_SetGroupFastTriggerThreshold(handle[b],15, 20934);          	/* Set the TRn input signal threshold */

	ret = CAEN_DGTZ_LoadDRS4CorrectionData(handle[b], CAEN_DGTZ_DRS4_5GHz);
	ret = CAEN_DGTZ_EnableDRS4Correction(handle[b]);

	ret = CAEN_DGTZ_SetMaxNumEventsBLT(handle[b],1000000);                              /* Set the max number of events to transfer in a sigle readout */
        ret = CAEN_DGTZ_SetAcquisitionMode(handle[b],CAEN_DGTZ_SW_CONTROLLED);          /* Set the acquisition mode */






	
        if(ret != CAEN_DGTZ_Success) {
            printf("\n !!!!!!!!!!!!  Error during Digitizer Configuration: %d     !!!!!!!!!!!!! \n", ret);
	    Sleep(100);
            goto QuitProgram;
	    
        } else { 
	  
	     printf("\n|----->  Successful Connection \n");
	     Sleep(100);
	}

    }




/* ******************************************* Acquisition Control ***************************************************/


    printf("\n\n|-----> Press 'S' to start the acquisition |\n");
        printf("|-----> Press 'K' to stop  the acquisition |\n");
        printf("|-----> Press 'Q' to quit  the application |\n\n");

    while (1) {
		c = checkCommand();
		if (c == 9) break;
		if (c == 2) return -1;
		Sleep(100);
    }




    
    /*************************************************************** Memory Readout Buffer Allocation *****************************************************/
    printf("|-----> Reading |\n");

    ret = CAEN_DGTZ_MallocReadoutBuffer(handle[b],&buffer,&size);

 
     /*************************************************************** Start Acquisition ********************************************************************/
    
    for(b=0; b<MAXNB; b++){
            ret = CAEN_DGTZ_SWStartAcquisition(handle[b]);
    
    
    /************************************************************** Start acquisition loop ****************************************************************/
    
    
	while(1) {
        for(b=0; b<MAXNB; b++) {
		    

		    ret = CAEN_DGTZ_SendSWtrigger(handle[b]); /* ???? Send a SW Trigger ???? */
		    ret = CAEN_DGTZ_ReadData(handle[b],CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,buffer,&bsize); /* Read the buffer from the digitizer */

		    ret = CAEN_DGTZ_GetNumEvents(handle[b],buffer,bsize,&numEvents);

		   
		    
	


             
		    for (i=0;i<numEvents;i++) {
               				 /* Get the Infos and pointer to the event */
			    ret = CAEN_DGTZ_GetEventInfo(handle[b],buffer,bsize,i,&eventInfo,&evtptr);

               				 /* Decode the event to get the data */

			    //*************************************
			    // Event Elaboration
			    //*************************************

			      ret = CAEN_DGTZ_AllocateEvent(handle[b],(void**) &Evt);
			      ret = CAEN_DGTZ_DecodeEvent(handle[b],evtptr,(void**) &Evt);
			      Nevent +=1 ;

			
				threshold_sampl[0] = get_thr_x((Evt->DataGroup[0]).DataChannel[0]);
				bsln[0] = get_baseline((Evt->DataGroup[0]).DataChannel[0]).value;
				CFD_sampl[0] = get_cft((Evt->DataGroup[0]).DataChannel[0], 20, threshold_sampl[0], 0.45, bsln[0]);

				threshold_sampl[1] = get_thr_x((Evt->DataGroup[0]).DataChannel[1]);
				bsln[1] = get_baseline((Evt->DataGroup[0]).DataChannel[1]).value;
				CFD_sampl[1] = get_cft((Evt->DataGroup[0]).DataChannel[1], 20, threshold_sampl[1], 0.45, bsln[1]);

				threshold_sampl[17] = get_thr_x((Evt->DataGroup[2]).DataChannel[1]);
				bsln[17] = get_baseline((Evt->DataGroup[2]).DataChannel[1]).value;
				CFD_sampl[17] = get_cft((Evt->DataGroup[2]).DataChannel[1], 20, threshold_sampl[17], 0.45, bsln[17]);

				threshold_sampl[16] = get_thr_x((Evt->DataGroup[2]).DataChannel[0]);
				bsln[16] = get_baseline((Evt->DataGroup[2]).DataChannel[0]).value;
				CFD_sampl[16] = get_cft((Evt->DataGroup[2]).DataChannel[0], 20, threshold_sampl[16], 0.45, bsln[16]);


		
					
			
				TOF_CFD1 = (CFD_sampl[1]-CFD_sampl[0])*.2;    // Arm1 (ch1-ch0)
				TOF_CFD2 = (CFD_sampl[17]-CFD_sampl[16])*.2;    // Arm2 (ch17-ch16)
				
				

				printf("---> EVENT #%d | TOF ARM1: %f | TOF ARM2: %f \n",Nevent,TOF_CFD1,TOF_CFD2);




			      		/* LIVE PLOT */

					#ifdef WFRMPLOT

					
					
					auto g1 = new TGraph(1023,timearr,((Evt->DataGroup[0]).DataChannel[0]));
					g1->SetLineColor(4);
					g1->SetLineWidth(2);
					
					auto g2 = new TGraph(1023,timearr,((Evt->DataGroup[0]).DataChannel[1]));
					g2->SetLineWidth(2);
					g2->SetLineColor(0);

					auto g3 = new TGraph(1023,timearr,((Evt->DataGroup[1]).DataChannel[0]));
					g3->SetLineWidth(2);
					g3->SetLineColor(1);

					auto g4 = new TGraph(1023,timearr,((Evt->DataGroup[1]).DataChannel[1]));
					g4->SetLineWidth(2);
					g4->SetLineColor(3);
					
                                      

					g1->Draw("L");
					g2->Draw("same");
					g3->Draw("same");
					g4->Draw("same");
					auto legend = new TLegend();
					legend->AddEntry(g4,"ARM2: STOP");
					legend->AddEntry(g2,"ARM1: STOP");
					legend->AddEntry(g3,"ARM2: START");
					legend->AddEntry(g1,"ARM1: START");
					legend->Draw();
					//serv->Register("graph/subfolder",g1);	
					canvas->Update();
					//getchar();
					//usleep(1000000);
					delete g1;
					delete g2;
					delete g3;
					delete g4;
					delete legend;
					
					#endif






					/* Plot Matrices */

					#ifdef MATRICES

					toftof->Fill(TOF_CFD2,TOF_CFD1);
					toftof->Draw("COLZ");
					canvas2->Update();

					
					
					

					#endif






			    		
					/* Write data to file */
					#ifdef WRITEFILE
					
//					for (int dp = 0; dp < 1023; dp++) {     
//						      DataTr0[dp] = ((Evt->DataGroup[0]).DataChannel[8])[dp];
//							Ch0[dp] = ((Evt->DataGroup[0]).DataChannel[0])[dp];
//							Ch1[dp] = ((Evt->DataGroup[0]).DataChannel[1])[dp];     
//							}
//					
					T->Fill();
					

		 			      //for (int k = 0; k < 2; k++){
						//for (int j = 0; j < 9; j++) { 
						  //fprintf(fptr,"\nGroup %d | Ch %d : ",k,j);
						    //for (int i = 0; i < 1023; i++) {     
						      //fprintf(fptr," %f ", ((Evt->DataGroup[k]).DataChannel[j])[i]);     
							//}      
						      //}
						    //}


					#endif
			

			    ret = CAEN_DGTZ_FreeEvent(handle[b],(void**) &Evt);
			    //delete f1;
					

		   }



  
		    c = checkCommand();
		    if (c == 1) goto Continue;
		    if (c == 2) goto Continue;
        } // end of loop on boards
    } // end of readout loop

Continue:
    for(b=0; b<MAXNB; b++){
	printf("\n|-----> Board %d: Retrieved %d Events\n",b, Nevent);
	printf("\n|----->  Output in TOSCA_%ld.root \n",now);
	}
    goto QuitProgram;



/******************************************* Quit program routine *************************************************/
QuitProgram:

    // Free the buffers and close the digitizers
	ret = CAEN_DGTZ_SWStopAcquisition(handle[b]);
	ret = CAEN_DGTZ_FreeReadoutBuffer(&buffer);
	fclose(fptr);

// Close Root file
       T->Write();
       rf->Write();
       rf->Close();
        
         

    for(b=0; b<MAXNB; b++){
        ret = CAEN_DGTZ_CloseDigitizer(handle[b]);
    }
    printf("\n|-----> Press ENTER to quit... \n");
    getchar();
    printf("\n|-----> Quitting Program... \n");
   
 }
}












