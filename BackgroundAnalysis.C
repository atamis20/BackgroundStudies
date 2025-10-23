//
//   root -l './BackgroundAnalysis.C()' for default naming
//
//    or
//
//   root -l './BackgroundAnalysis.C("DIS_file","electron beam-gas file", "proton-beam gas file")'
//

void BackgroundAnalysis(const char *name_DIS ="DIS_275_10.edm4hep.tracking_and_pfrich.root", const char *name_e = "beam_gas_e.edm4hep.tracking_and_pfrich.root", const char *name_p = "beam_gas_proton_275.edm4hep.tracking_and_pfrich.root")
{
  auto fdata = new TFile(name_DIS);
  auto fdata2 = new TFile(name_e); 
  auto fdata3 = new TFile(name_p);

  //Adjust these when you're changing your simulation files.
  //luminosity and cross section for DIS
  int luminosity =1*pow(10,7);//1 x 10^34 cm^-2 -1 = 1*10^7  Cancelling out units in cross section, so that we can enter cross-section in mb
  double DIS_crosssec = .000328369; //mb  //mb = 10^(-31) m^2 = 10^(-27) cm^2

  //Interaction rate for beam gas files (from )
  double bg_elec_rate = 3.177*1000000;
  double bg_proton_rate = 22.5*1000;



  //Feed in DIS file
  TTree *t_DIS = (TTree*) (fdata->Get("events"));
  t_DIS->SetName("events_DIS");
  auto event_DIS = new vector<edm4hep::SimTrackerHitData>;
  t_DIS->SetBranchAddress("PFRICHHits", &event_DIS);
  int nEvents_DIS = t_DIS->GetEntries();

  //Feed in electron beam+gas file
  
  TTree *t_bg_elec = (TTree*) (fdata2->Get("events"));
  t_bg_elec->SetName("events_bg_elec");
  auto event_bg_elec = new vector<edm4hep::SimTrackerHitData>;
  t_bg_elec->SetBranchAddress("PFRICHHits", &event_bg_elec);
  int nEvents_bg_elec = t_bg_elec->GetEntries();

  //Feed in proton beam_gas file
  
  TTree *t_bg_proton = (TTree*) (fdata3->Get("events"));
  t_bg_proton->SetName("events_bg_proton");
  auto event_bg_proton = new vector<edm4hep::SimTrackerHitData>;
  t_bg_proton->SetBranchAddress("PFRICHHits", &event_bg_proton);
  int nEvents_bg_proton = t_bg_proton->GetEntries();
  
  
  printf("%d event(s) total\n", nEvents_DIS);

 

  //Declare histos

  //Quantum efficiency, and its wavelength ranges, implimented as a hardcoded scaling factor (rather than actually dropping photons) for now in order to not cut out statistics unnessesarily for ease of re-running, will impliment with IRT directly in the future.
  
  vector<double> quantumefficiency = {/*160nm*/0.25,/*180nm*/0.26,/*200nm*/0.27,/*220nm*/0.30,/*240nm*/0.32,/*260nm*/0.35,/*280nm*/0.36, /*300nm*/0.36,/*320nm*/0.36, /*340nm*/0.36,  /*360nm*/0.37,/*380nm*/0.35,/*400nm*/0.30,/*420nm*/0.27,/*440nm*/0.24,/*460nm*/0.20,/*480nm*/0.18,/*500nm*/0.15,/*520nm*/0.13,/*540nm*/0.11, /*560nm*/0.10, /*580nm*/0.09,/*600nm*/0.08,/*620nm*/0.07,/*640nm*/0.05, /*660nm*/0.05};
    
  vector<double> wavelengthV = {160,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660};

  //Number of x,y bins in each direction and the spanning range, in cm
  int nBins = 200;
  int xmin = -70;
  int xmax = 70;

  //scaling factor that converts /s normalized histograms into /year.  Assuming 6 months of 24hr running.
  double scalingfactor = 1*(365/2.0)*24*60*60;


  //Decare Histograms used to record photon counts and other properties

  //DIS
  TH2D* DIS_XY = new TH2D("DIS_XY","",nBins,xmin,xmax,nBins,xmin,xmax);
  TH2D* DIS_momentumXY = new TH2D("DIS_momentumXY","",100,-0.00000002,0.000000002,100,-0.00000002,0.000000002);
  TH1D* DIS_momentumZ = new TH1D("DIS_momentumZ","",100,-0.00000001,0.000000002);
  TH1D* DIS_Energy = new TH1D("DIS_Energy","",100,0,.00000003);
  TH1D* DIS_Wavelength = new TH1D("DIS_Wavelength","",100,0,800);
  TH1D* DIS_Wavelength_eff = new TH1D("DIS_Wavelength_eff","",100,0,800);
  TH2D* DIS_WavelengthVsR_2D = new TH2D("DISWavelengthVsR_2D","",nBins/2,0,xmax,100,0,800);
  TH2D* DIS_TimeVsWavelength_2D = new TH2D("DISTimeVsWavelength_2D","",100,0,100,100,0,800);
  TH1D* DIS_WavelengthVsR = new TH1D("DISWavelengthVsR","",nBins/2,0,xmax);
  TH1D* DIS_TimeVsR = new TH1D("DISTimeVsR","",nBins/2,0,xmax);
  TH2D* DIS_TimeVsR_2D = new TH2D("DISTimeVsR_2D","",nBins/2,0,xmax,100,0,100);
  TH1D* DIS_CountsVsR = new TH1D("DISCountsVsR","",nBins/2,0,xmax);
  TH1D* DIS_TimeHist = new TH1D("DISTimeHist","",100,0,100);
  TH1D* DIS_TimeVsWavelength = new TH1D("DISTimeVsWavelength","",100,0,100);
  
  //bg_elec
  TH2D* bg_elec_XY = new TH2D("bg_elec_XY","",nBins,xmin,xmax,nBins,xmin,xmax);
  TH2D* bg_elec_momentumXY = new TH2D("bg_elec_momentumXY","",100,-0.00000002,0.000000002,100,-0.00000002,0.000000002);
  TH1D* bg_elec_momentumZ = new TH1D("bg_elec_momentumZ","",100,-0.00000001,0.000000002);
  TH1D* bg_elec_Energy = new TH1D("bg_elec_Energy","",100,0,.00000003);
  TH1D* bg_elec_Wavelength = new TH1D("bg_elec_Wavelength","",100,0,800);
  TH1D* bg_elec_Wavelength_eff = new TH1D("bg_elec_Wavelength_eff","",100,0,800);
  TH2D* bg_elec_WavelengthVsR_2D = new TH2D("bg_elecWavelengthVsR_2D","",nBins/2,0,xmax,100,0,800);
  TH2D* bg_elec_TimeVsWavelength_2D = new TH2D("bg_elecTimeVsWavelength_2D","",100,0,100,100,0,800);
  TH1D* bg_elec_WavelengthVsR = new TH1D("bg_elecWavelengthVsR","",nBins/2,0,xmax);
  TH1D* bg_elec_TimeVsR = new TH1D("bg_elecTimeVsR","",nBins/2,0,xmax);
  TH2D* bg_elec_TimeVsR_2D = new TH2D("bg_elecTimeVsR_2D","",nBins/2,0,xmax,100,0,100);
  TH1D* bg_elec_CountsVsR = new TH1D("bg_elecCountsVsR","",nBins/2,0,xmax);
  TH1D* bg_elec_TimeHist = new TH1D("bg_elecTimeHist","",100,0,100);
  TH1D* bg_elec_TimeVsWavelength = new TH1D("bg_elecTimeVsWavelength","",100,0,100);
  
  //bg_proton
  TH2D* bg_proton_XY = new TH2D("bg_proton_XY","",nBins,xmin,xmax,nBins,xmin,xmax);
  TH2D* bg_proton_momentumXY = new TH2D("bg_proton_momentumXY","",100,-0.00000002,0.000000002,100,-0.00000002,0.000000002);
  TH1D* bg_proton_momentumZ = new TH1D("bg_proton_momentumZ","",100,-0.00000001,0.000000002);
  TH1D* bg_proton_Energy = new TH1D("bg_proton_Energy","",100,0,.00000003);
  TH1D* bg_proton_Wavelength = new TH1D("bg_proton_Wavelength","",100,0,800);
  TH1D* bg_proton_Wavelength_eff = new TH1D("bg_proton_Wavelength_eff","",100,0,800);
  TH2D* bg_proton_WavelengthVsR_2D = new TH2D("bg_protonWavelengthVsR_2D","",nBins/2,0,xmax,100,0,800);
  TH2D* bg_proton_TimeVsWavelength_2D = new TH2D("bg_protonTimeVsWavelength_2D","",100,0,100,100,0,800);
  TH1D* bg_proton_WavelengthVsR = new TH1D("bg_protonWavelengthVsR","",nBins/2,0,xmax);
  TH1D* bg_proton_TimeVsR = new TH1D("bg_protonTimeVsR","",nBins/2,0,xmax);
  TH2D* bg_proton_TimeVsR_2D = new TH2D("bg_protonTimeVsR_2D","",nBins/2,0,xmax,100,0,100);
  TH1D* bg_proton_CountsVsR = new TH1D("bg_protonCountsVsR","",nBins/2,0,xmax);
  TH1D* bg_proton_TimeHist = new TH1D("bg_protonTimeHist","",100,0,100);
  TH1D* bg_proton_TimeVsWavelength = new TH1D("bg_protonTimeVsWavelength","",100,0,100);
  

  //total collection
  TH2D* total_XY = new TH2D("total_XY","",nBins,xmin,xmax,nBins,xmin,xmax);

  
  //Normalize histogram by bin area and cross section, into units of cm^(-2)s^(-1)
  double binScaling = 1/(pow(((xmax-xmin)/((double)nBins)),2));
  double eventweight_DIS = (luminosity*DIS_crosssec/nEvents_DIS);

  //Analyze DIS events
  for(int i = 1; i < nEvents_DIS; i++){

    
    if ( !(i%((nEvents_DIS/20))) ) cout << "DIS: " << 100*(double)i/(double)nEvents_DIS << " percent done "<< endl;

    t_DIS->GetEntry(i);
    int Hits = event_DIS->size();
    
    for(int j = 1; j < Hits; j++){
      auto xsim = event_DIS->at(j);
      auto xvec = xsim.position.x;
      auto yvec = xsim.position.y;
      auto xmom = xsim.momentum.x;
      auto ymom = xsim.momentum.y;
      auto zmom = xsim.momentum.z;
      auto time = xsim.time;
      auto energy = sqrt(pow(xmom,2)+pow(ymom,2)+pow(zmom,2));
      auto wavelength = (4.136*3*pow(10,-7))/(energy);
      auto R = sqrt(pow(xvec/10,2)+pow(yvec/10,2));
      if(wavelength < 2){continue;}

      //determine efficiency to use dependenting on wavelength
      double eff = 0.05;
      for(int wav = 0; wav < wavelengthV.size()-1; wav++){
	if((wavelength > wavelengthV[wav]) && (wavelength < wavelengthV[wav+1])){eff=quantumefficiency[wav];}
      }
      eff = eff*0.81;
      
      
      DIS_WavelengthVsR->Fill(R,wavelength*eventweight_DIS*eff);
      DIS_WavelengthVsR_2D->Fill(R,wavelength,eventweight_DIS*eff);
      DIS_TimeVsR->Fill(R,time*eventweight_DIS*eff);
      DIS_TimeVsR_2D->Fill(R,time,eventweight_DIS*eff);
      DIS_TimeHist->Fill(time,eventweight_DIS*eff);
      DIS_TimeVsWavelength->Fill(time,wavelength*eventweight_DIS*eff);
      DIS_TimeVsWavelength_2D->Fill(time,wavelength,eventweight_DIS*eff);
      DIS_CountsVsR->Fill(R,eventweight_DIS*eff);
      DIS_momentumXY->Fill(xmom,ymom,eventweight_DIS*eff);
      DIS_momentumZ->Fill(zmom,eventweight_DIS*eff);
      DIS_Energy->Fill(energy,eventweight_DIS*eff);
      DIS_Wavelength->Fill(wavelength,eventweight_DIS);
      DIS_Wavelength_eff->Fill(wavelength,eventweight_DIS*eff);
      DIS_XY->Fill(xvec/10,yvec/10,eventweight_DIS*binScaling*eff); 
    }
  }
  
  DIS_TimeVsWavelength->Divide(DIS_TimeHist);
  DIS_WavelengthVsR->Divide(DIS_CountsVsR);
  DIS_TimeVsR->Divide(DIS_CountsVsR);


  double eventweight_bg_elec = (bg_elec_rate/(nEvents_bg_elec));

  //analyze electron beamgas events
  for(int i = 1; i < nEvents_bg_elec; i++){

    
    if ( !(i%(nEvents_bg_elec/20)) ) cout << "electron beam gas: " << 100*(double)i/(double)nEvents_bg_elec << " percent done "<< endl;

    
    t_bg_elec->GetEntry(i);
    int Hits = event_bg_elec->size();
    
    for(int j = 1; j < Hits; j++){

      auto xsim = event_bg_elec->at(j);
      auto xvec = xsim.position.x;
      auto yvec = xsim.position.y;
      auto xmom = xsim.momentum.x;
      auto ymom = xsim.momentum.y;
      auto zmom = xsim.momentum.z;
      auto time = xsim.time;
  
      auto energy = sqrt(pow(xmom,2)+pow(ymom,2)+pow(zmom,2));
      auto wavelength = (4.136*3*pow(10,-7))/(energy);
      auto R = sqrt(pow(xvec/10,2)+pow(yvec/10,2));
      if(wavelength < 2){continue;}
      
      double eff = 0.05;
      for(int wav = 0; wav < wavelengthV.size()-1; wav++){
	if((wavelength > wavelengthV[wav]) && (wavelength < wavelengthV[wav+1])){eff=quantumefficiency[wav];}
      }
      eff = eff*0.81;
    
      bg_elec_WavelengthVsR->Fill(R,wavelength*eventweight_bg_elec*eff);
      bg_elec_WavelengthVsR_2D->Fill(R,wavelength,eventweight_bg_elec*eff);
      bg_elec_TimeVsR->Fill(R,time*eventweight_bg_elec*eff);
      bg_elec_TimeVsR_2D->Fill(R,time,eventweight_bg_elec*eff);
      bg_elec_TimeHist->Fill(time,eventweight_bg_elec*eff);
      bg_elec_TimeVsWavelength->Fill(time,wavelength*eventweight_bg_elec*eff);
      bg_elec_TimeVsWavelength_2D->Fill(time,wavelength,eventweight_bg_elec*eff);
      bg_elec_CountsVsR->Fill(R,eventweight_bg_elec*eff);
      bg_elec_momentumXY->Fill(xmom,ymom,eventweight_bg_elec*eff);
      bg_elec_momentumZ->Fill(zmom,eventweight_bg_elec*eff);
      bg_elec_Energy->Fill(energy,eventweight_bg_elec*eff);
      bg_elec_Wavelength->Fill(wavelength,eventweight_bg_elec);
      bg_elec_Wavelength_eff->Fill(wavelength,eventweight_bg_elec*eff);
      bg_elec_XY->Fill(xvec/10,yvec/10,eventweight_bg_elec*binScaling*eff); 
    }
  }
  
  bg_elec_TimeVsWavelength->Divide(bg_elec_TimeHist);
  bg_elec_WavelengthVsR->Divide(bg_elec_CountsVsR);
  bg_elec_TimeVsR->Divide(bg_elec_CountsVsR);



  
  printf("%d event(s) total\n", nEvents_bg_proton);

  
  double eventweight_bg_proton = (bg_proton_rate/nEvents_bg_proton);
  for(int i = 1; i < nEvents_bg_proton; i++){
    
    if ( !(i%(nEvents_bg_proton/20)) ) cout << "proton beam gas: " << 100*(double)i/(double)nEvents_bg_proton << " percent done "<< endl;
    
    t_bg_proton->GetEntry(i);
    int Hits = event_bg_proton->size();
    
    for(int j = 1; j < Hits; j++){
      auto xsim = event_bg_proton->at(j);
      auto xvec = xsim.position.x;
      auto yvec = xsim.position.y;
      auto xmom = xsim.momentum.x;
      auto ymom = xsim.momentum.y;
      auto zmom = xsim.momentum.z;
      auto time = xsim.time;

      auto energy = sqrt(pow(xmom,2)+pow(ymom,2)+pow(zmom,2));
      auto wavelength = (4.136*3*pow(10,-7))/(energy);
      auto R = sqrt(pow(xvec/10,2)+pow(yvec/10,2));
      if(wavelength < 2){continue;}
      
      
      double eff = 0.05;
      for(int wav = 0; wav < wavelengthV.size()-1; wav++){
	if((wavelength > wavelengthV[wav]) && (wavelength < wavelengthV[wav+1])){eff=quantumefficiency[wav];}
      }
      eff = eff*0.81;

      
      bg_proton_WavelengthVsR->Fill(R,wavelength*eventweight_bg_proton*eff);
      bg_proton_WavelengthVsR_2D->Fill(R,wavelength,eventweight_bg_proton*eff);
      bg_proton_TimeVsR->Fill(R,time*eventweight_bg_proton*eff);
      bg_proton_TimeVsR_2D->Fill(R,time,eventweight_bg_proton*eff);
      bg_proton_TimeHist->Fill(time,eventweight_bg_proton*eff);
      bg_proton_TimeVsWavelength->Fill(time,wavelength*eventweight_bg_proton*eff);
      bg_proton_TimeVsWavelength_2D->Fill(time,wavelength,eventweight_bg_proton*eff);
      bg_proton_CountsVsR->Fill(R,eventweight_bg_proton*eff);
      bg_proton_momentumXY->Fill(xmom,ymom,eventweight_bg_proton*eff);
      bg_proton_momentumZ->Fill(zmom,eventweight_bg_proton*eff);
      bg_proton_Energy->Fill(energy,eventweight_bg_proton*eff);
      bg_proton_Wavelength->Fill(wavelength,eventweight_bg_proton);
      bg_proton_Wavelength_eff->Fill(wavelength,eventweight_bg_proton*eff);
      bg_proton_XY->Fill(xvec/10,yvec/10,eventweight_bg_proton*binScaling*eff); 
    }
  }
  bg_proton_TimeVsWavelength->Divide(bg_proton_TimeHist);
  bg_proton_WavelengthVsR->Divide(bg_proton_CountsVsR);
  bg_proton_TimeVsR->Divide(bg_proton_CountsVsR);
  
  gStyle->SetOptStat(0);

  
  gStyle->SetTitleSize(0.2,"t");
  TEllipse* InnerR = new TEllipse(0,0,9.5,9.5);
  TEllipse* OuterR = new TEllipse(0,0,65,65);
  InnerR->SetLineColor(kRed);
  OuterR->SetLineColor(kRed);
  InnerR->SetLineWidth(2);
  OuterR->SetLineWidth(2);
  InnerR->SetFillStyle(0);
  OuterR->SetFillStyle(0);
  TCanvas* RadHardC_XY_DIS = new TCanvas("RadHardC_XY_DIS","RadHardC_XY_DIS",50,50,1000,1000);
  RadHardC_XY_DIS->SetBottomMargin(0.15);
  RadHardC_XY_DIS->SetTopMargin(0.15);
  RadHardC_XY_DIS->SetRightMargin(0.15);
  RadHardC_XY_DIS->SetLeftMargin(0.15);
  RadHardC_XY_DIS->SetLogz();
  //DIS_XY->SetTitle("Hi");
  DIS_XY->SetTitle("10x275GeV e+p DIS");
  DIS_XY->GetZaxis()->SetTitle("Photon Rate [cm^{-2}s^{-1}]");
  DIS_XY->GetZaxis()->SetRangeUser(1,100000);
  DIS_XY->GetZaxis()->SetLabelSize(0.02);
  DIS_XY->GetXaxis()->SetTitle("X [cm]");
  DIS_XY->GetYaxis()->SetTitle("Y [cm]");
  
  DIS_XY->Draw("colz");
  // OuterR->Draw();
  //InnerR->Draw();
  
  RadHardC_XY_DIS->SaveAs("Plots/DIS/PhotonRate_DIS.pdf");

    TCanvas* DIS_MomentumXYC = new TCanvas("DIS_MomentumXYC","MomentumXYC");
  DIS_MomentumXYC->SetLogz();
  DIS_momentumXY->GetXaxis()->SetTitle("DIS_x Momentum [GeV]");
  DIS_momentumXY->GetYaxis()->SetTitle("DIS_y Momentum [GeV]");
  DIS_momentumXY->Draw("colz");

  TCanvas* DIS_EnergyC = new TCanvas("DIS_EnergyC","EnergyC");
  DIS_EnergyC->SetLogy();
  DIS_Energy->GetXaxis()->SetTitle("Photon Energy [GeV]");
  DIS_Energy->GetYaxis()->SetTitle("Rate [s^{-1}]");
  DIS_Energy->Draw("p E1");
  DIS_EnergyC->SaveAs("Plots/DIS/PhotonEnergy.pdf");

  TCanvas* DIS_WavelengthC = new TCanvas("DIS_WavelengthC","WavelengthC");
  DIS_WavelengthC->SetLogy();
  DIS_Wavelength->GetXaxis()->SetTitle("Photon Wavelength [nm]");
  DIS_Wavelength->GetYaxis()->SetTitle("Rate [s^{-1}]");
  DIS_Wavelength->Draw("p E1");
  DIS_Wavelength_eff->SetLineColor(kRed);
  DIS_Wavelength_eff->Draw("p E1 same");
  DIS_WavelengthC->SaveAs("Plots/DIS/PhotonWavelength_DIS.pdf");

  TLegend* EffL = new TLegend();
  EffL->SetBorderSize(0);
  EffL->AddEntry(DIS_Wavelength,"No Quantum Efficiency");
  EffL->AddEntry(DIS_Wavelength_eff,"Quantum Efficiency Applied");
  EffL->Draw();

  TCanvas* DIS_TimeVsR_2DC = new TCanvas("DIS_Time Vs R 2D","Time Vs R 2D");
  DIS_TimeVsR_2DC->SetLogz();
  DIS_TimeVsR_2D->GetXaxis()->SetTitle("Time");
  DIS_TimeVsR_2D->GetYaxis()->SetTitle("Radius [cm]");
  DIS_TimeVsR_2D->Draw("colz");
  DIS_TimeVsR_2D->SaveAs("Plots/DIS/TimeVsR_DIS.pdf");

  TCanvas* DIS_TimeC = new TCanvas("DIS_Time","Time");
  DIS_TimeHist->GetXaxis()->SetTitle("Time");
  DIS_TimeHist->GetYaxis()->SetTitle("Rate [s^{-1}]");
  DIS_TimeHist->Draw("p E1");
  DIS_TimeC->SaveAs("Plots/DIS/ArrivalTime_DIS.pdf");
  
  TCanvas* DIS_TimeVsWavelengthC = new TCanvas("DIS_Time Vs Wavelength","Time Vs Wavelength");
  DIS_TimeVsWavelength->GetXaxis()->SetTitle("Time");
  DIS_TimeVsWavelength->GetYaxis()->SetTitle("Average Wavelength [nm]");
  DIS_TimeVsWavelength->Draw("p E1");
  DIS_TimeVsWavelength->SaveAs("Plots/DIS/TimeVsWavelength_DIS.pdf");
  
  TCanvas* DIS_TimeVsWavelength_2DC = new TCanvas("DIS_Time Vs Wavelength 2D","Time Vs Wavelength 2D");
  DIS_TimeVsWavelength_2DC->SetLogz();
  DIS_TimeVsWavelength_2D->GetXaxis()->SetTitle("Time");
  DIS_TimeVsWavelength_2D->GetYaxis()->SetTitle("Photon Wavelength [nm]");
  DIS_TimeVsWavelength_2D->Draw("colz");
  DIS_TimeVsWavelength_2DC->SaveAs("Plots/DIS/TimeVsWavelength_DIS.pdf");
  
  TCanvas* DIS_WavelengthVsRC = new TCanvas("DIS_Wavelength Vs R","Wavelength Vs R");
  DIS_WavelengthVsR->GetXaxis()->SetTitle("Radius [cm]");
  DIS_WavelengthVsR->GetYaxis()->SetTitle("Average Wavelength [nm]");
  DIS_WavelengthVsR->GetYaxis()->SetRangeUser(320,440);
  DIS_WavelengthVsR->Draw("p E1");
  DIS_WavelengthVsRC->SaveAs("Plots/DIS/WavelengthVsR.pdf");
  
  TCanvas* DIS_WavelengthVsR_2DC = new TCanvas("DIS_Wavelength Vs R 2D","Wavelength Vs R 2D");
  DIS_WavelengthVsR_2DC->SetLogz();
  DIS_WavelengthVsR_2D->GetXaxis()->SetTitle("Radius [cm]");
  DIS_WavelengthVsR_2D->GetYaxis()->SetTitle("Photon Wavelength [nm]");
  DIS_WavelengthVsR_2D->Draw("colz");
  DIS_WavelengthVsR_2DC->SaveAs("Plots/DIS/WavelengthVsR_2D.pdf");


  
  //beam+gas elec plotting
  TCanvas* RadHardC_XY_bg_elec = new TCanvas("RadHardC_XY_bg_elec","RadHardC_XY_bg_elec",50,50,1000,1000);
  RadHardC_XY_bg_elec->SetBottomMargin(0.15);
  RadHardC_XY_bg_elec->SetTopMargin(0.15);
  RadHardC_XY_bg_elec->SetRightMargin(0.15);
  RadHardC_XY_bg_elec->SetLeftMargin(0.15);
  RadHardC_XY_bg_elec->SetLogz();
  //bg_elec_XY->SetTitle("Hi");
  bg_elec_XY->SetTitle("10GeV electron beam+gas");
  bg_elec_XY->GetZaxis()->SetTitle("Photon Rate [cm^{-2}s^{-1}]");
  bg_elec_XY->GetZaxis()->SetRangeUser(1,100000);
  bg_elec_XY->GetZaxis()->SetLabelSize(0.02);
  bg_elec_XY->GetXaxis()->SetTitle("X [cm]");
  bg_elec_XY->GetYaxis()->SetTitle("Y [cm]");
  bg_elec_XY->Draw("colz");

  RadHardC_XY_bg_elec->SaveAs("Plots/bg_elec/PhotonRate_bg_elec.pdf");

  TCanvas* bg_elec_MomentumXYC = new TCanvas("bg_elec_MomentumXYC","MomentumXYC");
  bg_elec_MomentumXYC->SetLogz();
  bg_elec_momentumXY->GetXaxis()->SetTitle("bg_elec_x Momentum [GeV]");
  bg_elec_momentumXY->GetYaxis()->SetTitle("bg_elec_y Momentum [GeV]");
  bg_elec_momentumXY->Draw("colz");

  TCanvas* bg_elec_EnergyC = new TCanvas("bg_elec_EnergyC","EnergyC");
  bg_elec_EnergyC->SetLogy();
  bg_elec_Energy->GetXaxis()->SetTitle("Photon Energy [GeV]");
  bg_elec_Energy->GetYaxis()->SetTitle("Rate [s^{-1}]");
  bg_elec_Energy->Draw("p E1");
  bg_elec_EnergyC->SaveAs("Plots/bg_elec/PhotonEnergy.pdf");

  TCanvas* bg_elec_WavelengthC = new TCanvas("bg_elec_WavelengthC","WavelengthC");
  bg_elec_WavelengthC->SetLogy();
  bg_elec_Wavelength->GetXaxis()->SetTitle("Photon Wavelength [nm]");
  bg_elec_Wavelength->GetYaxis()->SetTitle("Rate [s^{-1}]");
  bg_elec_Wavelength->Draw("p E1");
  bg_elec_Wavelength_eff->SetLineColor(kRed);
  bg_elec_Wavelength_eff->Draw("p E1 same");
  bg_elec_WavelengthC->SaveAs("Plots/bg_elec/PhotonWavelength_bg_elec.pdf");

  TCanvas* bg_elec_TimeVsR_2DC = new TCanvas("bg_elec_Time Vs R 2D","Time Vs R 2D");
  bg_elec_TimeVsR_2DC->SetLogz();
  bg_elec_TimeVsR_2D->GetXaxis()->SetTitle("Time");
  bg_elec_TimeVsR_2D->GetYaxis()->SetTitle("Radius [cm]");
  bg_elec_TimeVsR_2D->Draw("colz");
  bg_elec_TimeVsR_2D->SaveAs("Plots/bg_elec/TimeVsR_bg_elec.pdf");

  TCanvas* bg_elec_TimeC = new TCanvas("bg_elec_Time","Time");
  bg_elec_TimeHist->GetXaxis()->SetTitle("Time");
  bg_elec_TimeHist->GetYaxis()->SetTitle("Rate [s^{-1}]");
  bg_elec_TimeHist->Draw("p E1");
  bg_elec_TimeC->SaveAs("Plots/bg_elec/ArrivalTime_bg_elec.pdf");
  
  TCanvas* bg_elec_TimeVsWavelengthC = new TCanvas("bg_elec_Time Vs Wavelength","Time Vs Wavelength");
  bg_elec_TimeVsWavelength->GetXaxis()->SetTitle("Time");
  bg_elec_TimeVsWavelength->GetYaxis()->SetTitle("Average Wavelength [nm]");
  bg_elec_TimeVsWavelength->Draw("p E1");
  bg_elec_TimeVsWavelength->SaveAs("Plots/bg_elec/TimeVsWavelength_bg_elec.pdf");
  
  TCanvas* bg_elec_TimeVsWavelength_2DC = new TCanvas("bg_elec_Time Vs Wavelength 2D","Time Vs Wavelength 2D");
  bg_elec_TimeVsWavelength_2DC->SetLogz();
  bg_elec_TimeVsWavelength_2D->GetXaxis()->SetTitle("Time");
  bg_elec_TimeVsWavelength_2D->GetYaxis()->SetTitle("Photon Wavelength [nm]");
  bg_elec_TimeVsWavelength_2D->Draw("colz");
  bg_elec_TimeVsWavelength_2DC->SaveAs("Plots/bg_elec/TimeVsWavelength_bg_elec.pdf");
  
  TCanvas* bg_elec_WavelengthVsRC = new TCanvas("bg_elec_Wavelength Vs R","Wavelength Vs R");
  bg_elec_WavelengthVsR->GetXaxis()->SetTitle("Radius [cm]");
  bg_elec_WavelengthVsR->GetYaxis()->SetTitle("Average Wavelength [nm]");
  bg_elec_WavelengthVsR->GetYaxis()->SetRangeUser(320,440);
  bg_elec_WavelengthVsR->Draw("p E1");
  bg_elec_WavelengthVsRC->SaveAs("Plots/bg_elec/WavelengthVsR.pdf");
  
  TCanvas* bg_elec_WavelengthVsR_2DC = new TCanvas("bg_elec_Wavelength Vs R 2D","Wavelength Vs R 2D");
  bg_elec_WavelengthVsR_2DC->SetLogz();
  bg_elec_WavelengthVsR_2D->GetXaxis()->SetTitle("Radius [cm]");
  bg_elec_WavelengthVsR_2D->GetYaxis()->SetTitle("Photon Wavelength [nm]");
  bg_elec_WavelengthVsR_2D->Draw("colz");
  bg_elec_WavelengthVsR_2DC->SaveAs("Plots/bg_elec/WavelengthVsR_2D.pdf");
  
  
  
  //OuterR->Draw();
  //InnerR->Draw();

  
  TCanvas* RadHardC_XY_bg_proton = new TCanvas("RadHardC_XY_bg_proton","RadHardC_XY_bg_proton",50,50,1000,1000);
  RadHardC_XY_bg_proton->SetBottomMargin(0.15);
  RadHardC_XY_bg_proton->SetTopMargin(0.15);
  RadHardC_XY_bg_proton->SetRightMargin(0.15);
  RadHardC_XY_bg_proton->SetLeftMargin(0.15);
  RadHardC_XY_bg_proton->SetLogz();
  //bg_proton_XY->SetTitle("Hi");
  bg_proton_XY->SetTitle("275GeV proton beam+gas");
  bg_proton_XY->GetZaxis()->SetTitle("Photon Rate [cm^{-2}s^{-1}]");
  //bg_proton_XY->GetZaxis()->SetRangeUser(100,10000000);
  bg_proton_XY->GetZaxis()->SetLabelSize(0.02);
  bg_proton_XY->GetXaxis()->SetTitle("X [cm]");
  bg_proton_XY->GetYaxis()->SetTitle("Y [cm]");

  
  bg_proton_XY->GetZaxis()->SetRangeUser(1,100000);
  bg_proton_XY->Draw("colz");
  //OuterR->Draw();
  //InnerR->Draw();
  RadHardC_XY_bg_proton->SaveAs("Plots/bg_proton/PhotonRate_bg_proton.pdf");

    TCanvas* bg_proton_MomentumXYC = new TCanvas("bg_proton_MomentumXYC","MomentumXYC");
  bg_proton_MomentumXYC->SetLogz();
  bg_proton_momentumXY->GetXaxis()->SetTitle("bg_proton_x Momentum [GeV]");
  bg_proton_momentumXY->GetYaxis()->SetTitle("bg_proton_y Momentum [GeV]");
  bg_proton_momentumXY->Draw("colz");

  TCanvas* bg_proton_EnergyC = new TCanvas("bg_proton_EnergyC","EnergyC");
  bg_proton_EnergyC->SetLogy();
  bg_proton_Energy->GetXaxis()->SetTitle("Photon Energy [GeV]");
  bg_proton_Energy->GetYaxis()->SetTitle("Rate [s^{-1}]");
  bg_proton_Energy->Draw("p E1");
  bg_proton_EnergyC->SaveAs("Plots/bg_proton/PhotonEnergy.pdf");

  TCanvas* bg_proton_WavelengthC = new TCanvas("bg_proton_WavelengthC","WavelengthC");
  bg_proton_WavelengthC->SetLogy();
  bg_proton_Wavelength->GetXaxis()->SetTitle("Photon Wavelength [nm]");
  bg_proton_Wavelength->GetYaxis()->SetTitle("Rate [s^{-1}]");
  bg_proton_Wavelength->Draw("p E1");
  bg_proton_Wavelength_eff->SetLineColor(kRed);
  bg_proton_Wavelength_eff->Draw("p E1 same");
  bg_proton_WavelengthC->SaveAs("Plots/bg_proton/PhotonWavelength_bg_proton.pdf");

  TCanvas* bg_proton_TimeVsR_2DC = new TCanvas("bg_proton_Time Vs R 2D","Time Vs R 2D");
  bg_proton_TimeVsR_2D->GetXaxis()->SetTitle("Time");
  bg_proton_TimeVsR_2D->GetYaxis()->SetTitle("Radius [cm]");
  bg_proton_TimeVsR_2D->Draw("colz");
  bg_proton_TimeVsR_2D->SaveAs("Plots/bg_proton/TimeVsR_bg_proton.pdf");

  TCanvas* bg_proton_TimeC = new TCanvas("bg_proton_Time","Time");
  bg_proton_TimeHist->GetXaxis()->SetTitle("Time");
  bg_proton_TimeHist->GetYaxis()->SetTitle("Rate [s^{-1}]");
  bg_proton_TimeHist->Draw("p E1");
  bg_proton_TimeC->SaveAs("Plots/bg_proton/ArrivalTime_bg_proton.pdf");
  
  TCanvas* bg_proton_TimeVsWavelengthC = new TCanvas("bg_proton_Time Vs Wavelength","Time Vs Wavelength");
  bg_proton_TimeVsWavelength->GetXaxis()->SetTitle("Time");
  bg_proton_TimeVsWavelength->GetYaxis()->SetTitle("Average Wavelength [nm]");
  bg_proton_TimeVsWavelength->Draw("p E1");
  bg_proton_TimeVsWavelength->SaveAs("Plots/bg_proton/TimeVsWavelength_bg_proton.pdf");
  
  TCanvas* bg_proton_TimeVsWavelength_2DC = new TCanvas("bg_proton_Time Vs Wavelength 2D","Time Vs Wavelength 2D");
  bg_proton_TimeVsWavelength_2D->GetXaxis()->SetTitle("Time");
  bg_proton_TimeVsWavelength_2D->GetYaxis()->SetTitle("Photon Wavelength [nm]");
  bg_proton_TimeVsWavelength_2D->Draw("colz");
  bg_proton_TimeVsWavelength_2DC->SaveAs("Plots/bg_proton/TimeVsWavelength_bg_proton.pdf");
  
  TCanvas* bg_proton_WavelengthVsRC = new TCanvas("bg_proton_Wavelength Vs R","Wavelength Vs R");
  bg_proton_WavelengthVsR->GetXaxis()->SetTitle("Radius [cm]");
  bg_proton_WavelengthVsR->GetYaxis()->SetTitle("Average Wavelength [nm]");
  bg_proton_WavelengthVsR->GetYaxis()->SetRangeUser(320,440);
  bg_proton_WavelengthVsR->Draw("p E1");
  bg_proton_WavelengthVsRC->SaveAs("Plots/bg_proton/WavelengthVsR.pdf");
  
  TCanvas* bg_proton_WavelengthVsR_2DC = new TCanvas("bg_proton_Wavelength Vs R 2D","Wavelength Vs R 2D");
  bg_proton_WavelengthVsR_2D->GetXaxis()->SetTitle("Radius [cm]");
  bg_proton_WavelengthVsR_2D->GetYaxis()->SetTitle("Photon Wavelength [nm]");
  bg_proton_WavelengthVsR_2D->Draw("colz");
  bg_proton_WavelengthVsR_2DC->SaveAs("Plots/bg_proton/WavelengthVsR_2D.pdf");

  
  total_XY->Add(DIS_XY);
  total_XY->Add(bg_elec_XY);
  total_XY->Add(bg_proton_XY);
  
  TCanvas* Totalc = new TCanvas("Totalc","Totalc",50,50,1000,1000);
  Totalc->SetBottomMargin(0.15);
  Totalc->SetTopMargin(0.15);
  Totalc->SetRightMargin(0.15);
  Totalc->SetLeftMargin(0.15);
  Totalc->SetLogz();
  //total_XY->SetTitle("Hi");
  total_XY->SetTitle("");
  total_XY->GetZaxis()->SetTitle("Photon Rate [cm^{-2}s^{-1}]");
  total_XY->GetZaxis()->SetRangeUser(1,100000);
  total_XY->GetZaxis()->SetLabelSize(0.02);
  total_XY->GetXaxis()->SetTitle("X [cm]");
  total_XY->GetYaxis()->SetTitle("Y [cm]");
  
  total_XY->Draw("colz");

  Totalc->SaveAs("Plots/Total/PhotonRate_Total.pdf");
  //OuterR->Draw();
  //InnerR->Draw();


  //Scale to one year
  TH2D* total_XY_OneYear = (TH2D*) total_XY->Clone("total_XY_OneYear");
  total_XY_OneYear->Scale(scalingfactor);
  
  TCanvas* Total_OneYearc = new TCanvas("Total_OneYearc","Total_OneYearc",50,50,1000,1000);
  Total_OneYearc->SetBottomMargin(0.15);
  Total_OneYearc->SetTopMargin(0.15);
  Total_OneYearc->SetRightMargin(0.15);
  Total_OneYearc->SetLeftMargin(0.15);
  Total_OneYearc->SetLogz();
  //total_XY_OneYear->SetTitle("Hi");
  total_XY_OneYear->SetTitle("");
  total_XY_OneYear->GetZaxis()->SetTitle("Photons per Year [cm^{-2}]");
  total_XY_OneYear->GetZaxis()->SetRangeUser(10000000,1000000000000);
  total_XY_OneYear->GetZaxis()->SetLabelSize(0.02);
  total_XY_OneYear->GetXaxis()->SetTitle("X [cm]");
  total_XY_OneYear->GetYaxis()->SetTitle("Y [cm]");
  
  total_XY_OneYear->Draw("colz");

  Total_OneYearc->SaveAs("Plots/Total/PhotonRate_OneYear_Total.pdf");

  TCanvas* Total_WavelengthC = new TCanvas("Total_WavelengthC","WavelengthC");
  Total_WavelengthC->SetLogy();
  DIS_Wavelength->GetXaxis()->SetTitle("Time");
  DIS_Wavelength->GetYaxis()->SetTitle("Rate [s^{-1}]");
  DIS_Wavelength->SetLineColor(kBlue);
  DIS_Wavelength->GetYaxis()->SetRangeUser(1000,1000000);
  DIS_Wavelength->Draw("p E1");
  bg_elec_Wavelength->SetLineColor(kBlack);
  bg_elec_Wavelength->Draw("p E1 same");
  bg_proton_Wavelength->SetLineColor(kRed);
  bg_proton_Wavelength->Draw("p E1 same");

  
  TLegend* TypeL = new TLegend();
  TypeL->SetBorderSize(0);
  TypeL->AddEntry(DIS_Wavelength,"DIS");
  TypeL->AddEntry(bg_elec_Wavelength,"electron beam gas");
  TypeL->AddEntry(bg_proton_Wavelength,"proton beam gas");
  TypeL->Draw();
  
  Total_WavelengthC->SaveAs("Plots/Total/PhotonWavelength_Total.pdf");

  
  TCanvas* Total_TimeC = new TCanvas("Total_Time","Time");
  DIS_TimeHist->GetXaxis()->SetTitle("Time");
  DIS_TimeHist->GetYaxis()->SetTitle("Rate [s^{-1}]");
  DIS_TimeHist->SetLineColor(kBlue);
  DIS_TimeHist->Draw("p E1");
  bg_elec_TimeHist->SetLineColor(kBlack);
  bg_elec_TimeHist->Draw("p E1 same");
  bg_proton_TimeHist->SetLineColor(kRed);
  bg_proton_TimeHist->Draw("p E1 same");
  TypeL->Draw();
  Total_TimeC->SaveAs("Plots/Total/ArrivalTime_Total.pdf");
  
  
  TCanvas* Total_WavelengthVsRC = new TCanvas("Total_Wavelength Vs R","Wavelength Vs R");
  DIS_WavelengthVsR->GetXaxis()->SetTitle("Radius [cm]");
  DIS_WavelengthVsR->GetYaxis()->SetTitle("Average Wavelength [nm]");
  DIS_WavelengthVsR->GetYaxis()->SetRangeUser(320,440);
  DIS_WavelengthVsR->SetLineColor(kBlue);
  bg_elec_WavelengthVsR->SetLineColor(kBlack);
  bg_proton_WavelengthVsR->SetLineColor(kRed);
  
  DIS_WavelengthVsR->Draw("p E1");
  bg_elec_WavelengthVsR->Draw("p E1 same");
  bg_proton_WavelengthVsR->Draw("p E1 same");
  
  Total_WavelengthVsRC->SaveAs("Plots/Total/WavelengthVsR.pdf");
  TypeL->Draw();

  TCanvas* Total_MomentumZC = new TCanvas("MomentumZC","MomentumZC");
  Total_MomentumZC->SetLogy();
  DIS_momentumZ->GetXaxis()->SetTitle("Z Momentum [eV]");
  DIS_momentumZ->GetYaxis()->SetTitle("Rate [s^{-1}]");
  DIS_momentumZ->SetLineColor(kBlue);
  DIS_momentumZ->GetYaxis()->SetRangeUser(10,100000);
  DIS_momentumZ->Draw("p E1");
  bg_elec_momentumZ->SetLineColor(kBlack);
  bg_elec_momentumZ->Draw("p E1 same");
  bg_proton_momentumZ->SetLineColor(kRed);
  bg_proton_momentumZ->Draw("p E1 same");
  Total_MomentumZC->SaveAs("Plots/Total/Total_MomentumZ.pdf");
  TypeL->Draw();

  

  
  
} // frich_hit_map()
