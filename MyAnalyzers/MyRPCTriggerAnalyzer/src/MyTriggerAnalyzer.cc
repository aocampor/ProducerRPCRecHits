         //2 -*- C++ -*-
//
// Package:    MyAnalyser/MyTriggerAnalyzer
// Class:      MyTriggerAnalyzer
// 
/**\class MyTriggerAnalyzer MyTriggerAnalyzer.cc MyAnalyser/MyTriggerAnalyzer/plugins/MyTriggerAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Piet and Ahmed Sayed
//         Created:  Thu, 20 Nov 2014 14:39:09 GMT
//
// class declaration
// system include files
#include <memory>
#include <fstream>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

// root include files
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// RPC Geometry
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include "DataFormats/Provenance/interface/Timestamp.h"

#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>
#include <Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h>

// #include "DQMServices/Core/interface/DQMStore.h"
// #include "DQMServices/Core/interface/MonitorElement.h"

#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>


// L1 Trigger
#include <DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h>
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include <DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h>
//
// class declaration
//
//

class MyTriggerAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MyTriggerAnalyzer(const edm::ParameterSet&);
      ~MyTriggerAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

    //  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
        
  std::vector<int> m_GMTcandidatesBx;
  std::vector<int> m_DTcandidatesBx;
  std::vector<int> m_RPCcandidatesBx;
  std::vector<int> m_CSCcandidatesBx;

  bool debug;
  std::string rootFileName;
  TFile * outputfile;
 
  
  TH1F * RPC_E_Triggers_ETA_All, * RPC_E_Triggers_ETA_Q0, * RPC_E_Triggers_ETA_Q1, * RPC_E_Triggers_ETA_Q2, * RPC_E_Triggers_ETA_Q3;
  TH1F * RPC_E_Triggers_PHI_All, * RPC_E_Triggers_PHI_Q0, * RPC_E_Triggers_PHI_Q1, * RPC_E_Triggers_PHI_Q2, * RPC_E_Triggers_PHI_Q3;
  TH1F * RPC_B_Triggers_ETA_All, * RPC_B_Triggers_ETA_Q0, * RPC_B_Triggers_ETA_Q1, * RPC_B_Triggers_ETA_Q2, * RPC_B_Triggers_ETA_Q3;
  TH1F * RPC_B_Triggers_PHI_All, * RPC_B_Triggers_PHI_Q0, * RPC_B_Triggers_PHI_Q1, * RPC_B_Triggers_PHI_Q2, * RPC_B_Triggers_PHI_Q3; 
  TH1F * CSC_Triggers_ETA_All, * CSC_Triggers_ETA_Q0, * CSC_Triggers_ETA_Q1, * CSC_Triggers_ETA_Q2, * CSC_Triggers_ETA_Q3;
  TH1F * CSC_Triggers_PHI_All, * CSC_Triggers_PHI_Q0, * CSC_Triggers_PHI_Q1, * CSC_Triggers_PHI_Q2, * CSC_Triggers_PHI_Q3;
  TH1F * DT_Triggers_ETA_All, * DT_Triggers_ETA_Q0, * DT_Triggers_ETA_Q1, * DT_Triggers_ETA_Q2, * DT_Triggers_ETA_Q3;
  TH1F * DT_Triggers_PHI_All, * DT_Triggers_PHI_Q0, * DT_Triggers_PHI_Q1, * DT_Triggers_PHI_Q2, * DT_Triggers_PHI_Q3;

  TH2F * RPC_E_Triggers_ETA_PHI_All;
  TH2F * RPC_B_Triggers_ETA_PHI_All;
  TH2F * CSC_Triggers_ETA_PHI_All;
  TH2F * DT_Triggers_ETA_PHI_All;

  // LumiBlock
  TH1F * MyLumiHistogram_RPCb;
  TH1F * MyLumiHistogram_RPCf;
  TH1F * MyLumiHistogram_DT;
  TH1F * MyLumiHistogram_CSC;

  TH2F * RPC_B_Tr_ETA_lum_All;
  TH2F * RPC_B_Tr_PHI_lum_All;

  TH2F * RPC_E_Tr_ETA_lum_All;
  TH2F * RPC_E_Tr_PHI_lum_All;
  

  TCanvas  * Canvas_RPC_CSC_PHI,  * Canvas_RPC_CSC_ETA, * Canvas_RPC_CSC_ETA_PHI;

  // edm::InputTag m_rpcDigiLabel;
  edm::InputTag m_gtReadoutLabel;
  edm::InputTag m_gmtReadoutLabel;
  
  int countTriggersInLumiSection_RPCb;
  int countTriggersInLumiSection_RPCf;
  int countTriggersInLumiSection_DT;
  int countTriggersInLumiSection_CSC;
  int myCurrentLumiSection;
};

//
// constants, enums and typedefs
//

	int n_phi = 144;
	int n_eta =  48;
	double n_phi_1 =  0.0000; 
	double n_phi_2 =  6.2832;
	double n_eta_1 = -2.40;
	double n_eta_2 =  2.40;
	double n_eta_exact = 39;
	double n_eta_vec[] = {-2.40, -2.30, -2.20, -2.10, -1.97, -1.85, -1.73, -1.61, -1.48, -1.36, -1.24, -1.14, -1.04, -0.93, -0.83, -0.72, -0.58, -0.44, -0.27, -0.07, 0.07, 0.27, 0.44, 0.58, 0.72, 0.83, 0.93, 1.04, 1.14, 1.24, 1.36, 1.48, 1.61, 1.73, 1.85, 1.97, 2.10, 2.20, 2.30, 2.40};


//
// static data member definitions
//

//
// constructors and destructor
//
MyTriggerAnalyzer::MyTriggerAnalyzer(const edm::ParameterSet& iConfig)

{

   if(debug) std::cout<<"MyRPCTriggerAnalyzer :: Constructor]"<<std::endl;

   //now do what ever initialization is needed
   m_gtReadoutLabel     = iConfig.getParameter<edm::InputTag>("GTReadoutRcd");
   m_gmtReadoutLabel    = iConfig.getParameter<edm::InputTag>("GMTReadoutRcd");
   rootFileName         = iConfig.getUntrackedParameter<std::string>("RootFileName");
   debug                = iConfig.getUntrackedParameter<bool>("Debug");

   outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

   // Lumonisity
   MyLumiHistogram_RPCb	= new TH1F("MyLumiHistogram_RPCb",  "MyLumiHistogram_RPCb", 5001 , 1, 5000);
   MyLumiHistogram_RPCf	= new TH1F("MyLumiHistogram_RPCf",  "MyLumiHistogram_RPCf", 5001 , 1, 5000);
   MyLumiHistogram_DT	= new TH1F("MyLumiHistogram_DT",  "MyLumiHistogram_DT", 5001 , 1, 5000);
   MyLumiHistogram_CSC	= new TH1F("MyLumiHistogram_CSC",  "MyLumiHistogram_CSC", 5001 , 1, 5000);

   RPC_B_Tr_ETA_lum_All	= new TH2F("RPC_B_Tr_ETA_lum_All",  "RPC_B_Tr_ETA_lum_All",5001, 1, 5000, n_eta_exact, n_eta_vec);
   RPC_B_Tr_PHI_lum_All	= new TH2F("RPC_B_Tr_PHI_lum_All",  "RPC_B_Tr_PHI_lum_All", 5001, 1, 5000, n_phi, n_phi_1, n_phi_2);
   RPC_E_Tr_ETA_lum_All	= new TH2F("RPC_E_Tr_ETA_lum_All",  "RPC_E_Tr_ETA_lum_All",5001, 1, 5000, n_eta_exact, n_eta_vec);
   RPC_E_Tr_PHI_lum_All	= new TH2F("RPC_E_Tr_PHI_lum_All",  "RPC_E_Tr_PHI_lum_All", 5001, 1, 5000, n_phi, n_phi_1, n_phi_2);

   //RPC_B
   RPC_B_Triggers_PHI_All = new TH1F("RPC_B_Triggers_PHI_All",  "RPC_B_Triggers_PHI_All", n_phi, n_phi_1, n_phi_2);
   RPC_B_Triggers_PHI_Q0  = new TH1F("RPC_B_Triggers_PHI_Q0",   "RPC_B_Triggers_PHI_Q0",  n_phi, n_phi_1, n_phi_2);
   RPC_B_Triggers_PHI_Q1  = new TH1F("RPC_B_Triggers_PHI_Q1",   "RPC_B_Triggers_PHI_Q1",  n_phi, n_phi_1, n_phi_2);
   RPC_B_Triggers_PHI_Q2  = new TH1F("RPC_B_Triggers_PHI_Q2",   "RPC_B_Triggers_PHI_Q2",  n_phi, n_phi_1, n_phi_2);
   RPC_B_Triggers_PHI_Q3  = new TH1F("RPC_B_Triggers_PHI_Q3",   "RPC_B_Triggers_PHI_Q3",  n_phi, n_phi_1, n_phi_2);   
   RPC_B_Triggers_ETA_All = new TH1F("RPC_B_Triggers_ETA_All",  "RPC_B_Triggers_ETA_All", n_eta_exact, n_eta_vec);
   RPC_B_Triggers_ETA_Q0  = new TH1F("RPC_B_Triggers_ETA_Q0",   "RPC_B_Triggers_ETA_Q0",  n_eta_exact, n_eta_vec);
   RPC_B_Triggers_ETA_Q1  = new TH1F("RPC_B_Triggers_ETA_Q1",   "RPC_B_Triggers_ETA_Q1",  n_eta_exact, n_eta_vec);
   RPC_B_Triggers_ETA_Q2  = new TH1F("RPC_B_Triggers_ETA_Q2",   "RPC_B_Triggers_ETA_Q2",  n_eta_exact, n_eta_vec);
   RPC_B_Triggers_ETA_Q3  = new TH1F("RPC_B_Triggers_ETA_Q3",   "RPC_B_Triggers_ETA_Q3",  n_eta_exact, n_eta_vec);
   RPC_B_Triggers_ETA_PHI_All = new TH2F("RPC_B_Triggers_ETA_PHI_All",  "RPC_B_Triggers_ETA_PHI_All", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);

   //RPC_E
   RPC_E_Triggers_PHI_All = new TH1F("RPC_E_Triggers_PHI_All",  "RPC_E_Triggers_PHI_All", n_phi, n_phi_1, n_phi_2);
   RPC_E_Triggers_PHI_Q0  = new TH1F("RPC_E_Triggers_PHI_Q0",   "RPC_E_Triggers_PHI_Q0",  n_phi, n_phi_1, n_phi_2);
   RPC_E_Triggers_PHI_Q1  = new TH1F("RPC_E_Triggers_PHI_Q1",   "RPC_E_Triggers_PHI_Q1",  n_phi, n_phi_1, n_phi_2);
   RPC_E_Triggers_PHI_Q2  = new TH1F("RPC_E_Triggers_PHI_Q2",   "RPC_E_Triggers_PHI_Q2",  n_phi, n_phi_1, n_phi_2);
   RPC_E_Triggers_PHI_Q3  = new TH1F("RPC_E_Triggers_PHI_Q3",   "RPC_E_Triggers_PHI_Q3",  n_phi, n_phi_1, n_phi_2);   
   RPC_E_Triggers_ETA_All = new TH1F("RPC_E_Triggers_ETA_All",  "RPC_E_Triggers_ETA_All", n_eta_exact, n_eta_vec);
   RPC_E_Triggers_ETA_Q0  = new TH1F("RPC_E_Triggers_ETA_Q0",   "RPC_E_Triggers_ETA_Q0",  n_eta_exact, n_eta_vec);
   RPC_E_Triggers_ETA_Q1  = new TH1F("RPC_E_Triggers_ETA_Q1",   "RPC_E_Triggers_ETA_Q1",  n_eta_exact, n_eta_vec);
   RPC_E_Triggers_ETA_Q2  = new TH1F("RPC_E_Triggers_ETA_Q2",   "RPC_E_Triggers_ETA_Q2",  n_eta_exact, n_eta_vec);
   RPC_E_Triggers_ETA_Q3  = new TH1F("RPC_E_Triggers_ETA_Q3",   "RPC_E_Triggers_ETA_Q3",  n_eta_exact, n_eta_vec);
   RPC_E_Triggers_ETA_PHI_All = new TH2F("RPC_E_Triggers_ETA_PHI_All",  "RPC_E_Triggers_ETA_PHI_All", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);

   //CSC
   CSC_Triggers_PHI_All = new TH1F("CSC_Triggers_PHI_All",  "CSC_Triggers_PHI_All", n_phi, n_phi_1, n_phi_2);
   CSC_Triggers_PHI_Q0  = new TH1F("CSC_Triggers_PHI_Q0",   "CSC_Triggers_PHI_Q0",  n_phi, n_phi_1, n_phi_2);
   CSC_Triggers_PHI_Q1  = new TH1F("CSC_Triggers_PHI_Q1",   "CSC_Triggers_PHI_Q1",  n_phi, n_phi_1, n_phi_2);
   CSC_Triggers_PHI_Q2  = new TH1F("CSC_Triggers_PHI_Q2",   "CSC_Triggers_PHI_Q2",  n_phi, n_phi_1, n_phi_2);
   CSC_Triggers_PHI_Q3  = new TH1F("CSC_Triggers_PHI_Q3",   "CSC_Triggers_PHI_Q3",  n_phi, n_phi_1, n_phi_2);   
   CSC_Triggers_ETA_All = new TH1F("CSC_Triggers_ETA_All",  "CSC_Triggers_ETA_All", n_eta_exact, n_eta_vec);
   CSC_Triggers_ETA_Q0  = new TH1F("CSC_Triggers_ETA_Q0",   "CSC_Triggers_ETA_Q0",  n_eta_exact, n_eta_vec);
   CSC_Triggers_ETA_Q1  = new TH1F("CSC_Triggers_ETA_Q1",   "CSC_Triggers_ETA_Q1",  n_eta_exact, n_eta_vec);
   CSC_Triggers_ETA_Q2  = new TH1F("CSC_Triggers_ETA_Q2",   "CSC_Triggers_ETA_Q2",  n_eta_exact, n_eta_vec);
   CSC_Triggers_ETA_Q3  = new TH1F("CSC_Triggers_ETA_Q3",   "CSC_Triggers_ETA_Q3",  n_eta_exact, n_eta_vec);
   CSC_Triggers_ETA_PHI_All = new TH2F("CSC_Triggers_ETA_PHI_All",  "CSC_Triggers_ETA_PHI_All", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);

   //DT
   DT_Triggers_PHI_All = new TH1F("DT_Triggers_PHI_All",  "DT_Triggers_PHI_All", n_phi, n_phi_1, n_phi_2);
   DT_Triggers_PHI_Q0  = new TH1F("DT_Triggers_PHI_Q0",   "DT_Triggers_PHI_Q0",  n_phi, n_phi_1, n_phi_2);
   DT_Triggers_PHI_Q1  = new TH1F("DT_Triggers_PHI_Q1",   "DT_Triggers_PHI_Q1",  n_phi, n_phi_1, n_phi_2);
   DT_Triggers_PHI_Q2  = new TH1F("DT_Triggers_PHI_Q2",   "DT_Triggers_PHI_Q2",  n_phi, n_phi_1, n_phi_2);
   DT_Triggers_PHI_Q3  = new TH1F("DT_Triggers_PHI_Q3",   "DT_Triggers_PHI_Q3",  n_phi, n_phi_1, n_phi_2);   
   DT_Triggers_ETA_All = new TH1F("DT_Triggers_ETA_All",  "DT_Triggers_ETA_All", n_eta_exact, n_eta_vec);
   DT_Triggers_ETA_Q0  = new TH1F("DT_Triggers_ETA_Q0",   "DT_Triggers_ETA_Q0",  n_eta_exact, n_eta_vec);
   DT_Triggers_ETA_Q1  = new TH1F("DT_Triggers_ETA_Q1",   "DT_Triggers_ETA_Q1",  n_eta_exact, n_eta_vec);
   DT_Triggers_ETA_Q2  = new TH1F("DT_Triggers_ETA_Q2",   "DT_Triggers_ETA_Q2",  n_eta_exact, n_eta_vec);
   DT_Triggers_ETA_Q3  = new TH1F("DT_Triggers_ETA_Q3",   "DT_Triggers_ETA_Q3",  n_eta_exact, n_eta_vec);
   DT_Triggers_ETA_PHI_All = new TH2F("DT_Triggers_ETA_PHI_All",  "DT_Triggers_ETA_PHI_All", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);
  //////////

   // Create Canvas
	Canvas_RPC_CSC_PHI = new TCanvas("Canvas_RPC_CSC_PHI","Canvas_RPC_CSC_PHI",1);
	Canvas_RPC_CSC_ETA = new TCanvas("Canvas_RPC_CSC_ETA","Canvas_RPC_CSC_ETA",1);
	Canvas_RPC_CSC_ETA_PHI = new TCanvas("Canvas_RPC_CSC_ETA_PHI","Canvas_RPC_CSC_ETA_PHI",1);

	Canvas_RPC_CSC_PHI->cd(); Canvas_RPC_CSC_PHI->Divide (2,1);
	Canvas_RPC_CSC_PHI->cd(1); RPC_E_Triggers_PHI_All->Draw();
	Canvas_RPC_CSC_PHI->cd(2); CSC_Triggers_PHI_All->Draw();

	Canvas_RPC_CSC_ETA->cd(); Canvas_RPC_CSC_ETA->Divide (2,1);
	Canvas_RPC_CSC_ETA->cd(1); RPC_E_Triggers_ETA_All->Draw(); 
	Canvas_RPC_CSC_ETA->cd(2); CSC_Triggers_ETA_All->Draw();
	
	Canvas_RPC_CSC_ETA_PHI->cd(); Canvas_RPC_CSC_ETA_PHI->Divide (2,1);
	Canvas_RPC_CSC_ETA_PHI->cd(1); 	RPC_E_Triggers_ETA_PHI_All->Draw(); 
	Canvas_RPC_CSC_ETA_PHI->cd(2); 	CSC_Triggers_ETA_PHI_All->Draw();


	//labels
	RPC_B_Triggers_PHI_All->GetXaxis()->SetTitle("#phi"); RPC_B_Triggers_PHI_All->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_B_Triggers_PHI_All->SetTitle("RPC_B_Triggers_PHI_All");
	RPC_B_Triggers_PHI_Q0->GetXaxis()->SetTitle("#phi"); RPC_B_Triggers_PHI_Q0->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_B_Triggers_PHI_Q0->SetTitle("RPCb_Triggers_PHI_Q0");
	RPC_B_Triggers_PHI_Q1->GetXaxis()->SetTitle("#phi"); RPC_B_Triggers_PHI_Q1->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_B_Triggers_PHI_Q1->SetTitle("RPCb_Triggers_PHI_Q1");
	RPC_B_Triggers_PHI_Q2->GetXaxis()->SetTitle("#phi"); RPC_B_Triggers_PHI_Q2->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_B_Triggers_PHI_Q2->SetTitle("RPCb_Triggers_PHI_Q2");
	RPC_B_Triggers_PHI_Q3->GetXaxis()->SetTitle("#phi"); RPC_B_Triggers_PHI_Q3->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_B_Triggers_PHI_Q3->SetTitle("RPCb_Triggers_PHI_Q3");
	RPC_B_Triggers_ETA_All->GetXaxis()->SetTitle("#eta"); RPC_B_Triggers_ETA_All->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_B_Triggers_ETA_All->SetTitle("RPC_B_Triggers_ETA_All");
	RPC_B_Triggers_ETA_Q0->GetXaxis()->SetTitle("#eta"); RPC_B_Triggers_ETA_Q0->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_B_Triggers_ETA_Q0->SetTitle("RPCb_Triggers_ETA_Q0");
	RPC_B_Triggers_ETA_Q1->GetXaxis()->SetTitle("#eta"); RPC_B_Triggers_ETA_Q1->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_B_Triggers_ETA_Q1->SetTitle("RPCb_Triggers_ETA_Q1");
	RPC_B_Triggers_ETA_Q2->GetXaxis()->SetTitle("#eta"); RPC_B_Triggers_ETA_Q2->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_B_Triggers_ETA_Q2->SetTitle("RPCb_Triggers_ETA_Q2");
	RPC_B_Triggers_ETA_Q3->GetXaxis()->SetTitle("#eta"); RPC_B_Triggers_ETA_Q3->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_B_Triggers_ETA_Q3->SetTitle("RPCb_Triggers_ETA_Q3");
	RPC_B_Triggers_ETA_PHI_All->GetXaxis()->SetTitle("#eta"); RPC_B_Triggers_ETA_PHI_All->GetYaxis()->SetTitle("#phi"); 	RPC_B_Triggers_ETA_PHI_All->SetTitle("RPCb Triggers ETA PHI All #Run 240850 Express Stream Data");
	
	RPC_E_Triggers_PHI_All->GetXaxis()->SetTitle("#phi"); RPC_E_Triggers_PHI_All->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_E_Triggers_PHI_All->SetTitle("RPCf_Triggers_PHI_All");
	RPC_E_Triggers_PHI_Q0->GetXaxis()->SetTitle("#phi"); RPC_E_Triggers_PHI_Q0->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_E_Triggers_PHI_Q0->SetTitle("RPCf_Triggers_PHI_Q0");
	RPC_E_Triggers_PHI_Q1->GetXaxis()->SetTitle("#phi"); RPC_E_Triggers_PHI_Q1->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_E_Triggers_PHI_Q1->SetTitle("RPCf_Triggers_PHI_Q1");
	RPC_E_Triggers_PHI_Q2->GetXaxis()->SetTitle("#phi"); RPC_E_Triggers_PHI_Q2->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_E_Triggers_PHI_Q2->SetTitle("RPCf_Triggers_PHI_Q2");
	RPC_E_Triggers_PHI_Q3->GetXaxis()->SetTitle("#phi"); RPC_E_Triggers_PHI_Q3->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_E_Triggers_PHI_Q3->SetTitle("RPCf_Triggers_PHI_Q3");
	RPC_E_Triggers_ETA_All->GetXaxis()->SetTitle("#eta"); RPC_E_Triggers_ETA_All->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_E_Triggers_ETA_All->SetTitle("RPCf_Triggers_ETA_All");
	RPC_E_Triggers_ETA_Q0->GetXaxis()->SetTitle("#eta"); RPC_E_Triggers_ETA_Q0->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_E_Triggers_ETA_Q0->SetTitle("RPCf_Triggers_ETA_Q0");
	RPC_E_Triggers_ETA_Q1->GetXaxis()->SetTitle("#eta"); RPC_E_Triggers_ETA_Q1->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_E_Triggers_ETA_Q1->SetTitle("RPCf_Triggers_ETA_Q1");
	RPC_E_Triggers_ETA_Q2->GetXaxis()->SetTitle("#eta"); RPC_E_Triggers_ETA_Q2->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_E_Triggers_ETA_Q2->SetTitle("RPCf_Triggers_ETA_Q2");
	RPC_E_Triggers_ETA_Q3->GetXaxis()->SetTitle("#eta"); RPC_E_Triggers_ETA_Q3->GetYaxis()->SetTitle("Triggers Rate"); 	RPC_E_Triggers_ETA_Q3->SetTitle("RPCf_Triggers_ETA_Q3");
	RPC_E_Triggers_ETA_PHI_All->GetXaxis()->SetTitle("#eta"); RPC_E_Triggers_ETA_PHI_All->GetYaxis()->SetTitle("#phi"); 	RPC_E_Triggers_ETA_PHI_All->SetTitle("RPCf Triggers ETA PHI All #Run 240850 Express Stream Data");
	
	CSC_Triggers_PHI_All->GetXaxis()->SetTitle("#phi"); CSC_Triggers_PHI_All->GetYaxis()->SetTitle("Triggers Rate"); 	CSC_Triggers_PHI_All->SetTitle("CSC_Triggers_PHI_All");
	CSC_Triggers_PHI_Q0->GetXaxis()->SetTitle("#phi"); CSC_Triggers_PHI_Q0->GetYaxis()->SetTitle("Triggers Rate"); 	CSC_Triggers_PHI_Q0->SetTitle("CSC_Triggers_PHI_Q0");
	CSC_Triggers_PHI_Q1->GetXaxis()->SetTitle("#phi"); CSC_Triggers_PHI_Q1->GetYaxis()->SetTitle("Triggers Rate"); 	CSC_Triggers_PHI_Q1->SetTitle("CSC_Triggers_PHI_Q1");
	CSC_Triggers_PHI_Q2->GetXaxis()->SetTitle("#phi"); CSC_Triggers_PHI_Q2->GetYaxis()->SetTitle("Triggers Rate"); 	CSC_Triggers_PHI_Q2->SetTitle("CSC_Triggers_PHI_Q2");
	CSC_Triggers_PHI_Q3->GetXaxis()->SetTitle("#phi"); CSC_Triggers_PHI_Q3->GetYaxis()->SetTitle("Triggers Rate"); 	CSC_Triggers_PHI_Q3->SetTitle("CSC_Triggers_PHI_Q3");
	CSC_Triggers_ETA_All->GetXaxis()->SetTitle("#eta"); CSC_Triggers_ETA_All->GetYaxis()->SetTitle("Triggers Rate"); 	CSC_Triggers_ETA_All->SetTitle("CSC_Triggers_ETA_All");
	CSC_Triggers_ETA_Q0->GetXaxis()->SetTitle("#eta"); CSC_Triggers_ETA_Q0->GetYaxis()->SetTitle("Triggers Rate"); 	CSC_Triggers_ETA_Q0->SetTitle("CSC_Triggers_ETA_Q0");
	CSC_Triggers_ETA_Q1->GetXaxis()->SetTitle("#eta"); CSC_Triggers_ETA_Q1->GetYaxis()->SetTitle("Triggers Rate"); 	CSC_Triggers_ETA_Q1->SetTitle("CSC_Triggers_ETA_Q1");
	CSC_Triggers_ETA_Q2->GetXaxis()->SetTitle("#eta"); CSC_Triggers_ETA_Q2->GetYaxis()->SetTitle("Triggers Rate"); 	CSC_Triggers_ETA_Q2->SetTitle("CSC_Triggers_ETA_Q2");
	CSC_Triggers_ETA_Q3->GetXaxis()->SetTitle("#eta"); CSC_Triggers_ETA_Q3->GetYaxis()->SetTitle("Triggers Rate"); 	CSC_Triggers_ETA_Q3->SetTitle("CSC_Triggers_ETA_Q3");
	CSC_Triggers_ETA_PHI_All->GetXaxis()->SetTitle("#eta"); CSC_Triggers_ETA_PHI_All->GetYaxis()->SetTitle("#phi"); 	CSC_Triggers_ETA_PHI_All->SetTitle("CSC Triggers ETA PHI All #Run 240850 Express Stream Data");

	DT_Triggers_PHI_All->GetXaxis()->SetTitle("#phi"); DT_Triggers_PHI_All->GetYaxis()->SetTitle("Triggers Rate"); 	DT_Triggers_PHI_All->SetTitle("DT_Triggers_PHI_All");
	DT_Triggers_PHI_Q0->GetXaxis()->SetTitle("#phi"); DT_Triggers_PHI_Q0->GetYaxis()->SetTitle("Triggers Rate"); 	DT_Triggers_PHI_Q0->SetTitle("DT_Triggers_PHI_Q0");
	DT_Triggers_PHI_Q1->GetXaxis()->SetTitle("#phi"); DT_Triggers_PHI_Q1->GetYaxis()->SetTitle("Triggers Rate"); 	DT_Triggers_PHI_Q1->SetTitle("DT_Triggers_PHI_Q1");
	DT_Triggers_PHI_Q2->GetXaxis()->SetTitle("#phi"); DT_Triggers_PHI_Q2->GetYaxis()->SetTitle("Triggers Rate"); 	DT_Triggers_PHI_Q2->SetTitle("DT_Triggers_PHI_Q2");
	DT_Triggers_PHI_Q3->GetXaxis()->SetTitle("#phi"); DT_Triggers_PHI_Q3->GetYaxis()->SetTitle("Triggers Rate"); 	DT_Triggers_PHI_Q3->SetTitle("DT_Triggers_PHI_Q3");
	DT_Triggers_ETA_All->GetXaxis()->SetTitle("#eta"); DT_Triggers_ETA_All->GetYaxis()->SetTitle("Triggers Rate"); 	DT_Triggers_ETA_All->SetTitle("DT_Triggers_ETA_All");
	DT_Triggers_ETA_Q0->GetXaxis()->SetTitle("#eta"); DT_Triggers_ETA_Q0->GetYaxis()->SetTitle("Triggers Rate"); 	DT_Triggers_ETA_Q0->SetTitle("DT_Triggers_ETA_Q0");
	DT_Triggers_ETA_Q1->GetXaxis()->SetTitle("#eta"); DT_Triggers_ETA_Q1->GetYaxis()->SetTitle("Triggers Rate"); 	DT_Triggers_ETA_Q1->SetTitle("DT_Triggers_ETA_Q1");
	DT_Triggers_ETA_Q2->GetXaxis()->SetTitle("#eta"); DT_Triggers_ETA_Q2->GetYaxis()->SetTitle("Triggers Rate"); 	DT_Triggers_ETA_Q2->SetTitle("DT_Triggers_ETA_Q2");
	DT_Triggers_ETA_Q3->GetXaxis()->SetTitle("#eta"); DT_Triggers_ETA_Q3->GetYaxis()->SetTitle("Triggers Rate"); 	DT_Triggers_ETA_Q3->SetTitle("DT_Triggers_ETA_Q3");
	DT_Triggers_ETA_PHI_All->GetXaxis()->SetTitle("#eta"); DT_Triggers_ETA_PHI_All->GetYaxis()->SetTitle("#phi"); 	DT_Triggers_ETA_PHI_All->SetTitle("DT Triggers ETA PHI All #Run 240850 Express Stream Data");

	MyLumiHistogram_RPCb->GetXaxis()->SetTitle("Luminosity Section"); MyLumiHistogram_RPCb->GetYaxis()->SetTitle("RPCb Triggers Rate [Hz]"); MyLumiHistogram_RPCb->SetTitle("Investigate the stability of the RPCb Trigger #Run 240850 Express Stream Data");
	MyLumiHistogram_RPCf->GetXaxis()->SetTitle("Luminosity Section"); MyLumiHistogram_RPCf->GetYaxis()->SetTitle("RPCf Triggers Rate [Hz]"); MyLumiHistogram_RPCf->SetTitle("Investigate the stability of the RPCf Trigger #Run 240850 Express Stream Data");
	MyLumiHistogram_DT->GetXaxis()->SetTitle("Luminosity Section"); MyLumiHistogram_DT->GetYaxis()->SetTitle("DT Triggers Rate [Hz]"); MyLumiHistogram_DT->SetTitle("Investigate the stability of the DT Trigger #Run 240850 Express Stream Data");
	MyLumiHistogram_CSC->GetXaxis()->SetTitle("Luminosity Section"); MyLumiHistogram_CSC->GetYaxis()->SetTitle("CSC Triggers Rate [Hz]"); MyLumiHistogram_CSC->SetTitle("Investigate the stability of the DT Trigger #Run 240850 Express Stream Data");

	RPC_B_Tr_ETA_lum_All->GetXaxis()->SetTitle("Luminosity Section"); RPC_B_Tr_ETA_lum_All->GetYaxis()->SetTitle("#eta"); 	RPC_B_Tr_ETA_lum_All->SetTitle("ETA vs luminosity section RPCb #Run 240850 Express Stream Data");
	RPC_B_Tr_PHI_lum_All->GetXaxis()->SetTitle("Luminosity Section"); RPC_B_Tr_PHI_lum_All->GetYaxis()->SetTitle("#phi"); 	RPC_B_Tr_PHI_lum_All->SetTitle("PHI vs luminosity section RPCb #Run 240850 Express Stream Data");
	RPC_E_Tr_ETA_lum_All->GetXaxis()->SetTitle("Luminosity Section"); RPC_E_Tr_ETA_lum_All->GetYaxis()->SetTitle("#eta");	RPC_E_Tr_ETA_lum_All->SetTitle("ETA vs luminosity section RPCf #Run 240850 Express Stream Data");
	RPC_E_Tr_PHI_lum_All->GetXaxis()->SetTitle("Luminosity Section"); RPC_E_Tr_PHI_lum_All->GetYaxis()->SetTitle("#phi"); 	RPC_E_Tr_PHI_lum_All->SetTitle("PHI vs luminosity section RPCf #Run 240850 Express Stream Data");
	
}


MyTriggerAnalyzer::~MyTriggerAnalyzer()
{
 

   if(debug) std::cout<<"MyTriggerAnalyzer :: Destructor :: begin]"<<std::endl; 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

   outputfile->cd();   
 /*

   
	CSC_Triggers_PHI_All->Write();
	CSC_Triggers_PHI_Q0->Write();
	CSC_Triggers_PHI_Q1->Write();
	CSC_Triggers_PHI_Q2->Write();
	CSC_Triggers_PHI_Q3->Write();  
	CSC_Triggers_ETA_All->Write();
	CSC_Triggers_ETA_Q0->Write();
	CSC_Triggers_ETA_Q1->Write();
	CSC_Triggers_ETA_Q2->Write();
	CSC_Triggers_ETA_Q3->Write();   
	

	DT_Triggers_PHI_All->Write();
	DT_Triggers_PHI_Q0->Write();
	DT_Triggers_PHI_Q1->Write();
	DT_Triggers_PHI_Q2->Write();
	DT_Triggers_PHI_Q3->Write();  
	DT_Triggers_ETA_All->Write();
	DT_Triggers_ETA_Q0->Write();
	DT_Triggers_ETA_Q1->Write();
	DT_Triggers_ETA_Q2->Write();
	DT_Triggers_ETA_Q3->Write();   
	
	*/
   CSC_Triggers_ETA_PHI_All->Write();
   DT_Triggers_ETA_PHI_All->Write();

	Canvas_RPC_CSC_PHI->Write();
	Canvas_RPC_CSC_ETA->Write();
	Canvas_RPC_CSC_ETA_PHI->Write();


	RPC_E_Triggers_PHI_All->Write();
	RPC_E_Triggers_PHI_Q0->Write();
	RPC_E_Triggers_PHI_Q1->Write();
	RPC_E_Triggers_PHI_Q2->Write();
	RPC_E_Triggers_PHI_Q3->Write(); 
	RPC_E_Triggers_ETA_All->Write();
	RPC_E_Triggers_ETA_Q0->Write();
	RPC_E_Triggers_ETA_Q1->Write();
	RPC_E_Triggers_ETA_Q2->Write();
	RPC_E_Triggers_ETA_Q3->Write();   

	RPC_B_Triggers_PHI_All->Write();
	RPC_B_Triggers_PHI_Q0->Write();
	RPC_B_Triggers_PHI_Q1->Write();
	RPC_B_Triggers_PHI_Q2->Write();
	RPC_B_Triggers_PHI_Q3->Write();  
	RPC_B_Triggers_ETA_All->Write();
	RPC_B_Triggers_ETA_Q0->Write();
	RPC_B_Triggers_ETA_Q1->Write();
	RPC_B_Triggers_ETA_Q2->Write();
	RPC_B_Triggers_ETA_Q3->Write();

	RPC_E_Triggers_ETA_PHI_All->Write();
	RPC_B_Triggers_ETA_PHI_All->Write();

	MyLumiHistogram_RPCb->Write();
	MyLumiHistogram_RPCf->Write();
	MyLumiHistogram_DT->Write();
	MyLumiHistogram_CSC->Write();

	RPC_B_Tr_ETA_lum_All->Write();
	RPC_B_Tr_PHI_lum_All->Write();
	RPC_E_Tr_ETA_lum_All->Write();
	RPC_E_Tr_PHI_lum_All->Write();


   outputfile->Close();
   if(debug) std::cout<<"MyTriggerAnalyzer :: Destructor :: end]"<<std::endl; 
   
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyTriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	/*
using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   */
	int evNum = (iEvent.id()).event();
	int rnNum = (iEvent.id()).run();	
	int luNum = (iEvent.id()).luminosityBlock(); myCurrentLumiSection = luNum; 
	//double t = 23.31;

	/*
	y-axis is not "Quality", it is "RPC Trigger Candidates". Or better ...
	Can you divide the entries in each bin by 23,31s * total number of Lumisections of each run (so different numbers for run 229684 and run 229713!). 
	Then you can write on the y-axis "Trigger Rate". Please do this. Then make your plots a bit larger in the presentation.
	*/
	//double Run229684 = 1598;  // Lumisections
	//double Run229713 = 1045;  // Lumisections
	//double Const229684 = 37260;  // seconds

	//double Const232956 = 47063;  // seconds
	//double Const233121 = 41473;  // seconds
	//double Const234304 = 37742;  // seconds
	//double Const234407 = 12279;  // seconds

	double Constxxxyyy = 18696 ; 
	//#Run 240850 Express Stream Data


	//float T_unix = (iEvent.time()).unixTime();

	edm::Handle<L1MuGMTReadoutCollection> pCollection;
	iEvent.getByLabel(m_gmtReadoutLabel,pCollection);

	if ( ! pCollection.isValid() ) {
    //edm::LogError("discriminateGMT") << "can't find L1MuGMTReadoutCollection with label "<< m_gmtReadoutLabel ;
    std::cout<<"can't find L1MuGMTReadoutCollection with label "<< m_gmtReadoutLabel<<std::endl;
    //return -1; 
	}

	  // get GMT readout collection
	const L1MuGMTReadoutCollection * gmtRC = pCollection.product();
  
     // get record vector
	std::vector<L1MuGMTReadoutRecord>::const_iterator RRItr;
	std::vector<L1MuGMTReadoutRecord> gmt_records = gmtRC->getRecords();

	for( RRItr = gmt_records.begin(); RRItr != gmt_records.end(); ++RRItr ) {
   
	// int BxCounter = RRItr->getBxCounter();
	//int BxInEvent = RRItr->getBxInEvent();
    int BxInEventNew = RRItr->getBxNr();
	//std::cout<<" BxCounter= "<<BxCounter<<" BxInEvent= "<<BxInEvent<<" BxInEventNew= "<<BxInEventNew<<std::endl;
	
	
	int nrpcB = 0;
	int nrpcE = 0;
    int ndtB  = 0;
	int nCSC = 0;

	std::vector<L1MuRegionalCand> BrlRpcCands = RRItr->getBrlRPCCands();
	std::vector<L1MuRegionalCand> FwdRPCCands = RRItr->getFwdRPCCands();
	std::vector<L1MuRegionalCand> BrlDtCands  = RRItr->getDTBXCands ();
	std::vector<L1MuRegionalCand> CSCCands = RRItr->getCSCCands();

	std::vector<L1MuRegionalCand>::const_iterator RCItr;
	
		//RPC_B Triggers  
		for( RCItr = BrlRpcCands.begin(); RCItr !=BrlRpcCands.end(); ++RCItr) {
			if ( !(*RCItr).empty() ) {
			m_GMTcandidatesBx.push_back( BxInEventNew );
			nrpcB++;
				if(debug) {
						 std::cout<<"Run :: "<<rnNum<<" Event :: "<<evNum<<" | "; 
						 //std::cout<<"RPC Barrel Trigger :: q = "<<RCItr->quality()<<" pt = "<<RCItr->ptValue()<<" eta = "<<RCItr->etaValue()<<" phi = "<<RCItr->phiValue()<<std::endl;
				} 
				if(!debug) {

					int quality_RPC_B = RCItr->quality();
					double eta_RPC_B = RCItr->etaValue();	
					double phi_RPC_B = RCItr->phiValue();	

					if(debug) std::cout<<"Fill All Histos"<<std::endl;


					
//if ( T_unix>= 1427179140 && T_unix<= 1427180352)
//{

					
					// Triggers Rate
					RPC_B_Triggers_ETA_PHI_All->Fill(eta_RPC_B, phi_RPC_B,1/Constxxxyyy);
					RPC_B_Triggers_ETA_All->Fill(eta_RPC_B,1/Constxxxyyy);
					RPC_B_Triggers_PHI_All->Fill(phi_RPC_B,1/Constxxxyyy);

					// Luminosity Section
					RPC_B_Tr_ETA_lum_All->Fill(myCurrentLumiSection, eta_RPC_B);
					RPC_B_Tr_PHI_lum_All->Fill(myCurrentLumiSection, phi_RPC_B);

						switch (quality_RPC_B) {
							 case 0: {
								if(debug) std::cout<<"Fill Quality 0 Histos"<<std::endl;
								RPC_B_Triggers_ETA_Q0->Fill(eta_RPC_B,1/Constxxxyyy);
								RPC_B_Triggers_PHI_Q0->Fill(phi_RPC_B,1/Constxxxyyy);
								} break;
							case 1: {
								if(debug) std::cout<<"Fill Quality 1 Histos"<<std::endl;
								RPC_B_Triggers_ETA_Q1->Fill(eta_RPC_B,1/Constxxxyyy);
								RPC_B_Triggers_PHI_Q1->Fill(phi_RPC_B,1/Constxxxyyy);

								} break;
							case 2: {
								if(debug) std::cout<<"Fill Quality 2 Histos"<<std::endl;
								RPC_B_Triggers_ETA_Q2->Fill(eta_RPC_B,1/Constxxxyyy);
								RPC_B_Triggers_PHI_Q2->Fill(phi_RPC_B,1/Constxxxyyy);

								} break;
							case 3: {
								if(debug) std::cout<<"Fill Quality 3 Histos"<<std::endl;
								RPC_B_Triggers_ETA_Q3->Fill(eta_RPC_B,1/Constxxxyyy);
								RPC_B_Triggers_PHI_Q3->Fill(phi_RPC_B,1/Constxxxyyy);

								} break;
							default : std::cout<<"Quality_RPC_B = "<<quality_RPC_B<<std::endl;
						}
					++countTriggersInLumiSection_RPCb;
//}

				}
						
			}

		}		
		
		
		//RPC_E Triggers  
		for( RCItr = FwdRPCCands.begin(); RCItr !=FwdRPCCands.end(); ++RCItr) {
			if ( !(*RCItr).empty() ) {
			m_GMTcandidatesBx.push_back( BxInEventNew );
			nrpcE++;
				if(debug) {
						 std::cout<<"Run :: "<<rnNum<<" Event :: "<<evNum<<" | "; 
						 //std::cout<<"RPC Barrel Trigger :: q = "<<RCItr->quality()<<" pt = "<<RCItr->ptValue()<<" eta = "<<RCItr->etaValue()<<" phi = "<<RCItr->phiValue()<<std::endl;
				} 
				if(!debug) {
					
					int quality_RPC_E = RCItr->quality();
					double eta_RPC_E = RCItr->etaValue();	
					double phi_RPC_E = RCItr->phiValue();	

					if(debug) std::cout<<"Fill All Histos"<<std::endl;
//if ( T_unix>= 1427179140 && T_unix<= 1427180352)
//{
					
					// Triggers Rate
					RPC_E_Triggers_ETA_PHI_All->Fill(eta_RPC_E, phi_RPC_E,1/Constxxxyyy);
					RPC_E_Triggers_ETA_All->Fill(eta_RPC_E,1/Constxxxyyy);
					RPC_E_Triggers_PHI_All->Fill(phi_RPC_E,1/Constxxxyyy);

					// Luminosity Section
					RPC_E_Tr_ETA_lum_All->Fill(myCurrentLumiSection, eta_RPC_E);
					RPC_E_Tr_PHI_lum_All->Fill(myCurrentLumiSection, phi_RPC_E);

						switch (quality_RPC_E) {
							 case 0: {
								if(debug) std::cout<<"Fill Quality 0 Histos"<<std::endl;
								RPC_E_Triggers_ETA_Q0->Fill(eta_RPC_E,1/Constxxxyyy);
								RPC_E_Triggers_PHI_Q0->Fill(phi_RPC_E,1/Constxxxyyy);
								} break;
							case 1: {
								if(debug) std::cout<<"Fill Quality 1 Histos"<<std::endl;
								RPC_E_Triggers_ETA_Q1->Fill(eta_RPC_E,1/Constxxxyyy);
								RPC_E_Triggers_PHI_Q1->Fill(phi_RPC_E,1/Constxxxyyy);
								} break;
							case 2: {
								if(debug) std::cout<<"Fill Quality 2 Histos"<<std::endl;
								RPC_E_Triggers_ETA_Q2->Fill(eta_RPC_E,1/Constxxxyyy);
								RPC_E_Triggers_PHI_Q2->Fill(phi_RPC_E,1/Constxxxyyy);
								} break;
							case 3: {
								if(debug) std::cout<<"Fill Quality 3 Histos"<<std::endl;
								RPC_E_Triggers_ETA_Q3->Fill(eta_RPC_E,1/Constxxxyyy);
								RPC_E_Triggers_PHI_Q3->Fill(phi_RPC_E,1/Constxxxyyy);

								} break;
							default : std::cout<<"Quality_RPC_E = "<<quality_RPC_E<<std::endl;
						}
					++countTriggersInLumiSection_RPCf;
//}

				}
						
			}

	
		}
		


		
		// DT Trigger
		for( RCItr = BrlDtCands.begin(); RCItr !=BrlDtCands.end(); ++RCItr) {
			if ( !(*RCItr).empty() ) {
		m_DTcandidatesBx.push_back(BxInEventNew);
		ndtB++;

				if(debug) {
				std::cout<<"Run :: "<<rnNum<<" Event :: "<<evNum<<" | "; 
				//std::cout<<"DT Trigger :: q = "<<RCItr->quality()<<" pt = "<<RCItr->ptValue()<<" eta = "<<RCItr->etaValue()<<" phi = "<<RCItr->phiValue()<<std::endl;
				}
				if(!debug) {

				int quality_DT = RCItr->quality();
				double eta_DT = RCItr->etaValue();
				double phi_DT = RCItr->phiValue();

				if(debug) std::cout<<"Fill All Histos"<<std::endl;
//if ( T_unix>= 1427179140 && T_unix<= 1427180352)
//{
				DT_Triggers_ETA_PHI_All->Fill(eta_DT, phi_DT,1/Constxxxyyy);
				DT_Triggers_ETA_All->Fill(eta_DT,1/Constxxxyyy);
				DT_Triggers_PHI_All->Fill(phi_DT,1/Constxxxyyy);
	  
					switch (quality_DT) {
						 case 0: {
							if(debug) std::cout<<"Fill Quality 0 Histos"<<std::endl;
							DT_Triggers_ETA_Q0->Fill(eta_DT,1/Constxxxyyy);
							DT_Triggers_PHI_Q0->Fill(phi_DT,1/Constxxxyyy);
							} break;
						case 1: {
							if(debug) std::cout<<"Fill Quality 1 Histos"<<std::endl;
							DT_Triggers_ETA_Q1->Fill(eta_DT,1/Constxxxyyy);
							DT_Triggers_PHI_Q1->Fill(phi_DT,1/Constxxxyyy);
							} break;
						case 2: {
							if(debug) std::cout<<"Fill Quality 2 Histos"<<std::endl;
							DT_Triggers_ETA_Q2->Fill(eta_DT,1/Constxxxyyy);
							DT_Triggers_PHI_Q2->Fill(phi_DT,1/Constxxxyyy);
							} break;
						case 3: {
							if(debug) std::cout<<"Fill Quality 3 Histos"<<std::endl;
							DT_Triggers_ETA_Q3->Fill(eta_DT,1/Constxxxyyy);
							DT_Triggers_PHI_Q3->Fill(phi_DT,1/Constxxxyyy);
							} break;
						default : std::cout<<"Quality_DT = "<<quality_DT<<std::endl;
					}
					++countTriggersInLumiSection_DT;
//}
				}
			}
		}

		// CSC Trigger
		for( RCItr = CSCCands.begin(); RCItr !=CSCCands.end(); ++RCItr) {
			if ( !(*RCItr).empty() ) {
		m_DTcandidatesBx.push_back(BxInEventNew);
		nCSC++;
				if(debug) {
				std::cout<<"Run :: "<<rnNum<<" Event :: "<<evNum<<" | "; 
				//std::cout<<"CSC Trigger :: q = "<<RCItr->quality()<<" pt = "<<RCItr->ptValue()<<" eta = "<<RCItr->etaValue()<<" phi = "<<RCItr->phiValue()<<std::endl;
				}
				if(!debug) {

				int quality_CSC = RCItr->quality();
				double eta_CSC = RCItr->etaValue();
				double phi_CSC = RCItr->phiValue();

				if(debug) std::cout<<"Fill All Histos"<<std::endl;
//if ( T_unix>= 1427179140 && T_unix<= 1427180352)
//{
				CSC_Triggers_ETA_PHI_All->Fill(eta_CSC, phi_CSC,1/Constxxxyyy);
				CSC_Triggers_ETA_All->Fill(eta_CSC,1/Constxxxyyy);
				CSC_Triggers_PHI_All->Fill(phi_CSC,1/Constxxxyyy);
	  
					switch (quality_CSC) {
						 case 0: {
							if(debug) std::cout<<"Fill Quality 0 Histos"<<std::endl;
							CSC_Triggers_ETA_Q0->Fill(eta_CSC,1/Constxxxyyy);
							CSC_Triggers_PHI_Q0->Fill(phi_CSC,1/Constxxxyyy);
							} break;
						case 1: {
							if(debug) std::cout<<"Fill Quality 1 Histos"<<std::endl;
							CSC_Triggers_ETA_Q1->Fill(eta_CSC,1/Constxxxyyy);
							CSC_Triggers_PHI_Q1->Fill(phi_CSC,1/Constxxxyyy);
							} break;
						case 2: {
							if(debug) std::cout<<"Fill Quality 2 Histos"<<std::endl;
							CSC_Triggers_ETA_Q2->Fill(eta_CSC,1/Constxxxyyy);
							CSC_Triggers_PHI_Q2->Fill(phi_CSC,1/Constxxxyyy);
							} break;
						case 3: {
							if(debug) std::cout<<"Fill Quality 3 Histos"<<std::endl;
							CSC_Triggers_ETA_Q3->Fill(eta_CSC,1/Constxxxyyy);
							CSC_Triggers_PHI_Q3->Fill(phi_CSC,1/Constxxxyyy);
							} break;
						default : std::cout<<"Quality_CSC = "<<quality_CSC<<std::endl;
					}
					++countTriggersInLumiSection_CSC;
//}
				}

			}
		}



	
	}
}




// ------------ method called once each job just before starting event loop  ------------
void 
MyTriggerAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyTriggerAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MyTriggerAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MyTriggerAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------


void 
MyTriggerAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
	countTriggersInLumiSection_RPCb = 0;
	countTriggersInLumiSection_RPCf = 0;
	countTriggersInLumiSection_DT = 0;
	countTriggersInLumiSection_CSC = 0;
}


// ------------ method called when ending the processing of a luminosity block  ------------
void 
MyTriggerAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
	MyLumiHistogram_RPCb->SetBinContent(myCurrentLumiSection,countTriggersInLumiSection_RPCb/23.31);
	MyLumiHistogram_RPCf->SetBinContent(myCurrentLumiSection,countTriggersInLumiSection_RPCf/23.31);
	MyLumiHistogram_DT->SetBinContent(myCurrentLumiSection,countTriggersInLumiSection_DT/23.31);
	MyLumiHistogram_CSC->SetBinContent(myCurrentLumiSection,countTriggersInLumiSection_CSC/23.31);
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyTriggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyTriggerAnalyzer);
