// -*- C++ -*-
//
// Package:    genAnalyzer/genAnalyzer
// Class:      genAnalyzer
// 
/**\class genAnalyzer genAnalyzer.cc genAnalyzer/genAnalyzer/plugins/genAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ahmed Sayed Hamed Ali
//         Created:  Sun, 05 Jul 2015 14:19:59 GMT
//
//
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
#include "TMath.h"                // added

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
// why is this one not included by default
#include "FWCore/Framework/interface/ESHandle.h"

// Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include <DataFormats/GeometrySurface/interface/LocalError.h>
#include <DataFormats/GeometryVector/interface/LocalPoint.h>
#include <DataFormats/GeometryVector/interface/LocalPoint.h>
#include "DataFormats/Common/interface/Handle.h"

// GenParticle

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// deltaPhi
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CommonTools/Utils/interface/normalizedPhi.h"

// delta R
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//
//


class genAnalyzer : public edm::EDAnalyzer {
   public:
      explicit genAnalyzer(const edm::ParameterSet&);
      ~genAnalyzer();

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
 

  TH1F * phi_mu_beforeL1, * eta_mu_beforeL1, * pt_mu_beforeL1;
  TH2F * genparticles_ETA_PHI_beforeL1;

  TH1F * phi_mu_afterL1, * eta_mu_afterL1, * pt_mu_afterL1;
  TH2F * genparticles_ETA_PHI_afterL1;
  

  TH1F * RPC_B_Triggers_ETA_All, * RPC_B_Triggers_ETA_Q0, * RPC_B_Triggers_ETA_Q1, * RPC_B_Triggers_ETA_Q2, * RPC_B_Triggers_ETA_Q3;
  TH1F * RPC_B_Triggers_PHI_All, * RPC_B_Triggers_PHI_Q0, * RPC_B_Triggers_PHI_Q1, * RPC_B_Triggers_PHI_Q2, * RPC_B_Triggers_PHI_Q3; 

  TH2F * RPC_B_Triggers_ETA_PHI_All;


  // LumiBlock
  TH1F * MyLumiHistogram_RPCb;

  // edm::InputTag m_rpcDigiLabel;
  edm::InputTag m_gtReadoutLabel;
  edm::InputTag m_gmtReadoutLabel;
  
  
  int countTriggersInLumiSection_RPCb;
  int myCurrentLumiSection;
};

//
//  constants, enums and typedefs
//

	int n_phi = 144;
	int n_eta =  48;
	double n_phi_1 =  -6.2832; 
	double n_phi_2 =  6.2832;
	double n_eta_1 = -2.40;
	double n_eta_2 =  2.40;
	double n_eta_exact = 39;
	double n_eta_vec[] = {-2.40, -2.30, -2.20, -2.10, -1.97, -1.85, -1.73, -1.61, -1.48, -1.36, -1.24, -1.14, -1.04, -0.93, -0.83, -0.72, -0.58, -0.44, -0.27, -0.07, 0.07, 0.27, 0.44, 0.58, 0.72, 0.83, 0.93, 1.04, 1.14, 1.24, 1.36, 1.48, 1.61, 1.73, 1.85, 1.97, 2.10, 2.20, 2.30, 2.40};
	int pt_n  = 100;	double pt_x1  = 0,    pt_x2 = 100;
//
// static data member definitions
//

//
// constructors and destructor
//
genAnalyzer::genAnalyzer(const edm::ParameterSet& iConfig)

{

   if(debug) std::cout<<"MyRPCTriggerAnalyzer :: Constructor]"<<std::endl;

   //now do what ever initialization is needed
   m_gtReadoutLabel     = iConfig.getParameter<edm::InputTag>("GTReadoutRcd");
   m_gmtReadoutLabel    = iConfig.getParameter<edm::InputTag>("GMTReadoutRcd");
   rootFileName         = iConfig.getUntrackedParameter<std::string>("RootFileName");
   debug                = iConfig.getUntrackedParameter<bool>("Debug");
   

   outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

	phi_mu_beforeL1    = new TH1F("phi_mu_beforeL1", "phi_mu_beforeL1", n_phi, n_phi_1, n_phi_2);
	eta_mu_beforeL1    = new TH1F("eta_mu_beforeL1", "eta_mu_beforeL1", n_eta_exact, n_eta_vec);
	pt_mu_beforeL1     = new TH1F("pt_mu_beforeL1", "pt_mu_beforeL1", pt_n, pt_x1, pt_x2);
	phi_mu_afterL1    = new TH1F("phi_mu_afterL1", "phi_mu_afterL1", n_phi, n_phi_1, n_phi_2);
	eta_mu_afterL1    = new TH1F("eta_mu_afterL1", "eta_mu_afterL1", n_eta_exact, n_eta_vec);
	pt_mu_afterL1     = new TH1F("pt_mu_afterL1", "pt_mu_afterL1", pt_n, pt_x1, pt_x2);

	genparticles_ETA_PHI_beforeL1 = new TH2F("genparticles_ETA_PHI_beforeL1",  "genparticles_ETA_PHI_beforeL1", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);
	genparticles_ETA_PHI_afterL1 = new TH2F("genparticles_ETA_PHI_afterL1",  "genparticles_ETA_PHI_afterL1", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);

   // Lumonisity
   MyLumiHistogram_RPCb	= new TH1F("MyLumiHistogram_RPCb",  "MyLumiHistogram_RPCb", 5001 , 1, 5000);

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
	

	MyLumiHistogram_RPCb->GetXaxis()->SetTitle("Luminosity Section"); MyLumiHistogram_RPCb->GetYaxis()->SetTitle("RPCb Triggers Rate [Hz]"); MyLumiHistogram_RPCb->SetTitle("Investigate the stability of the RPCb Trigger #Run 240850 Express Stream Data");

	
}


genAnalyzer::~genAnalyzer()
{
 

   if(debug) std::cout<<"genAnalyzer :: Destructor :: begin]"<<std::endl; 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

   outputfile->cd();

//	MyLumiHistogram_RPCb->Write();
	
	
  	/*RPC_B_Triggers_PHI_All->Write();
	RPC_B_Triggers_PHI_Q0->Write();
	RPC_B_Triggers_PHI_Q1->Write();
	RPC_B_Triggers_PHI_Q2->Write();
	RPC_B_Triggers_PHI_Q3->Write();  
	RPC_B_Triggers_ETA_All->Write();
	RPC_B_Triggers_ETA_Q0->Write();
	RPC_B_Triggers_ETA_Q1->Write();
	RPC_B_Triggers_ETA_Q2->Write();
	RPC_B_Triggers_ETA_Q3->Write();*/
	
   	phi_mu_beforeL1 ->Write();
	phi_mu_afterL1 ->Write();
	eta_mu_beforeL1 ->Write();
	eta_mu_afterL1 ->Write();
	pt_mu_beforeL1 ->Write();
	pt_mu_afterL1 ->Write();
	genparticles_ETA_PHI_beforeL1->Write();
	genparticles_ETA_PHI_afterL1->Write();
	RPC_B_Triggers_ETA_PHI_All->Write();
	
   outputfile->Close();
   if(debug) std::cout<<"genAnalyzer :: Destructor :: end]"<<std::endl; 
   
}


//
// member functions
//

// ------------ method called for each event  ------------
void
genAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	//float T_unix = (iEvent.time()).unixTime();
	int evNum = (iEvent.id()).event();
	int rnNum = (iEvent.id()).run();	
	int luNum = (iEvent.id()).luminosityBlock(); myCurrentLumiSection = luNum;
	double Constxxxyyy = 1 ; // for calculation of trigger rate 
	
	
	//////////////////////////////////////////////
	// Handling trigger and gen collections///////
	//////////////////////////////////////////////

	// Trigger collection 
	edm::Handle<L1MuGMTReadoutCollection> pCollection;
	iEvent.getByLabel(m_gmtReadoutLabel,pCollection);
	// GenParticles collection
	std::vector<reco::GenParticle> theGenParticles;
	edm::Handle<reco::GenParticleCollection> genParticles;
	iEvent.getByLabel("genParticles", genParticles);
	theGenParticles.insert(theGenParticles.end(),genParticles->begin(),genParticles->end());
	//std::cout << "This Event has " <<  theGenParticles.size() << " genParticles" << std::endl;
	using namespace reco;
	// loop over genparticles
	for(unsigned int i=0; i<genParticles->size(); ++i) {
		const GenParticle & g = (*genParticles)[i];
		std::cout<<"Gen mu Candidate | id = "<<std::setw(5)<<g.pdgId()<<" | st = "<<std::setw(5)<<g.status()<<" | pt = "<<std::setw(12)<<g.pt();
		std::cout<<" GeV/c | et = "<<std::setw(12)<<g.et()<<" GeV | eta = "<<std::setw(12)<<g.eta()<<" | phi = "<<std::setw(12)<<g.phi()<<std::endl;
		std::cout<<"!!! Filling Histograms befor applying the L1 trigger!!!"<<std::endl;
		//double newgenphi  = normalizedPhi(g.phi());
		 if ( g.eta() < 1 && g.eta() > -1){
			double genEta = g.eta();
			double genPhi = g.phi();
			phi_mu_beforeL1->Fill(g.phi()); 
			eta_mu_beforeL1->Fill(g.eta());
			pt_mu_beforeL1->Fill(g.et());
			genparticles_ETA_PHI_beforeL1->Fill(g.eta(),g.phi());
			
			//////////////////////
			// Applying L1 ///////
			//////////////////////
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
				int BxInEventNew = RRItr->getBxNr();
				// int BxCounter = RRItr->getBxCounter();
				//int BxInEvent = RRItr->getBxInEvent();
				//std::cout<<" BxCounter= "<<BxCounter<<" BxInEvent= "<<BxInEvent<<" BxInEventNew= "<<BxInEventNew<<std::endl;
				int nrpcB = 0;
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
							// Triggers Rate
							double newgenphi  = normalizedPhi(phi_RPC_B);
							RPC_B_Triggers_ETA_PHI_All->Fill(eta_RPC_B, newgenphi,1/Constxxxyyy);
							RPC_B_Triggers_ETA_All->Fill(eta_RPC_B,1/Constxxxyyy);
							RPC_B_Triggers_PHI_All->Fill(phi_RPC_B,1/Constxxxyyy);
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
							} // end of trigger quality
							++countTriggersInLumiSection_RPCb;
							//calculate Delta R between gen and L1
							float dR = deltaR(eta_RPC_B, phi_RPC_B, genEta,genPhi);
							std::cout<<"Delta R = "<< dR <<std::endl;
							if(dR < 0.5){ //  Delta R condition
								std::cout<<"!!! Filling Histograms after applying the L1 trigger!!!"<<std::endl;
								phi_mu_afterL1->Fill(g.phi()); 
								eta_mu_afterL1->Fill(g.eta());
								pt_mu_afterL1->Fill(g.et());
								genparticles_ETA_PHI_afterL1->Fill(g.eta(),g.phi());
							} // end of Delta R condition
						} // end of if statment debug
					} // end of if stament 	
				} // end for RPCb trigger candidate loop
			} // end of GMT loop
		 }// end of ETA cut
	} // end of gen loop
}// end of the method




// ------------ method called once each job just before starting event loop  ------------
void 
genAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
genAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
genAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
genAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------


void 
genAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
	countTriggersInLumiSection_RPCb = 0;
}


// ------------ method called when ending the processing of a luminosity block  ------------
void 
genAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
	MyLumiHistogram_RPCb->SetBinContent(myCurrentLumiSection,countTriggersInLumiSection_RPCb/23.31);
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
genAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(genAnalyzer);
