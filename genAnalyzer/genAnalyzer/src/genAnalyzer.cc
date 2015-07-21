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

    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
        
  std::vector<int> m_GMTcandidatesBx;
  std::vector<int> m_DTcandidatesBx;
  std::vector<int> m_RPCcandidatesBx;
  std::vector<int> m_CSCcandidatesBx;

  std::string rootFileName;
  TFile * outputfile;
  // Before
  TH1F * phi_mu_beforeL1_leak, * eta_mu_beforeL1_leak, * pt_mu_beforeL1_leak;
  TH2F * genparticles_ETA_PHI_beforeL1;
  // After
  TH1F * phi_mu_afterL1_leak, * eta_mu_afterL1_leak, * pt_mu_afterL1_leak;
  TH2F * genparticles_ETA_PHI_afterL1_leak;

  // edm::InputTag m_rpcDigiLabel;
  edm::InputTag m_gtReadoutLabel;
  edm::InputTag m_gmtReadoutLabel;
  int countTriggersInLumiSection_RPCb;
  int myCurrentLumiSection;
  int Counter_nrpcB_loop_one;
  int Counter_nrpcB_loop_two;
  int Counter_nrpcB_loop_three;
  int Counter_nrpcB_loop_all;
  
  int n_outside_GMT_loop;
  int n_inside_GMT_outsideRPCb_loop;
  int n_inside_RPCb_loop;

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

    std::cout<<"MyRPCTriggerAnalyzer :: Constructor]"<<std::endl;

	//now do what ever initialization is needed
	m_gtReadoutLabel     = iConfig.getParameter<edm::InputTag>("GTReadoutRcd");
	m_gmtReadoutLabel    = iConfig.getParameter<edm::InputTag>("GMTReadoutRcd");
	rootFileName         = iConfig.getUntrackedParameter<std::string>("RootFileName");

	outputfile = new TFile(rootFileName.c_str(), "RECREATE" );

	phi_mu_beforeL1_leak    = new TH1F("phi_mu_beforeL1_leak", "phi_mu_beforeL1_leak", n_phi, n_phi_1, n_phi_2);
	eta_mu_beforeL1_leak    = new TH1F("eta_mu_beforeL1_leak", "eta_mu_beforeL1_leak", n_eta_exact, n_eta_vec);
	pt_mu_beforeL1_leak     = new TH1F("pt_mu_beforeL1_leak", "pt_mu_beforeL1_leak", pt_n, pt_x1, pt_x2);
	phi_mu_afterL1_leak    = new TH1F("phi_mu_afterL1_leak", "phi_mu_afterL1_leak", n_phi, n_phi_1, n_phi_2);
	eta_mu_afterL1_leak    = new TH1F("eta_mu_afterL1_leak", "eta_mu_afterL1_leak", n_eta_exact, n_eta_vec);
	pt_mu_afterL1_leak     = new TH1F("pt_mu_afterL1_leak", "pt_mu_afterL1_leak", pt_n, pt_x1, pt_x2);
	genparticles_ETA_PHI_beforeL1 = new TH2F("genparticles_ETA_PHI_beforeL1",  "genparticles_ETA_PHI_beforeL1", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);
	genparticles_ETA_PHI_afterL1_leak = new TH2F("genparticles_ETA_PHI_afterL1_leak",  "genparticles_ETA_PHI_afterL1_leak", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);
    
}


genAnalyzer::~genAnalyzer()
{
 
    //std::cout<<"genAnalyzer :: Destructor :: begin]"<<std::endl; 
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    outputfile->cd();
   	phi_mu_beforeL1_leak ->Write();
	phi_mu_afterL1_leak ->Write();
	eta_mu_beforeL1_leak ->Write();
	eta_mu_afterL1_leak ->Write();
	pt_mu_beforeL1_leak ->Write();
	pt_mu_afterL1_leak ->Write();
	genparticles_ETA_PHI_beforeL1->Write();
	genparticles_ETA_PHI_afterL1_leak->Write();
    outputfile->Close();
    //std::cout<<"genAnalyzer :: Destructor :: end]"<<std::endl; 
}


//
// member functions
//

// ------------ method called for each event  ------------
void
genAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	//////////////////////////////////////////////
	// Handling trigger and gen collections///////
	//////////////////////////////////////////////

	// Trigger collection 
	edm::Handle<L1MuGMTReadoutCollection> pCollection;
	iEvent.getByLabel(m_gmtReadoutLabel,pCollection);
	// get GMT readout collection
	const L1MuGMTReadoutCollection * gmtRC = pCollection.product();
	// get record vector
	std::vector<L1MuGMTReadoutRecord>::const_iterator RRItr;
	std::vector<L1MuGMTReadoutRecord> gmt_records = gmtRC->getRecords();
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
		//std::cout<<"Gen mu Candidate | id = "<<std::setw(5)<<g.pdgId()<<" | st = "<<std::setw(5)<<g.status()<<" | pt = "<<std::setw(12)<<g.pt();
		//std::cout<<" GeV/c | et = "<<std::setw(12)<<g.et()<<" GeV | eta = "<<std::setw(12)<<g.eta()<<" | phi = "<<std::setw(12)<<g.phi()<<std::endl;
		//double newgenphi  = normalizedPhi(g.phi());
		//if ( g.eta() < 1.0 && g.eta() > -1.0){
			/* ok std::cout<<"!!! Filling Histograms Befor applying the L1 trigger!!!"<<std::endl;
			 std::cout<<"######The Candidate | id = "<<std::setw(5)<<g.pdgId();
			 std::cout<<"  | eta = "<<std::setw(12)<<g.eta()<<" | phi = "<<std::setw(12)<<g.phi()<<std::endl;*/
			 double genEta = g.eta();
			 double genPhi = g.phi();
			 phi_mu_beforeL1_leak->Fill(g.phi()); 
			 eta_mu_beforeL1_leak->Fill(g.eta());
			 pt_mu_beforeL1_leak->Fill(g.et());
			 genparticles_ETA_PHI_beforeL1->Fill(g.eta(),g.phi());
			 //////////////////////
			 // Applying L1 ///////
			 //////////////////////
			 n_outside_GMT_loop++;
			 std::cout<<"==>==> n outside L1 loop = "<< n_outside_GMT_loop <<std::endl;
			 for( RRItr = gmt_records.begin(); RRItr != gmt_records.end(); ++RRItr ) {
				 int BxInEventNew = RRItr->getBxNr();
				 int nrpcB = 0;
				 std::vector<L1MuRegionalCand> BrlRpcCands = RRItr->getBrlRPCCands();
				 std::vector<L1MuRegionalCand>::const_iterator RCItr;
				 n_inside_GMT_outsideRPCb_loop++;
				 std::cout<<"	==>==>==>==> n inside GMT and outside RPCb loop = "<< n_inside_GMT_outsideRPCb_loop <<std::endl;
				 //RPC_B Triggers  
				 for( RCItr = BrlRpcCands.begin(); RCItr !=BrlRpcCands.end(); ++RCItr) {
					 if ( !(*RCItr).empty() ) {
						 m_GMTcandidatesBx.push_back( BxInEventNew );
						 m_GMTcandidatesBx.size();
						 nrpcB++;
						 n_inside_RPCb_loop++;
						 std::cout<<"		==>==>==>==>==>==> n inside RPCb loop = "<< n_inside_RPCb_loop <<std::endl;

						if (BrlRpcCands.size() < 5 ){
							 //std::cout<<"##### BrlRpcCands size = "<< BrlRpcCands.size()<<std::endl;
							 //std::cout<<"nrpcB = "<< nrpcB <<std::endl;
							 if (nrpcB == 1 ){
								 Counter_nrpcB_loop_one++;
							 }
							 if (nrpcB == 2 ){
								 Counter_nrpcB_loop_two++;
							 }
							 if (nrpcB == 3 ){
								 Counter_nrpcB_loop_three++;
							 }
							 if (nrpcB < 10 ){
								 Counter_nrpcB_loop_all++;
							 }
							 //std::cout<<"Counter_nrpcB_loop_one= "<<Counter_nrpcB_loop_one <<std::endl;
							 //std::cout<<"Counter_nrpcB_loop_two= "<<Counter_nrpcB_loop_two <<std::endl;
							 //std::cout<<"Counter_nrpcB_loop_three= "<<Counter_nrpcB_loop_three <<std::endl;
							 //std::cout<<"Counter_nrpcB_loop_all= "<<Counter_nrpcB_loop_all <<std::endl;
						}


						 //ok std::cout<<"Applying the L1 trigger!!!"<< " nrpcB++ = "<< nrpcB++<<std::endl;
						 double eta_RPC_B = RCItr->etaValue();	
						 double phi_RPC_B = RCItr->phiValue();
						 //calculate Delta R between gen and L1
						 float dR = deltaR(eta_RPC_B, phi_RPC_B, genEta,genPhi);
						 // ok std::cout<<"Delta R = "<< dR <<std::endl;
						 if(dR < 0.5){ //  Delta R condition
							/* ok std::cout<<"!!! Filling Histograms After applying the L1 trigger!!!"<<std::endl;
							std::cout<<"######The Trigger Candidate | id = "<<std::setw(5)<<g.pdgId();
							std::cout<<"  | eta = "<<std::setw(12)<<g.eta()<<" | phi = "<<std::setw(12)<<g.phi()<<std::endl; */
							phi_mu_afterL1_leak->Fill(g.phi()); 
							eta_mu_afterL1_leak->Fill(g.eta());
							pt_mu_afterL1_leak->Fill(g.et());
							genparticles_ETA_PHI_afterL1_leak->Fill(g.eta(),g.phi());
						 } // end of Delta R condition
					 } // end of if stament 	
				 } // end for RPCb trigger candidate loop
			 } // end of GMT loop
		//}// end of ETA cut
	} // end of gen loop
}// end of the method




// ------------ method called once each job just before starting event loop  ------------
void 
genAnalyzer::beginJob()
{
 Counter_nrpcB_loop_one = 0;
 Counter_nrpcB_loop_two = 0;
 Counter_nrpcB_loop_three = 0;
 Counter_nrpcB_loop_all = 0;
 
 n_outside_GMT_loop = 0;
 n_inside_GMT_outsideRPCb_loop = 0;
 n_inside_RPCb_loop = 0;
 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
genAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void 
genAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a run  ------------

void 
genAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{

}


// ------------ method called when starting to processes a luminosity block  ------------


void 
genAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{

}


// ------------ method called when ending the processing of a luminosity block  ------------
void 
genAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{

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
