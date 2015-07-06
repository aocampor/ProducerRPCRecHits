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


// system include files
#include <memory>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

// root include files
#include <TRandom.h>
#include "TROOT.h"
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
#include "TMultiGraph.h"
#include "TMath.h"                // added

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
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
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include <DataFormats/GeometrySurface/interface/LocalError.h>
#include "DataFormats/Common/interface/Handle.h"

// GenParticle

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// deltaPhi
#include "DataFormats/Math/interface/deltaPhi.h"

//
// class declaration
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
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  
  std::string rootFileName;
  TFile * outputfile;
  
  TH1F * phi_mu, * eta_mu, * pt_mu;
  TH2F * genparticles_ETA_PHI;


};

//
// constants, enums and typedefs
//


int pt_n  = 100;	double pt_x1  = 0,    pt_x2 = 100;

	int n_phi = 144;
	int n_eta =  48;
	double n_phi_1 =  -6.2832; 
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
genAnalyzer::genAnalyzer(const edm::ParameterSet& iConfig)
  
{
  
  //now do what ever initialization is needed
  
  rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");
  
  outputfile = new TFile(rootFileName.c_str(), "RECREATE" );
  
  phi_mu    = new TH1F("phi_mu", "phi_mu", n_phi, n_phi_1, n_phi_2);
  eta_mu    = new TH1F("eta_mu", "eta_mu", n_eta_exact, n_eta_vec);
  pt_mu     = new TH1F("Pt_mu", "Pt_mu", pt_n, pt_x1, pt_x2);
  genparticles_ETA_PHI = new TH2F("genparticles_ETA_PHI",  "genparticles_ETA_PHI", n_eta_exact, n_eta_vec, n_phi, n_phi_1, n_phi_2);
  

 }


genAnalyzer::~genAnalyzer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
  outputfile->cd();
  
  phi_mu ->Write();
  eta_mu ->Write();
  pt_mu ->Write();
  genparticles_ETA_PHI->Write();
  
}

// member functions
//

// ------------ method called for each event  ------------
void
genAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
  

  
  // ===========================================================
  //      new analyzer
  // ===========================================================
  
  std::vector<reco::GenParticle> theGenParticles;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  theGenParticles.insert(theGenParticles.end(),genParticles->begin(),genParticles->end());
  //std::cout << "This Event has " <<  theGenParticles.size() << " genParticles" << std::endl;
  using namespace reco;
  for(unsigned int i=0; i<genParticles->size(); ++i) {
	  const GenParticle & g = (*genParticles)[i];
	  std::cout<<"Gen mu Candidate | id = "<<std::setw(5)<<g.pdgId()<<" | st = "<<std::setw(5)<<g.status()<<" | pt = "<<std::setw(12)<<g.pt();
      std::cout<<" GeV/c | et = "<<std::setw(12)<<g.et()<<" GeV | eta = "<<std::setw(12)<<g.eta()<<" | phi = "<<std::setw(12)<<g.phi()<<std::endl;
	 
	  std::cout<<"!!! Filling Histograms !!!"<<std::endl;

	  phi_mu->Fill(g.phi()); 
	  eta_mu->Fill(g.eta());
	  pt_mu->Fill(g.et());
	  genparticles_ETA_PHI->Fill(g.eta(),g.phi());
	  
  }
}


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
  iSetup.get<MuonGeometryRecord>().get(rpcGeom);
  iSetup.get<MuonGeometryRecord>().get(cscGeom);
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
/*
void 
genAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
genAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowPted parameters for the module  ------------
void
genAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The followPting says we do not know what parameters are allowPted so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(genAnalyzer);
