// -*- C++ -*-
//
// Package:    Producer/MyRPCProducer
// Class:      MyRPCProducer
// 
/**\class MyRPCProducer MyRPCProducer.cc Producer/MyRPCProducer/plugins/MyRPCProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ahmed Sayed Hamed Ali
//         Created:  Sat, 16 May 2015 07:40:20 GMT
//
//
// system include files
#include <memory>
#include <fstream>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>

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
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
 
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include "DataFormats/Provenance/interface/Timestamp.h"


#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>


#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>
#include <Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h>
#include "DataFormats/Common/interface/Handle.h"


#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "SimMuon/RPCDigitizer/src/RPCSimSetUp.h"
#include "SimMuon/RPCDigitizer/src/RPCDigiProducer.h"
#include "SimMuon/RPCDigitizer/src/RPCDigitizer.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimMuon/RPCDigitizer/src/RPCSynchronizer.h"
 
 //Random Number
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DataFormats/RPCDigi/interface/RPCDigi.h>
#include <DataFormats/MuonData/interface/MuonDigiCollection.h>

//
// class declaration
//

class MyRPCProducer : public edm::EDProducer {
   public:
      explicit MyRPCProducer(const edm::ParameterSet&);
      ~MyRPCProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  //edm::ESHandle <RPCGeometry> rpcGeom;
  edm::InputTag src_ ;
  //edm::InputTag RPCDigiTagSig_ ;
  typedef MuonDigiCollection<RPCDetId, RPCDigi> RPCDigiCollection;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MyRPCProducer::MyRPCProducer(const edm::ParameterSet& iConfig)
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
	//now do what ever initialization is needed
	src_          = iConfig.getParameter<edm::InputTag>("src");
	//RPCDigiTagSig_          = iConfig.getParameter<edm::InputTag>("RPCDigiTagSig");
	produces<RPCDigiCollection>();
}


MyRPCProducer::~MyRPCProducer()
{
 
// do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MyRPCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
// RPCDigiTagSig_ = src_
	
  //////////////////////////////////
  //// RPC DigiCollection producer//
  //////////////////////////////////

	using namespace edm; 
	using namespace std;
	// Loop over RPC digis, copying them from our own local storage
	// retrieve the digis
	Handle<RPCDigiCollection> pRPCdigis; 
	if (iEvent.getByLabel( src_, pRPCdigis )) { 
		auto_ptr<RPCDigiCollection> RPCDigiMerge( new RPCDigiCollection() );
		RPCDigiCollection::DigiRangeIterator RLayerIt;
		for (RLayerIt = pRPCdigis->begin(); RLayerIt != pRPCdigis->end(); ++RLayerIt) {
			// The layerId
			const RPCDetId& layerId = (*RLayerIt).first;
			// Get the iterators over the digis associated with this LayerId
			const RPCDigiCollection::Range& range = (*RLayerIt).second;
			//RPCPointVector.push_back(*RLayerIt);
			RPCDigiMerge->put(range, layerId);
		}
		iEvent.put(RPCDigiMerge);
	}
	
}
// ------------ method called once each job just before starting event loop  ------------
void 
MyRPCProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyRPCProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MyRPCProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MyRPCProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MyRPCProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MyRPCProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyRPCProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyRPCProducer);
