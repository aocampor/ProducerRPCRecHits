 // -*- C++ -*- Space modification
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
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
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
  edm::ESHandle <RPCGeometry> rpcGeom;
  edm::InputTag src_;
  //  typedef std::vector< const TrackingRecHit * >  Point;
  //  typedef std::vector<Point> RPC_RecHit_Colliection;
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
  src_  = iConfig.getParameter<edm::InputTag>( "src" );
  produces<RPCRecHitCollection>();//.setBranchAlias( "AhRecHits");
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
  ///////////////////////////
  //// RPC recHits producer//
  ///////////////////////////
  
  using namespace edm; 
  using namespace std;
  // retrieve the rpcrechits
  Handle<RPCRecHitCollection> rpcRecHits; 
  //iEvent.getByToken(src_,rpcRecHits);
  iEvent.getByLabel( src_, rpcRecHits );
  //std::cout << "I am in the producer" << std::endl;
  // create the vectors. Use auto_ptr, as these pointers will automatically
  // delete when they go out of scope, a very efficient way to reduce memory leaks.
  auto_ptr<RPCRecHitCollection> recHits( new RPCRecHitCollection() );
  // and already reserve some space for the new data, to control the size
  // of your executible's memory use.
  //const int size = rpcRecHits->size();
  //  recHits->reserve( size );
  //std::cout << "1" << std::endl;
  ESHandle<RPCGeometry> rpcGeo;
  iSetup.get<MuonGeometryRecord>().get(rpcGeo);
  //std::cout << "2" << std::endl;
  const std::vector< const RPCRoll * > & rolls = rpcGeo->rolls();

  //  RPCRecHit rh;
  edm::OwnVector<RPCRecHit> RPCPointVector;
  RPCRecHitCollection::const_iterator recHit;
  // loop to fill the vector 
  //std::cout << "3" << std::endl;
  for (int i = 0; i < (int) rolls.size() ; ++i){
    RPCDetId did =  rolls[i]->id();
    RPCPointVector.clear();
    //std::cout << "4" << std::endl;
    for (recHit = rpcRecHits->begin(); recHit != rpcRecHits->end(); recHit++) {
      //std::cout << "a5" << std::endl;
      RPCDetId rollId = (RPCDetId)(*recHit).rpcId();
      //std::cout << "b5" << std::endl;
      if(rollId != did)
	continue;
      //std::cout << "c5" << std::endl;
      RPCGeomServ rpcsrv(rollId);
      //std::cout << "d5" << std::endl;
      //LocalPoint recHitPos=recHit->localPosition();
      std::cout << "e5 " << rollId.region() << " " << rollId.ring() << " " << rollId.sector() << " " << rollId.subsector() << std::endl;
      std::cout << "e5 " << did.region() << " " << did.ring() << " " << did.sector() << " " << did.subsector() << std::endl;
      
      //const RPCRoll * rollasociated =  rpcGeom->roll( rollId);
      //std::cout << "f5" << std::endl;
      //const BoundPlane & RPCSurface = rollasociated->surface(); 
      //std::cout << "g5" << std::endl;
      //GlobalPoint RPCGlobalPoint = RPCSurface.toGlobal(recHitPos);
      //std::cout << "5" << std::endl;
      //  W_0 is Masked 
      if (rollId.region() == 0 && rollId.ring() != 0) {
	//std::cout<<"RPC Rec Hit in "<<rpcsrv.name();
	//// std::cout<<" Local Position = "<<recHitPos<<" Local Uncrt = "<<
	//std::cout<<" Global Position = "<<RPCGlobalPoint; // <<" Global Uncrt = "<<;
	//std::cout<<" Clustersize = "<<recHit->clusterSize()<<" First strip of cluster = "<<recHit->firstClusterStrip()<<" BX = "<<recHit->BunchX()<<std::endl;
	// fill the RPCRecHits in the vectors
	//rh = (*recHit);
	RPCPointVector.push_back(*recHit);
	
      }
    }
    //std::cout << "6" << std::endl;
    recHits->put(did, RPCPointVector.begin(),RPCPointVector.end()  );
  }

  // and save the vectors
  iEvent.put( recHits );//, "AhmedRecHits" );
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
