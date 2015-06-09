// -*- C++ -*-
//
// Package:    MyProducer/MyRPCProducer
// Class:      MyRPCProducer
// 
/**\class MyRPCProducer MyRPCProducer.cc MyProducer/MyRPCProducer/plugins/MyRPCProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alberto Andres Ocampo Rios
//         Created:  Mon, 08 Jun 2015 10:11:25 GMT
//
//


// system include files
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCRollSpecs.h"
#include "Geometry/RPCGeometry/interface/RPCChamber.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

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
 
  edm::InputTag src_;
     
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
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
  src_ = iConfig.getParameter<edm::InputTag>("src");
  produces<RPCDigiCollection>();
   //now do what ever other initialization is needed
  
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
   using namespace edm;
   using namespace std;

   Handle<RPCDigiCollection> pRPCdigis; 
   iEvent.getByLabel( src_, pRPCdigis );

   ESHandle<RPCGeometry> rpcGeo;
   iSetup.get<MuonGeometryRecord>().get(rpcGeo);

   auto_ptr<RPCDigiCollection> RPCDigiMerge( new RPCDigiCollection() );
   //TrackingGeometry::DetContainer::const_iterator it=rpcGeo->dets().begin();
   for (RPCGeometry::DetIdContainer::const_iterator it=rpcGeo->detIds().begin();
	it != rpcGeo->detIds().end(); it++){
     const RPCChamber* ch =  rpcGeo->chamber( *it );
     if(ch != 0){
       std::vector< const RPCRoll*> roles = (ch->rolls());
       for(std::vector<const RPCRoll*>::const_iterator r = roles.begin();
	   r != roles.end(); ++r){
	 RPCDetId rpcId = (*r)->id();
	
	 // Filter out sector10!
	 if( rpcId.region() == 0 && rpcId.sector() == 10 ){
	   //cout << "Filtering out sector 10!!!!!!!!!!!!!!!!!!! " << endl;
	   continue;
	 }
	
/*
	// Filter out Leaky RPC in station 2!

	if (
	(rpcId.region()!= 0&& rpcId.ring() != -1&& rpcId.station() != 2 && rpcId.sector() != 4 && rpcId.layer()!= 1)
	||(rpcId.region()!= 0&& rpcId.ring() != -1&& rpcId.station() != 2 && rpcId.sector() != 7 && rpcId.layer()!= 2)
	||(rpcId.region()!= 0&& rpcId.ring() != 0&& rpcId.station() != 2 && rpcId.sector() != 4 && rpcId.layer()!= 1)
	||(rpcId.region()!= 0&& rpcId.ring() != 0&& rpcId.station() != 2 && rpcId.sector() != 7 && rpcId.layer()!= 1)
	||(rpcId.region()!= 0&& rpcId.ring() != 2&& rpcId.station() != 2 && rpcId.sector() != 5 && rpcId.layer()!= 2)
	||(rpcId.region()!= 0&& rpcId.ring() != 2&& rpcId.station() != 2 && rpcId.sector() != 8 && rpcId.layer()!= 1)
	){
	   cout << "Filtering out leaky RPC in station 2! " << endl;
	   continue;
	 }

*/
	 RPCDigiCollection::Range rpcRangeDigi= pRPCdigis->get(rpcId);
	 std::vector<long> dig;
	 std::vector< std::vector <long> > prev;
	 prev.clear();
	 for (RPCDigiCollection::const_iterator digiIt = rpcRangeDigi.first;
	      digiIt!=rpcRangeDigi.second;++digiIt){
	   //std::cout << "Bx: " << (*digiIt).bx() << " rawId: " << rpcId.rawId() << " strip: " << (*digiIt).strip() << std::endl;
	   dig.clear();
	   bool verheit = false;
	   dig.push_back((*digiIt).bx());
	   dig.push_back(rpcId.rawId());
	   dig.push_back((*digiIt).strip());
	   for(unsigned int i=0; i<prev.size(); i++){
	     if(prev[i][0]==dig[0] and prev[i][1]==dig[1] and prev[i][2] == dig[2]) verheit = true;
	     //cout << "\t prev: " << prev[i][0] << " " << prev[i][1] << " " << prev[i][2] << endl;
	   }
	   if(!verheit){
	     prev.push_back(dig);
	     (*RPCDigiMerge).insertDigi(rpcId,(*digiIt));
	   }
	   
	 }
       }
     }
   }

   iEvent.put( RPCDigiMerge );
   
 
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
