// Adapted from PUPPI code https://github.com/violatingcp/Dummy
// Modified by L CORPE to use the flashgg legacy vertex selector rather than 0th vetrex in collection as Primary Vertex.
// November 2014. 

#ifndef CommonTools_Puppi_PuppiProducer_h_
#define CommonTools_Puppi_PuppiProducer_h_
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/FwdPtr.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/PtEtaPhiMass.h"
#include "Dummy/Puppi/interface/PuppiContainer.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
using namespace edm;
// ------------------------------------------------------------------------------------------
class PuppiProducer : public edm::EDProducer {

	public:
		explicit PuppiProducer(const edm::ParameterSet&);
		~PuppiProducer();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		typedef math::XYZTLorentzVector LorentzVector;

	private:
		virtual void beginJob() ;
		virtual void produce(edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run&, edm::EventSetup const&);
		virtual void endRun(edm::Run&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
		reco::PFCandidate::ParticleType translatePdgIdToType(int pdgid) const;

		edm::InputTag fPFCands;
		edm::InputTag fVertices;
		std::string     fPuppiName;
		std::string     fPFName;	
		std::string     fPVName;
		edm::EDGetTokenT< flashgg::VertexCandidateMap > vertexCandidateMapToken_;
		edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
		edm::EDGetTokenT<View<reco::Vertex> > vertexToken_;
		edm::EDGetTokenT<View<pat::PackedCandidate> > pfcandidateToken_;
		bool            fUseDZ;
		float           fDZCut;
		bool 						useFlashggVertex;
		PuppiContainer *fPuppiContainer;
		std::vector<RecoObj> fRecoObjCollection;
		//std::auto_ptr <vector<edm::FwdPtr< reco::PFCandidate > > >    fPuppiCandidates;
		std::auto_ptr <vector< reco::PFCandidate > >    fPuppiCandidates;
};
#endif
