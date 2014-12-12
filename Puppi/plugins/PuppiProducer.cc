// Adapted from PUPPI code https://github.com/violatingcp/Dummy
// Modified by L CORPE to use the flashgg legacy vertex selector rather than 0th vetrex in collection as Primary Vertex.
// November 2014. 

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
//Main File
#include "fastjet/PseudoJet.hh"
#include "Dummy/Puppi/plugins/PuppiProducer.h"
using namespace edm;
typedef edm::View<reco::Candidate> CandidateView;

// ------------------------------------------------------------------------------------------
PuppiProducer::PuppiProducer(const edm::ParameterSet& iConfig) :
	vertexCandidateMapToken_(consumes<flashgg::VertexCandidateMap>(iConfig.getParameter<edm::InputTag> ("VertexCandidateMapTag"))),
	diPhotonToken_(consumes<View<flashgg::DiPhotonCandidate> >(iConfig.getUntrackedParameter<edm::InputTag> ("DiPhotonTag", InputTag("flashggDiPhotons")))),
	vertexToken_(consumes<View<reco::Vertex> >(iConfig.getUntrackedParameter<InputTag> ("VertexTag", InputTag("offlineSlimmedPrimaryVertices")))),
	pfcandidateToken_(consumes<View<pat::PackedCandidate> >(iConfig.getUntrackedParameter<InputTag> ("PFCandidatesTag", InputTag("packedPFCandidates")))),
	useFlashggVertex(iConfig.getUntrackedParameter<bool>("UseFlashggVertex",false))
{
	fPuppiName = iConfig.getUntrackedParameter<std::string>("PuppiName");
	fUseDZ     = iConfig.getUntrackedParameter<bool>("UseDeltaZCut");
	fDZCut     = iConfig.getUntrackedParameter<double>("DeltaZCut");
	fPuppiContainer = new PuppiContainer(iConfig);
	fPFName    = iConfig.getUntrackedParameter<std::string>("candName"  ,"particleFlow");
	fPVName    = iConfig.getUntrackedParameter<std::string>("vertexName","offlinePrimaryVertices");
	produces<edm::ValueMap<float> > ("PuppiWeights");
	produces<reco::PFCandidateCollection>(fPuppiName);
}
// ------------------------------------------------------------------------------------------
PuppiProducer::~PuppiProducer(){
}
// ------------------------------------------------------------------------------------------
void PuppiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	// Get PFCandidate Collection
	/*	edm::Handle<CandidateView> hPFProduct;
			iEvent.getByLabel(fPFName,hPFProduct);
			assert(hPFProduct.isValid());
			const reco::CandidateView *PFCol = hPFProduct.product();*/

	// Get vertex collection w/PV as the first entry?
	/*	edm::Handle<reco::VertexCollection> hVertexProduct;
			iEvent.getByLabel(fPVName,hVertexProduct);
			assert(hVertexProduct.isValid());
			const reco::VertexCollection *pvCol = hVertexProduct.product(); //no need for PV coll anymore.
			if (pvCol->size());*/
	Handle<View<pat::PackedCandidate> > pfCandidates;
	iEvent.getByToken(pfcandidateToken_,pfCandidates);
	const PtrVector<pat::PackedCandidate>& PFCol = pfCandidates->ptrVector();
	//Fill the reco objects
	fRecoObjCollection.clear();
	Handle<View<reco::Vertex> > primaryVertices;
	iEvent.getByToken(vertexToken_,primaryVertices);
	const PtrVector<reco::Vertex>& pvPtrs = primaryVertices->ptrVector();

	Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
	iEvent.getByToken(diPhotonToken_,diPhotons);
	const PtrVector<flashgg::DiPhotonCandidate>& diPhotonPointers = diPhotons->ptrVector();

	edm::Ptr<reco::Vertex> lpv;
	int counter0=0;
	int counter1=0;
	int counter2=0;
	int counter3=0;
	int counter4=0;
	int counterFindVtx=0;

	if(useFlashggVertex){
		if(diPhotonPointers.size() >0){
			lpv = diPhotonPointers[0]->getVertex(); // The idea is to eventually loop over ALL diphoton canidates in future commit. FIXME
		} else {
			lpv = pvPtrs[0]; 
		} 
	} else {
		lpv = pvPtrs[0]; 
	}

	//std::cout << "DEBUG diphoton[0] position " << lpv->z() << ", 0th vtx pos " << pvPtrs[0]->z() << std::endl;	
	//for(CandidateView::const_iterator itPF = PFCol->begin(); itPF!=PFCol->end(); itPF++) {
	for(unsigned int pfLoop = 0; pfLoop < PFCol.size() ; pfLoop ++) {
		Ptr<pat::PackedCandidate> itPF = (PFCol[pfLoop]); 
		//	pat::PackedCandidate *itPF;
		//	itPF = *(PFCol[pfLoop]);
		RecoObj pReco;
		pReco.pt  = itPF->pt();
		pReco.eta = itPF->eta();
		pReco.phi = itPF->phi();
		pReco.m   = itPF->mass();
		pReco.charge = itPF->charge(); 
		const reco::Vertex *closestVtx = 0;
		double pDZ    = -9999; 
		double pD0    = -9999; 
		int    pVtxId = -9999;
		//bool lFirst = true;
		const pat::PackedCandidate *lPack = dynamic_cast<const pat::PackedCandidate* >(&(*itPF));
		const edm::Ptr<pat::PackedCandidate > lPackEdmPtr = PFCol[pfLoop];


		Handle<flashgg::VertexCandidateMap> vtxmap; //access the flashgg  vtx->PackedCandidate Map 
		iEvent.getByToken(vertexCandidateMapToken_,vtxmap);

		// Now loop over the vertices, and for each, we will check if the considered Packed Candidate is associated. 
		for (std::map<edm::Ptr<reco::Vertex>,edm::PtrVector<pat::PackedCandidate> >::const_iterator vi = vtxmap->begin() ; vi != vtxmap->end() ; vi++) {

			const edm::Ptr<reco::Vertex>  currentVertex = (vi->first); 

			// the arugment of the if returns 1 if lPack is in the vector of PackedCandidates corresponding to currentVertex in the Map.
			if (std::count((vtxmap->at(currentVertex).begin()),vtxmap->at(currentVertex).end(),itPF)) 
			{	
				// by definition, this is therefore the closest vertex.
				closestVtx = &(*currentVertex);
				//Now check if the currentVertex is the same as the legacy PV.
				//Trying to do if (lpv == *currentVertex) gave weird compilation errors, so matching positions is next best thing.
				if (std::fabs(lpv->position().Z() - currentVertex->position().Z())< 0.001 )
				{
					pVtxId =0; //confusingly, pVtxId == 0 when the PackedCandidate's associated vertex matches the Primary Vertex.
					// This is a hangover from the code, which was originally designed such that the 0th vertex in the vertex collection was the PV
					// In that case it would have made sense for pVtxId to be 0 when the vetrex was indeed the PV (ie 0th).
					// Now we are left with slightly weird notation... to be cleared up someday, but works for now.
					counterFindVtx++;
					break;
				}
				pVtxId =1; // conversely, pVtxId >0 when the particle comes from a vertex other than the PV.
			}
		}
		pDZ = lPack->dz(lpv->position()); //Unlike the original code, these quantities need to be calculated from the new PV, not 0th vtx.
		//pDZ = lPack->dz(); //Unlike the original code, these quantities need to be calculated from the new PV, not 0th vtx.
		pD0 = lPack->dxy(lpv->position());
		pReco.dZ = pDZ;
		pReco.d0 = pD0;
		
		if(closestVtx == 0) pReco.vtxId = -1;// means no matchign vertex from map to lpv is found.
		if(closestVtx != 0) pReco.vtxId = pVtxId;
		//if(closestVtx != 0) pReco.vtxChi2 = closestVtx->trackWeight(itPF->trackRef());
		//Set the id for Puppi Algo: 0 is neutral pfCandidate, id = 1 for particles coming from PV and id = 2 for charged particles from non-leading vertex
		pReco.id = 0;
		if(closestVtx != 0 && pVtxId == 0 && fabs(pReco.charge) > 0) {pReco.id = 1; counter1++;}// std::cout<<"DEBUG closestvtx !=0, pvtxid =0, chargedi, set pRedo.id=1"<<std::endl;}
		if(closestVtx != 0 && pVtxId > 0 && fabs(pReco.charge) > 0) {pReco.id = 2;  counter2++;}// std::cout<<"DEBUG closestvtx !=0, pvtxid >0, chargedi, set pRedo.id=2"<<std::endl;}
		//Add a dZ cut if wanted (this helps)
		if(fUseDZ && pDZ > -9999 && closestVtx == 0 && (fabs(pDZ) < fDZCut) && fabs(pReco.charge) > 0) {pReco.id = 1; counter3++ ;}//std::cout<<"DEBUG closestvtx =0, fdz cut, set pRedo.id=1"<<std::endl; }
		if(fUseDZ && pDZ > -9999 && closestVtx == 0 && (fabs(pDZ) > fDZCut) && fabs(pReco.charge) > 0) {pReco.id = 2; counter4++; }//std::cout<<"DEBUG closestvtx =0,  fdz cut set pRedo.id=2"<<std::endl; }
		//std::cout << "pVtxId = " << pVtxId << ", and charge = " << itPF->charge() << ", and closestVtx = " << closestVtx << ", and id = " << pReco.id << std::endl;
		if (pReco.id==0) counter0++;
		fRecoObjCollection.push_back(pReco);
		}

fPuppiContainer->initialize(fRecoObjCollection);
//Compute the weights and the candidates
const std::vector<double> lWeights = fPuppiContainer->puppiWeights();
//Fill it into the event
std::auto_ptr<edm::ValueMap<float> > lPupOut(new edm::ValueMap<float>());
edm::ValueMap<float>::Filler lPupFiller(*lPupOut);
//lPupFiller.insert(hPFProduct,lWeights.begin(),lWeights.end());
lPupFiller.insert(pfCandidates,lWeights.begin(),lWeights.end());
lPupFiller.fill();
iEvent.put(lPupOut,"PuppiWeights");

//Fill a new PF Candidate Collection
const std::vector<fastjet::PseudoJet> lCandidates = fPuppiContainer->puppiParticles();

fPuppiCandidates.reset( new vector<reco::PFCandidate> );    
for(unsigned int i0 = 0; i0 < lCandidates.size(); i0++) {
	reco::PFCandidate pCand(PFCol[lCandidates[i0].user_index()]->charge(),
			PFCol[lCandidates[i0].user_index()]->p4(),
			translatePdgIdToType(PFCol[lCandidates[i0].user_index()]->pdgId()));
	LorentzVector pVec; //pVec.SetPtEtaPhiM(lCandidates[i0].pt(),lCandidates[i0].eta(),lCandidates[i0].phi(),lCandidates[i0].Mass());
	pVec.SetPxPyPzE(lCandidates[i0].px(),lCandidates[i0].py(),lCandidates[i0].pz(),lCandidates[i0].E());
	pCand.setP4(pVec);
	fPuppiCandidates->push_back(pCand);
}

iEvent.put(fPuppiCandidates,fPuppiName);

}
// ------------------------------------------------------------------------------------------
reco::PFCandidate::ParticleType PuppiProducer::translatePdgIdToType(int pdgid) const {
	switch (std::abs(pdgid)) {
		case 211: return reco::PFCandidate::h;
		case 11:  return reco::PFCandidate::e;
		case 13:  return reco::PFCandidate::mu;
		case 22:  return reco::PFCandidate::gamma;
		case 130: return reco::PFCandidate::h0;
		case 1:   return reco::PFCandidate::h_HF;
		case 2:   return reco::PFCandidate::egamma_HF;
		case 0:   return reco::PFCandidate::X;  
		default: return reco::PFCandidate::X;
	}
}
// ------------------------------------------------------------------------------------------
void PuppiProducer::beginJob() {
}
// ------------------------------------------------------------------------------------------
void PuppiProducer::endJob() {
}
// ------------------------------------------------------------------------------------------
void PuppiProducer::beginRun(edm::Run&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void PuppiProducer::endRun(edm::Run&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void PuppiProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void PuppiProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void PuppiProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PuppiProducer);
