#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToOne.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "KrAFT/GeneratorTools/interface/Types.h"

#include "DataFormats/Common/interface/Wrapper.h"
#include "KrAFT/GeneratorTools/interface/GenTTbbCandidate.h"
#include <vector>

namespace pat {
  edm::RefProd<std::vector<pat::Jet> > dummy00;

  // pat::Jet -> reco::GenJet mapping
  edm::Wrapper<edm::AssociationMap<edm::OneToOne<std::vector<pat::Jet>, std::vector<reco::GenJet> ,unsigned int> > > dummy10;
  RecoToGenJetMap dummy11;
}

namespace reco {
  edm::RefProd<std::vector<reco::GenJet> > dummy01;

  // reco::GenJet -> reco::GenParticle mapping
  edm::Wrapper<edm::AssociationMap<edm::OneToMany<std::vector<reco::GenJet>,std::vector<reco::GenParticle>,unsigned int> > > dummy12;

}
namespace {
  struct KrAFT_GenericNtuple {

    vallot::GenTTbbCandidate dummyGenTTbbCandidate;
    edm::Wrapper<vallot::GenTTbbCandidate> dummyGenTTbbCandidateWrapper;
    std::vector<vallot::GenTTbbCandidate> dummyGenTTbbCandidateCollection;
    edm::Wrapper<std::vector<vallot::GenTTbbCandidate> > dummyGenTTbbCandidateCollectionWrapper;
    edm::Ptr<vallot::GenTTbbCandidate> dummyGenTTbbCandidatePtr;

  };

}

