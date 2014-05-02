#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "KrAFT/GeneratorTools/interface/GenTTbbCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"
#include "TFile.h"
#include <vector>
#include <TH1F.h>
#include <TH2F.h>

using namespace edm;
using namespace std;

class TTbar2bGenFilter : public edm::EDFilter
{
public:
  TTbar2bGenFilter(const edm::ParameterSet& pset);
  ~TTbar2bGenFilter() {};

  void beginJob();
  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);
  void endJob() {};
  bool isLastQuark(const reco::GenParticle&, const int&);
  bool isFromtop(const reco::GenParticle&);

  string debug;

private:
  bool applyFilter_;
  edm::InputTag genParticlesLabel_;
  edm::InputTag genJetsLabel_;
  int type_;

  //TTree* tree;
  /*
  TH1F* b_status3_daughterid;

  TH1F* b_from_top_pt;
  TH1F* b_from_top_multi;
  TH1F* b_from_top_motherid;
  TH1F* b_from_top_status;

  TH1F* b_from_nontop_pt;
  TH1F* b_from_nontop_multi;
  TH1F* b_from_nontop_motherid;
  TH1F* b_from_nontop_status;

  TH1F* b_multiplicity;
  */

  TH1F* h_multiplicity_bQuarks;
  TH1F* h_multiplicity_bGenJets;

  TH1F* h_multiplicity_bQuarks20;
  TH1F* h_multiplicity_bQuarks20DILVIS;
  TH1F* h_multiplicity_bQuarks20DILVISTTBB;

  TH1F* h_multiplicity_bGenJets20;
  TH1F* h_multiplicity_bGenJets40;

  TH1F* h_multiplicity_bGenJets20DILVIS;
  TH1F* h_multiplicity_bGenJets20DILVISTTBB;

  TH1F* h_multiplicity_addbGenJets20;
  TH1F* h_multiplicity_addbGenJets40;
  TH1F* h_multiplicity_addGenJets20;
  TH1F* h_multiplicity_addGenJets40;

  TH1F* h_multiplicity_GenJets;
  TH1F* h_multiplicity_GenJets10;
  TH1F* h_multiplicity_GenJets15;
  TH1F* h_multiplicity_GenJets20;
  TH1F* h_multiplicity_GenJets25;
  TH1F* h_multiplicity_GenJets30;

  TH1F* h_multiplicity_GenJetsDIL;
  TH1F* h_multiplicity_GenJets10DIL;
  TH1F* h_multiplicity_GenJets15DIL;
  TH1F* h_multiplicity_GenJets20DIL;
  TH1F* h_multiplicity_GenJets25DIL;
  TH1F* h_multiplicity_GenJets30DIL;

  TH1F* h_multiplicity_GenJetsDILVIS;
  TH1F* h_multiplicity_GenJets10DILVIS;
  TH1F* h_multiplicity_GenJets15DILVIS;
  TH1F* h_multiplicity_GenJets20DILVIS;
  TH1F* h_multiplicity_GenJets25DILVIS;
  TH1F* h_multiplicity_GenJets30DILVIS;

  TH1F* h_multiplicity_GenJetsDILVISTTBB;
  TH1F* h_multiplicity_GenJets10DILVISTTBB;
  TH1F* h_multiplicity_GenJets15DILVISTTBB;
  TH1F* h_multiplicity_GenJets20DILVISTTBB;
  TH1F* h_multiplicity_GenJets25DILVISTTBB;
  TH1F* h_multiplicity_GenJets30DILVISTTBB;

  TH1F* h_multiplicity_GenJetsDILVISTTCC;
  TH1F* h_multiplicity_GenJets10DILVISTTCC;
  TH1F* h_multiplicity_GenJets15DILVISTTCC;
  TH1F* h_multiplicity_GenJets20DILVISTTCC;
  TH1F* h_multiplicity_GenJets25DILVISTTCC;
  TH1F* h_multiplicity_GenJets30DILVISTTCC;

  TH1F* h_dR_additional_bquarks;
  TH1F* h_pt_additional_bquarks_1;
  TH1F* h_eta_additional_bquarks_1;
  TH1F* h_pt_additional_bquarks_2;
  TH1F* h_eta_additional_bquarks_2;
  TH2F* h_eta_mass_additional_bquarks_1;
  TH2F* h_eta_mass_additional_bquarks_2;
  TH1F* h_dR_additional_bquarks_20GeV;
  TH1F* h_dR_additional_bquarks_40GeV;

  TH1F* h_nEvents;
  TH1F* h_nEvents_addjets;
  TH1F* h_nEvents_parton;
  TH1F* h_nEvents_addparton;

  TH1F* h_nEvents_40GeV;
  TH1F* h_nEvents_40GeV_addjets;
  TH1F* h_nEvents_40GeV_parton;
  TH1F* h_nEvents_40GeV_addparton;

  TH1F* h_nEvents_inclusive;
  TH1F* h_nEvents_parton_inclusive;

  std::vector<vallot::GenTTbbCandidate>* ttbarGen;

};

TTbar2bGenFilter::TTbar2bGenFilter(const edm::ParameterSet& pset)
{
  applyFilter_= pset.getUntrackedParameter<bool>("applyFilter",false);
  genParticlesLabel_= pset.getParameter<edm::InputTag>("genParticlesLabel");
  genJetsLabel_= pset.getParameter<edm::InputTag>("genJetsLabel");
  type_ = pset.getUntrackedParameter<int>("type",0);

  edm::Service<TFileService> fs;
  //tree = fs->make<TTree>("tree", "Tree for Top quark study");

  /*
  b_from_top_pt  = fs->make<TH1F>( "b_from_top_pt"  , "p_{T}", 100,  0., 150. );
  b_from_top_multi  = fs->make<TH1F>( "b_from_top_multi"  , "Multiplicity", 10,  0, 10 );
  b_from_top_motherid  = fs->make<TH1F>( "b_from_top_motherid"  , "Mother PdgId", 3000,  0, 3000 );
  b_from_top_status  = fs->make<TH1F>( "b_from_top_status"  , "status", 10,  0, 10 );

  b_from_nontop_pt  = fs->make<TH1F>( "b_from_nontop_pt"  , "p_{T}", 100,  0., 150. );
  b_from_nontop_multi  = fs->make<TH1F>( "b_from_nontop_multi"  , "Multiplicity", 10,  0, 10 );
  b_from_nontop_motherid  = fs->make<TH1F>( "b_from_nontop_motherid"  , "Mother PdgId", 3000,  0, 3000 );
  b_from_nontop_status  = fs->make<TH1F>( "b_from_nontop_status"  , "status", 10,  0, 10 );

  b_multiplicity  = fs->make<TH1F>( "b_multiplicity"  , "Multiplicity", 10,  0, 10 );
  */

  h_multiplicity_bQuarks  = fs->make<TH1F>( "h_multiplicity_bQuarks"  , "Multiplicity", 10,  0, 10 );
  h_multiplicity_bGenJets  = fs->make<TH1F>( "h_multiplicity_bGenJets"  , "Multiplicity", 10,  0, 10 );

  h_dR_additional_bquarks = fs->make<TH1F>("h_dR_additional_bquarks", "#Delta R", 500, 0, 5);
  h_pt_additional_bquarks_1 = fs->make<TH1F>("h_pt_additional_bquarks_1", "p_{T}", 200, 0, 200);
  h_eta_additional_bquarks_1 = fs->make<TH1F>("h_eta_additional_bquarks_1", "#Eta", 600, -3, 3);
  h_pt_additional_bquarks_2 = fs->make<TH1F>("h_pt_additional_bquarks_2", "p_{T}", 200, 0, 200);
  h_eta_additional_bquarks_2 = fs->make<TH1F>("h_eta_additional_bquarks_2", "#Eta", 600, -3, 3);
  h_eta_mass_additional_bquarks_1 = fs->make<TH2F>("h_eta_mass_additional_bquarks_1", ";Mass;#eta", 1000, 0, 10, 600, -3, 3);
  h_eta_mass_additional_bquarks_2 = fs->make<TH2F>("h_eta_mass_additional_bquarks_2", ";Mass;#eta", 1000, 0, 10, 600, -3, 3);
  h_dR_additional_bquarks_20GeV = fs->make<TH1F>("h_dR_additional_bquarks_20GeV", "#Delta R", 500, 0, 5);
  h_dR_additional_bquarks_40GeV = fs->make<TH1F>("h_dR_additional_bquarks_40GeV", "#Delta R", 500, 0, 5);

  h_multiplicity_bQuarks20  = fs->make<TH1F>( "h_multiplicity_bQuarks20"  , "Multiplicity", 10,  0, 10 );
  h_multiplicity_bQuarks20DILVIS  = fs->make<TH1F>( "h_multiplicity_bQuarks20DILVIS"  , "Multiplicity", 10,  0, 10 );
  h_multiplicity_bQuarks20DILVISTTBB  = fs->make<TH1F>( "h_multiplicity_bQuarks20DILVISTTBB"  , "Multiplicity", 10,  0, 10 );

  h_multiplicity_bGenJets20  = fs->make<TH1F>( "h_multiplicity_bGenJets20"  , "Multiplicity", 10,  0, 10 );
  h_multiplicity_bGenJets40  = fs->make<TH1F>( "h_multiplicity_bGenJets40"  , "Multiplicity", 10,  0, 10 );
  h_multiplicity_bGenJets20DILVIS  = fs->make<TH1F>( "h_multiplicity_bGenJets20DILVIS"  , "Multiplicity", 10,  0, 10 );
  h_multiplicity_bGenJets20DILVISTTBB  = fs->make<TH1F>( "h_multiplicity_bGenJets20DILVISTTBB"  , "Multiplicity", 10,  0, 10 );

  h_multiplicity_addbGenJets20  = fs->make<TH1F>( "h_multiplicity_addbGenJets20"  , "b-Jet Multiplicity", 10,  0, 10 );
  h_multiplicity_addbGenJets40  = fs->make<TH1F>( "h_multiplicity_addbGenJets40"  , "b-Jet Multiplicity", 10,  0, 10 );
  h_multiplicity_addGenJets20  = fs->make<TH1F>( "h_multiplicity_addGenJets20"  , "Jet Multiplicity", 10,  0, 10 );
  h_multiplicity_addGenJets40  = fs->make<TH1F>( "h_multiplicity_addGenJets40"  , "Jet Multiplicity", 10,  0, 10 );

  h_multiplicity_GenJets  = fs->make<TH1F>( "h_multiplicity_GenJets"  , "Multiplicity", 30,  0, 30 );
  h_multiplicity_GenJets10  = fs->make<TH1F>( "h_multiplicity_GenJets10"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets15  = fs->make<TH1F>( "h_multiplicity_GenJets15"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets20  = fs->make<TH1F>( "h_multiplicity_GenJets20"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets25  = fs->make<TH1F>( "h_multiplicity_GenJets25"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets30  = fs->make<TH1F>( "h_multiplicity_GenJets30"  , "Multiplicity", 12,  0, 12 );
  
  h_multiplicity_GenJetsDIL  = fs->make<TH1F>( "h_multiplicity_GenJetsDIL"  , "Multiplicity", 30,  0, 30 );
  h_multiplicity_GenJets10DIL  = fs->make<TH1F>( "h_multiplicity_GenJets10DIL"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets15DIL  = fs->make<TH1F>( "h_multiplicity_GenJets15DIL"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets20DIL  = fs->make<TH1F>( "h_multiplicity_GenJets20DIL"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets25DIL  = fs->make<TH1F>( "h_multiplicity_GenJets25DIL"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets30DIL  = fs->make<TH1F>( "h_multiplicity_GenJets30DIL"  , "Multiplicity", 12,  0, 12 );

  h_multiplicity_GenJetsDILVIS  = fs->make<TH1F>( "h_multiplicity_GenJetsDILVIS"  , "Multiplicity", 30,  0, 30 );
  h_multiplicity_GenJets10DILVIS  = fs->make<TH1F>( "h_multiplicity_GenJets10DILVIS"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets15DILVIS  = fs->make<TH1F>( "h_multiplicity_GenJets15DILVIS"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets20DILVIS  = fs->make<TH1F>( "h_multiplicity_GenJets20DILVIS"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets25DILVIS  = fs->make<TH1F>( "h_multiplicity_GenJets25DILVIS"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets30DILVIS  = fs->make<TH1F>( "h_multiplicity_GenJets30DILVIS"  , "Multiplicity", 12,  0, 12 );

  h_multiplicity_GenJetsDILVISTTBB  = fs->make<TH1F>( "h_multiplicity_GenJetsDILVISTTBB"  , "Multiplicity", 30,  0, 30 );
  h_multiplicity_GenJets10DILVISTTBB  = fs->make<TH1F>( "h_multiplicity_GenJets10DILVISTTBB"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets15DILVISTTBB  = fs->make<TH1F>( "h_multiplicity_GenJets15DILVISTTBB"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets20DILVISTTBB  = fs->make<TH1F>( "h_multiplicity_GenJets20DILVISTTBB"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets25DILVISTTBB  = fs->make<TH1F>( "h_multiplicity_GenJets25DILVISTTBB"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets30DILVISTTBB  = fs->make<TH1F>( "h_multiplicity_GenJets30DILVISTTBB"  , "Multiplicity", 12,  0, 12 );

  h_multiplicity_GenJetsDILVISTTCC  = fs->make<TH1F>( "h_multiplicity_GenJetsDILVISTTCC"  , "Multiplicity", 30,  0, 30 );
  h_multiplicity_GenJets10DILVISTTCC  = fs->make<TH1F>( "h_multiplicity_GenJets10DILVISTTCC"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets15DILVISTTCC  = fs->make<TH1F>( "h_multiplicity_GenJets15DILVISTTCC"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets20DILVISTTCC  = fs->make<TH1F>( "h_multiplicity_GenJets20DILVISTTCC"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets25DILVISTTCC  = fs->make<TH1F>( "h_multiplicity_GenJets25DILVISTTCC"  , "Multiplicity", 12,  0, 12 );
  h_multiplicity_GenJets30DILVISTTCC  = fs->make<TH1F>( "h_multiplicity_GenJets30DILVISTTCC"  , "Multiplicity", 12,  0, 12 );

  h_nEvents = fs->make<TH1F>( "h_nEvents"  , "h_nEvents", 7,  0, 7 );
  h_nEvents_addjets = fs->make<TH1F>( "h_nEvents_addjets"  , "h_nEvents_addjets", 7,  0, 7 );
  h_nEvents_parton = fs->make<TH1F>( "h_nEvents_parton"  , "h_nEvents_parton", 6,  0, 6 );
  h_nEvents_addparton = fs->make<TH1F>( "h_nEvents_addparton"  , "h_nEvents_addparton", 6,  0, 6 );

  h_nEvents_40GeV = fs->make<TH1F>( "h_nEvents_40GeV"  , "h_nEvents_40GeV", 7,  0, 7 );
  h_nEvents_40GeV_addjets = fs->make<TH1F>( "h_nEvents_40GeV_addjets"  , "h_nEvents_40GeV_addjets", 7,  0, 7 );
  h_nEvents_40GeV_parton = fs->make<TH1F>( "h_nEvents_40GeV_parton"  , "h_nEvents_40GeV_parton", 6,  0, 6 );
  h_nEvents_40GeV_addparton = fs->make<TH1F>( "h_nEvents_40GeV_addparton"  , "h_nEvents_40GeV_addparton", 6,  0, 6 );

  h_nEvents_inclusive = fs->make<TH1F>( "h_nEvents_inclusive"  , "h_nEvents", 6,  0, 6 );
  h_nEvents_parton_inclusive = fs->make<TH1F>( "h_nEvents_parton_inclusive"  , "h_nEvents_parton", 6,  0, 6 );

  ttbarGen = new std::vector<vallot::GenTTbbCandidate>();
}

void TTbar2bGenFilter::beginJob(){
  //tree->Branch("ttbarGen","std::vector<vallot::CMGTTbarCandidate>", &ttbarGen);
}

bool TTbar2bGenFilter::filter(edm::Event& iEvent, const edm::EventSetup& eventSetup)
{
  const bool isRealData = iEvent.isRealData();

  if (isRealData )
    return true;

  ttbarGen->clear();

  bool accepted = false;

  using namespace std;
  using namespace edm;
  using namespace reco;
  debug = "";
  debug += "---------EVENT start---------\n" ;

  edm::Handle<std::vector<reco::GenJet> > genJets_;
  iEvent.getByLabel(genJetsLabel_,genJets_);

  edm::Handle<reco::GenParticleCollection> genParticles_;
  iEvent.getByLabel(genParticlesLabel_,genParticles_);

  const reco::GenParticleCollection* myGenParticles = 0;
  myGenParticles = &(*genParticles_);

  /*

  unsigned int nParticles = myGenParticles->size();

  int nb_from_top = 0;
  int nb_from_nontop = 0;
  int nb = 0;

  for ( unsigned int ip=0; ip<nParticles; ++ip ) { 

    const reco::GenParticle& p = (*myGenParticles)[ip];
   
    if ( abs(p.pdgId()) != 5 ) continue;
    bool isLast = isLastQuark(p, 5);
    if (isLast != true) continue;
 
    int status = p.status();

    const reco::GenParticle* mother = dynamic_cast<const reco::GenParticle*>(p.mother());
    int motherAbsPdgId = abs(mother->pdgId());

    bool isfromtop = isFromtop(p);

    if( !isfromtop ) {  
      b_from_nontop_status->Fill(status);
      b_from_nontop_motherid->Fill(motherAbsPdgId);
      b_from_nontop_pt->Fill(p.pt());
      nb_from_nontop++ ;
    }else{
      b_from_top_status->Fill(status);
      b_from_top_motherid->Fill(motherAbsPdgId);
      b_from_top_pt->Fill(p.pt());
      nb_from_top++  ;
    }

    nb++;
  }

  b_multiplicity->Fill(nb);
  */

  //gen information
  vallot::GenTTbbCandidate ttbarGenLevel;

  if(genJets_.isValid()){
    const std::vector<reco::GenJet>* myGenJets = 0;
    myGenJets = &(*genJets_);

    ttbarGenLevel.building(myGenJets, myGenParticles);

    ttbarGen->push_back(ttbarGenLevel);

    double dR_addb = deltaR(ttbarGenLevel.bquarks3().eta(), ttbarGenLevel.bquarks3().phi(),ttbarGenLevel.bquarks4().eta(), ttbarGenLevel.bquarks4().phi()); 
    if( dR_addb > 0 ){
      if( ttbarGenLevel.bquarks3().mass() > 0 && ttbarGenLevel.bquarks4().mass() > 0) { 
        h_dR_additional_bquarks->Fill(dR_addb);
        h_pt_additional_bquarks_1->Fill( ttbarGenLevel.bquarks3().pt() );
        h_eta_additional_bquarks_1->Fill(ttbarGenLevel.bquarks3().eta() );
        h_pt_additional_bquarks_2->Fill( ttbarGenLevel.bquarks4().pt() );
        h_eta_additional_bquarks_2->Fill(ttbarGenLevel.bquarks4().eta() );
        h_eta_mass_additional_bquarks_1->Fill(ttbarGenLevel.bquarks3().mass(), ttbarGenLevel.bquarks3().eta());
        h_eta_mass_additional_bquarks_2->Fill(ttbarGenLevel.bquarks4().mass(), ttbarGenLevel.bquarks4().eta());
        if( ttbarGenLevel.bquarks3().pt() > 20 && ttbarGenLevel.bquarks3().pt() > 20) h_dR_additional_bquarks_20GeV->Fill(dR_addb); 
        if( ttbarGenLevel.bquarks3().pt() > 40 && ttbarGenLevel.bquarks3().pt() > 40) h_dR_additional_bquarks_40GeV->Fill(dR_addb); 
      }
    }

    bool dil = ttbarGenLevel.diLeptonic() == 1 ;
    bool divis = ttbarGenLevel.lepton1().pt() > 20 && abs(ttbarGenLevel.lepton1().eta()) < 2.4 && ttbarGenLevel.lepton2().pt() > 20 && abs(ttbarGenLevel.lepton2().eta()) < 2.4 ;

    bool vis = ttbarGenLevel.lepton1().pt() > 20 && abs(ttbarGenLevel.lepton1().eta()) < 2.4 && ttbarGenLevel.lepton2().pt() > 20 && abs(ttbarGenLevel.lepton2().eta()) < 2.4 && ttbarGenLevel.NbJets20(type_) >= 2;
    bool nbjets2 = ttbarGenLevel.NbJets20(type_) >= 2;
    bool ttb = ttbarGenLevel.NbJets20(type_) == 3;
    bool njets4 = ttbarGenLevel.NJets20() >= 4;
    bool ttbb = ttbarGenLevel.NbJets20(type_) >= 4;
    bool ttcc = ttbarGenLevel.NcJets20(type_) >= 2;

    bool addvis = ttbarGenLevel.lepton1().pt() > 20 && abs(ttbarGenLevel.lepton1().eta()) < 2.4 && ttbarGenLevel.lepton2().pt() > 20 && abs(ttbarGenLevel.lepton2().eta()) < 2.4 && ttbarGenLevel.NaddJets20() >= 2;
    bool addttb = ttbarGenLevel.NaddbJets20(type_) == 1;
    bool addttbb = ttbarGenLevel.NaddbJets20(type_) >= 2;
 
    bool vis_40GeV = ttbarGenLevel.lepton1().pt() > 20 && abs(ttbarGenLevel.lepton1().eta()) < 2.4 && ttbarGenLevel.lepton2().pt() > 20 && abs(ttbarGenLevel.lepton2().eta()) < 2.4 && ttbarGenLevel.NbJets40(type_) >= 2;
    bool nbjets2_40GeV = ttbarGenLevel.NbJets40(type_) >= 2;
    bool ttb_40GeV = ttbarGenLevel.NbJets40(type_) == 3;
    bool njets4_40GeV = ttbarGenLevel.NJets40() >= 4;
    bool ttbb_40GeV = ttbarGenLevel.NbJets40(type_) >= 4;
    bool ttcc_40GeV = ttbarGenLevel.NcJets40(type_) >= 2;

    bool addvis_40GeV = ttbarGenLevel.lepton1().pt() > 20 && abs(ttbarGenLevel.lepton1().eta()) < 2.4 && ttbarGenLevel.lepton2().pt() > 20 && abs(ttbarGenLevel.lepton2().eta()) < 2.4 && ttbarGenLevel.NaddJets40() >= 2;
    bool addttb_40GeV = ttbarGenLevel.NaddbJets40(type_) == 1;
    bool addttbb_40GeV = ttbarGenLevel.NaddbJets40(type_) >= 2;

    bool nbpartons2 = ttbarGenLevel.NbQuarks20() >= 2;
    bool ttbb_parton = ttbarGenLevel.NbQuarks20() >= 4;

    bool nbpartons2_40GeV = ttbarGenLevel.NbQuarks40() >= 2;
    bool ttbb_parton_40GeV = ttbarGenLevel.NbQuarks40() >= 4;

    bool addttb_parton = ttbarGenLevel.NaddbQuarks20() == 1;
    bool addttbb_parton = ttbarGenLevel.NaddbQuarks20() >= 2;
 
    bool addttb_parton_40GeV = ttbarGenLevel.NaddbQuarks40() == 1;
    bool addttbb_parton_40GeV = ttbarGenLevel.NaddbQuarks40() >= 2;

    h_multiplicity_bQuarks->Fill( ttbarGenLevel.NbQuarks() ) ;
    h_multiplicity_bQuarks20->Fill( ttbarGenLevel.NbQuarks20() ) ;
    if( dil && vis ) h_multiplicity_bQuarks20DILVIS->Fill( ttbarGenLevel.NbQuarks20() );
    if( dil && vis && ttbb) h_multiplicity_bQuarks20DILVISTTBB->Fill( ttbarGenLevel.NbQuarks20() );
    h_multiplicity_bGenJets->Fill( ttbarGenLevel.NbJets(type_) );
    h_multiplicity_bGenJets20->Fill( ttbarGenLevel.NbJets20(type_) );
    h_multiplicity_bGenJets40->Fill( ttbarGenLevel.NbJets40(type_) );
    if( dil && vis ) h_multiplicity_bGenJets20DILVIS->Fill( ttbarGenLevel.NbJets20(type_) );
    if( dil && vis && ttbb ) h_multiplicity_bGenJets20DILVISTTBB->Fill( ttbarGenLevel.NbJets20(type_) );

    h_multiplicity_addbGenJets20->Fill( ttbarGenLevel.NaddbJets20(type_) );
    h_multiplicity_addbGenJets40->Fill( ttbarGenLevel.NaddbJets40(type_) );
    h_multiplicity_addGenJets20->Fill( ttbarGenLevel.NaddJets20() );
    h_multiplicity_addGenJets40->Fill( ttbarGenLevel.NaddJets40() );    

    h_nEvents->Fill(0);
    if( dil ) h_nEvents->Fill(1);
    if( dil && vis ) h_nEvents->Fill(2);
    if( dil && vis && njets4  ) h_nEvents->Fill(3);
    if( dil && vis && njets4 && ttbb  ) h_nEvents->Fill(4);
    if( dil && vis && njets4 && !ttbb && ttb  ) h_nEvents->Fill(5);
    if( dil && vis && njets4 && !ttbb && !ttb && ttcc  ) h_nEvents->Fill(6);

    h_nEvents_addjets->Fill(0);
    if( dil ) h_nEvents_addjets->Fill(1);
    if( dil && divis ) h_nEvents_addjets->Fill(2);
    if( dil && addvis  ) h_nEvents_addjets->Fill(3);
    if( dil && addvis && addttbb ) h_nEvents_addjets->Fill(4);
    if( dil && addvis && !addttbb && addttb ) h_nEvents_addjets->Fill(5);
    if( dil && addvis && !addttbb && !addttb && ttcc ) h_nEvents_addjets->Fill(6);

    h_nEvents_parton->Fill(0);
    if( dil ) h_nEvents_parton->Fill(1);
    if( dil && divis && nbpartons2 ) h_nEvents_parton->Fill(2);
    if( dil && divis && nbpartons2 && njets4  ) h_nEvents_parton->Fill(3);
    if( dil && divis && nbpartons2 && njets4 && ttbb_parton  ) h_nEvents_parton->Fill(4);
    if( dil && divis && nbpartons2 && njets4 && !ttbb_parton && ttcc  ) h_nEvents_parton->Fill(5);

    h_nEvents_addparton->Fill(0);
    if( dil ) h_nEvents_addparton->Fill(1);
    if( dil && divis ) h_nEvents_addparton->Fill(2);
    if( dil && addvis  ) h_nEvents_addparton->Fill(3);
    if( dil && addvis && addttbb_parton ) h_nEvents_addparton->Fill(4);
    if( dil && addvis && addttb_parton ) h_nEvents_addparton->Fill(5);

    h_nEvents_40GeV->Fill(0);
    if( dil ) h_nEvents_40GeV->Fill(1);
    if( dil && vis_40GeV ) h_nEvents_40GeV->Fill(2);
    if( dil && vis_40GeV && njets4_40GeV  ) h_nEvents_40GeV->Fill(3);
    if( dil && vis_40GeV && njets4_40GeV && ttbb_40GeV  ) h_nEvents_40GeV->Fill(4);
    if( dil && vis_40GeV && njets4_40GeV && !ttbb_40GeV && ttb_40GeV  ) h_nEvents_40GeV->Fill(5);
    if( dil && vis_40GeV && njets4_40GeV && !ttbb_40GeV && !ttb_40GeV && ttcc_40GeV  ) h_nEvents_40GeV->Fill(6);

    h_nEvents_40GeV_addjets->Fill(0);
    if( dil ) h_nEvents_40GeV_addjets->Fill(1);
    if( dil && divis ) h_nEvents_40GeV_addjets->Fill(2);
    if( dil && addvis_40GeV  ) h_nEvents_40GeV_addjets->Fill(3);
    if( dil && addvis_40GeV && addttbb_40GeV ) h_nEvents_40GeV_addjets->Fill(4);
    if( dil && addvis_40GeV && !addttbb_40GeV && addttb_40GeV ) h_nEvents_40GeV_addjets->Fill(5);
    if( dil && addvis_40GeV && !addttbb_40GeV && !addttb_40GeV && ttcc_40GeV ) h_nEvents_40GeV_addjets->Fill(6);

    h_nEvents_40GeV_parton->Fill(0);
    if( dil ) h_nEvents_40GeV_parton->Fill(1);
    if( dil && divis && nbpartons2_40GeV ) h_nEvents_40GeV_parton->Fill(2);
    if( dil && divis && nbpartons2_40GeV && njets4  ) h_nEvents_40GeV_parton->Fill(3);
    if( dil && divis && nbpartons2_40GeV && njets4 && ttbb_parton_40GeV  ) h_nEvents_40GeV_parton->Fill(4);
    if( dil && divis && nbpartons2_40GeV && njets4 && !ttbb_parton_40GeV && ttcc  ) h_nEvents_40GeV_parton->Fill(5);    

    h_nEvents_40GeV_addparton->Fill(0);
    if( dil ) h_nEvents_40GeV_addparton->Fill(1);
    if( dil && divis ) h_nEvents_40GeV_addparton->Fill(2);
    if( dil && addvis_40GeV  ) h_nEvents_40GeV_addparton->Fill(3);
    if( dil && addvis_40GeV && addttbb_parton_40GeV ) h_nEvents_40GeV_addparton->Fill(4);
    if( dil && addvis_40GeV && addttb_parton_40GeV ) h_nEvents_40GeV_addparton->Fill(5);
 
    if( njets4 && nbjets2 ) h_nEvents_inclusive->Fill(0);
    if( njets4 && nbjets2 && ttbb  ) h_nEvents_inclusive->Fill(1);
    if( dil && nbjets2 && njets4  ) h_nEvents_inclusive->Fill(2);
    if( dil && nbjets2 && njets4 && ttbb  ) h_nEvents_inclusive->Fill(3);
    if( vis && njets4  ) h_nEvents_inclusive->Fill(4);
    if( vis && njets4 && ttbb  ) h_nEvents_inclusive->Fill(5);

    if( njets4 && nbpartons2 ) h_nEvents_parton_inclusive->Fill(0);
    if( njets4 && nbpartons2 && ttbb_parton  ) h_nEvents_parton_inclusive->Fill(1);
    if( dil && nbpartons2 && njets4  ) h_nEvents_parton_inclusive->Fill(2);
    if( dil && nbpartons2 && njets4 && ttbb_parton  ) h_nEvents_parton_inclusive->Fill(3);
    if( divis && nbpartons2 && njets4  ) h_nEvents_parton_inclusive->Fill(4);
    if( divis && nbpartons2 && njets4 && ttbb_parton  ) h_nEvents_parton_inclusive->Fill(5);
 
    h_multiplicity_GenJets->Fill( ttbarGenLevel.NJets() );
    h_multiplicity_GenJets10->Fill( ttbarGenLevel.NJets10() );
    h_multiplicity_GenJets15->Fill( ttbarGenLevel.NJets15() );
    h_multiplicity_GenJets20->Fill( ttbarGenLevel.NJets20() );
    h_multiplicity_GenJets25->Fill( ttbarGenLevel.NJets25() );
    h_multiplicity_GenJets30->Fill( ttbarGenLevel.NJets30() );

    if( dil ){
      h_multiplicity_GenJetsDIL->Fill( ttbarGenLevel.NJets() );
      h_multiplicity_GenJets10DIL->Fill( ttbarGenLevel.NJets10() );
      h_multiplicity_GenJets15DIL->Fill( ttbarGenLevel.NJets15() );
      h_multiplicity_GenJets20DIL->Fill( ttbarGenLevel.NJets20() );
      h_multiplicity_GenJets25DIL->Fill( ttbarGenLevel.NJets25() );
      h_multiplicity_GenJets30DIL->Fill( ttbarGenLevel.NJets30() );
    }

   if( dil && vis ){
      h_multiplicity_GenJetsDILVIS->Fill( ttbarGenLevel.NJets() );
      h_multiplicity_GenJets10DILVIS->Fill( ttbarGenLevel.NJets10() );
      h_multiplicity_GenJets15DILVIS->Fill( ttbarGenLevel.NJets15() );
      h_multiplicity_GenJets20DILVIS->Fill( ttbarGenLevel.NJets20() );
      h_multiplicity_GenJets25DILVIS->Fill( ttbarGenLevel.NJets25() );
      h_multiplicity_GenJets30DILVIS->Fill( ttbarGenLevel.NJets30() );
    }

    if( dil && vis && ttbb ){
      h_multiplicity_GenJetsDILVISTTBB->Fill( ttbarGenLevel.NJets() );
      h_multiplicity_GenJets10DILVISTTBB->Fill( ttbarGenLevel.NJets10() );
      h_multiplicity_GenJets15DILVISTTBB->Fill( ttbarGenLevel.NJets15() );
      h_multiplicity_GenJets20DILVISTTBB->Fill( ttbarGenLevel.NJets20() );
      h_multiplicity_GenJets25DILVISTTBB->Fill( ttbarGenLevel.NJets25() );
      h_multiplicity_GenJets30DILVISTTBB->Fill( ttbarGenLevel.NJets30() );
    }

    if( dil && vis && ttcc ){
      h_multiplicity_GenJetsDILVISTTCC->Fill( ttbarGenLevel.NJets() );
      h_multiplicity_GenJets10DILVISTTCC->Fill( ttbarGenLevel.NJets10() );
      h_multiplicity_GenJets15DILVISTTCC->Fill( ttbarGenLevel.NJets15() );
      h_multiplicity_GenJets20DILVISTTCC->Fill( ttbarGenLevel.NJets20() );
      h_multiplicity_GenJets25DILVISTTCC->Fill( ttbarGenLevel.NJets25() );
      h_multiplicity_GenJets30DILVISTTCC->Fill( ttbarGenLevel.NJets30() );
    }
  }

  //if( nb_from_top <=2 && nb_from_nontop <=0){
  //  accepted = false;
  //}
  //else {
    //more than 2 b-quarks are from top. But it will be considered as tt+bb.
  //  accepted = true;
  //}

  //b_from_top_multi->Fill(nb_from_top);
  //b_from_nontop_multi->Fill(nb_from_nontop);

  if( ttbarGenLevel.NbJets20(type_) >= 4 ) accepted = true;

  //tree->Fill();

  if(!applyFilter_) return true;  
  else return accepted;

}

bool TTbar2bGenFilter::isLastQuark( const reco::GenParticle& p, const int & pdgId ){
   bool out = true;

   unsigned int nDaughters = p.numberOfDaughters();
   for ( unsigned iDaughter=0; iDaughter<nDaughters; ++iDaughter ) {
     const reco::Candidate* daugh = p.daughter(iDaughter);
     if( abs(daugh->pdgId()) == pdgId) {
       out = false;
       break;
     }
   }

   return out;
}

bool TTbar2bGenFilter::isFromtop( const reco::GenParticle& p){
  bool out = false;
  string tmp = "";
  tmp += "Let's study this b quark \n";
  string pt = Form("%f", p.pt());
  string pdgid = Form("%i",p.pdgId());
  tmp += "pt = " + pt + " id= " + pdgid + "\n";
  const reco::GenParticle* mother = dynamic_cast<const reco::GenParticle*>(p.mother());
  while( mother != 0 ){
    string id = Form("%i", mother->pdgId());
    string mopt = Form("%f", mother->pt());
    tmp += "mother pdgid= " + id + " pt= " + mopt +"\n";
    if( abs(mother->pdgId()) == 6 ) { 
      out = true;
    }
    mother = dynamic_cast<const reco::GenParticle*>(mother->mother());
  }
 
  if(out) debug += tmp; 
  return out;
}

DEFINE_FWK_MODULE(TTbar2bGenFilter);
