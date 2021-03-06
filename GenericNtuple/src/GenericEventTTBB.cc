#include "KrAFT/GenericNtuple/interface/GenericEventTTBB.h"

GenericEventTTBB::GenericEventTTBB(bool isMC)
{
  isMC_ = isMC;

  muons_pt_  = new doubles;
  muons_eta_ = new doubles;
  muons_phi_ = new doubles;
  muons_m_   = new doubles;
  muons_Q_      = new ints;
  muons_type_   = new uints;
  muons_relIso_ = new doubles;

  electrons_pt_  = new doubles;
  electrons_eta_ = new doubles;
  electrons_phi_ = new doubles;
  electrons_m_   = new doubles;
  electrons_Q_      = new ints;
  electrons_type_   = new uints;
  electrons_relIso_ = new doubles;

  electrons_mva_ = new doubles;
  electrons_scEta_ = new doubles;
  electrons_qConsistent_ = new uints;

  jets_pt_  = new doubles;
  jets_eta_ = new doubles;
  jets_phi_ = new doubles;
  jets_m_   = new doubles;
  jets_bTag_ = new doubles;
  jets_partonflavor_ = new doubles;
  jets_JESUp_ = new doubles;
  jets_JESDn_ = new doubles;

  jpsis_pt_  = new doubles;
  jpsis_eta_ = new doubles;
  jpsis_phi_ = new doubles;
  jpsis_m_   = new doubles;
  jpsis_lxy_ = new doubles;

  jpsis_pt1_  = new doubles;
  jpsis_eta1_ = new doubles;
  jpsis_phi1_ = new doubles;
  jpsis_pt2_  = new doubles;
  jpsis_eta2_ = new doubles;
  jpsis_phi2_ = new doubles;

  jpsis_nPixHits1_ = new ints;
  jpsis_nPixHits2_ = new ints;

  if ( isMC_ )
  {
    pdfWeights_ = new doubles;

    // JER
    jets_JER_   = new doubles;
    jets_JERUp_ = new doubles;
    jets_JERDn_ = new doubles;

    // GenJets
    genJets_pt_  = new doubles;
    genJets_eta_ = new doubles;
    genJets_phi_ = new doubles;
    genJets_m_   = new doubles;

    // Generator information
    genParticles_pt_  = new doubles;
    genParticles_eta_ = new doubles;
    genParticles_phi_ = new doubles;
    genParticles_m_   = new doubles;
    genParticles_pdgId_ = new ints;
    genParticles_mother_ = new ints;
    genJets_decayFromBHadron_ = new ints;
    genJets_decayFromCHadron_ = new ints;

    //matching genjet for b/c flavor
    bJets_ = new XYZTLorentzVectors;
    cJets_ = new XYZTLorentzVectors;
    bIDs_ = new ints;
    cIDs_ = new ints;

    bpIDs_ = new ints;
    cpIDs_ = new ints;
    bsDRs_ = new doubles;
    csDRs_ = new doubles;

    bprogenitor_pdgId_ = new ints;
    cprogenitor_pdgId_ = new ints;

    bprogenitor_pt_ = new doubles;
    bprogenitor_eta_ = new doubles;
    bprogenitor_phi_ = new doubles;
    cprogenitor_pt_ = new doubles;
    cprogenitor_eta_ = new doubles;
    cprogenitor_phi_ = new doubles;

  }
}

void GenericEventTTBB::book(TTree* tree)
{
  tree_ = tree;

  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");
  tree_->Branch("event", &event_, "event/I");

  tree_->Branch("puWeight", &puWeight_, "puWeight/D");
  tree_->Branch("puWeightUp", &puWeightUp_, "puWeightUp/D");
  tree_->Branch("puWeightDn", &puWeightDn_, "puWeightDn/D");
  tree_->Branch("nVertex", &nVertex_, "nVertex/I");
  tree_->Branch("nPileup", &nPileup_, "nPileup/I");

  tree_->Branch("muons_pt"  , muons_pt_  );
  tree_->Branch("muons_eta" , muons_eta_ );
  tree_->Branch("muons_phi" , muons_phi_ );
  tree_->Branch("muons_m"   , muons_m_   );
  tree_->Branch("muons_Q"   , muons_Q_   );
  tree_->Branch("muons_type", muons_type_);
  tree_->Branch("muons_relIso", muons_relIso_);

  tree_->Branch("electrons_pt"  , electrons_pt_  );
  tree_->Branch("electrons_eta" , electrons_eta_ );
  tree_->Branch("electrons_phi" , electrons_phi_ );
  tree_->Branch("electrons_m"   , electrons_m_   );
  tree_->Branch("electrons_Q"   , electrons_Q_   );
  tree_->Branch("electrons_type", electrons_type_);
  tree_->Branch("electrons_relIso", electrons_relIso_);

  tree_->Branch("electrons_mva", electrons_mva_);
  tree_->Branch("electrons_scEta", electrons_scEta_);
  tree_->Branch("electrons_qConsistent", electrons_qConsistent_);

  tree_->Branch("jets_pt"  , jets_pt_  );
  tree_->Branch("jets_eta" , jets_eta_ );
  tree_->Branch("jets_phi" , jets_phi_ );
  tree_->Branch("jets_m"   , jets_m_   );
  tree_->Branch("jets_bTag", jets_bTag_);
  tree_->Branch("jets_partonflavor", jets_partonflavor_);

  tree_->Branch("jets_JESUp", jets_JESUp_);
  tree_->Branch("jets_JESDn", jets_JESDn_);

  tree_->Branch("met_pt"     , &met_pt_     , "met_pt/D"     );
  tree_->Branch("metJESUp_pt", &metJESUp_pt_, "metJESUp_pt/D");
  tree_->Branch("metJESDn_pt", &metJESDn_pt_, "metJESDn_pt/D");

  tree_->Branch("met_phi"     , &met_phi_     , "met_phi/D"     );
  tree_->Branch("metJESUp_phi", &metJESUp_phi_, "metJESUp_phi/D");
  tree_->Branch("metJESDn_phi", &metJESDn_phi_, "metJESDn_phi/D");

  tree_->Branch("jpsis_pt" , jpsis_pt_ );
  tree_->Branch("jpsis_eta", jpsis_eta_);
  tree_->Branch("jpsis_phi", jpsis_phi_);
  tree_->Branch("jpsis_m"  , jpsis_m_  );
  tree_->Branch("jpsis_lxy", jpsis_lxy_);

  tree_->Branch("jpsis_pt1" , jpsis_pt1_ );
  tree_->Branch("jpsis_eta1", jpsis_eta1_);
  tree_->Branch("jpsis_phi1", jpsis_phi1_);
  tree_->Branch("jpsis_pt2" , jpsis_pt2_ );
  tree_->Branch("jpsis_eta2", jpsis_eta2_);
  tree_->Branch("jpsis_phi2", jpsis_phi2_);

  tree_->Branch("jpsis_nPixHits1", jpsis_nPixHits1_);
  tree_->Branch("jpsis_nPixHits2", jpsis_nPixHits2_);

  if ( isMC_ )
  {
    tree_->Branch("pdfWeights", pdfWeights_);
    tree_->Branch("gluon2t_N"     , &gluon2t_N_     , "gluon2t_N/I"     );

    tree_->Branch("jets_JER"  , jets_JER_);
    tree_->Branch("jets_JERUp", jets_JERUp_);
    tree_->Branch("jets_JERDn", jets_JERDn_);

    tree_->Branch("metJER_pt"  , &metJER_pt_  , "metJER_pt/D"  );
    tree_->Branch("metJERUp_pt", &metJERUp_pt_, "metJERUp_pt/D");
    tree_->Branch("metJERDn_pt", &metJERDn_pt_, "metJERDn_pt/D");

    tree_->Branch("metJER_phi"  , &metJER_phi_  , "metJER_phi/D"  );
    tree_->Branch("metJERUp_phi", &metJERUp_phi_, "metJERUp_phi/D");
    tree_->Branch("metJERDn_phi", &metJERDn_phi_, "metJERDn_phi/D");

    tree_->Branch("genWeight", &genWeight_, "genWeight/D");

    tree_->Branch("pdf_id1", &pdf_id1_, "pdf_id1/I");
    tree_->Branch("pdf_id2", &pdf_id2_, "pdf_id2/I");
    tree_->Branch("pdf_x1" , &pdf_x1_ , "pdf_x1/D" );
    tree_->Branch("pdf_x2" , &pdf_x2_ , "pdf_x2/D" );
    tree_->Branch("pdf_q"  , &pdf_q_  , "pdf_q/D"  );

    tree_->Branch("genJets_pt" , genJets_pt_ );
    tree_->Branch("genJets_eta", genJets_eta_);
    tree_->Branch("genJets_phi", genJets_phi_);
    tree_->Branch("genJets_m"  , genJets_m_  );
    tree_->Branch("genJets_decayFromBHadron"  , genJets_decayFromBHadron_  );
    tree_->Branch("genJets_decayFromCHadron"  , genJets_decayFromCHadron_  );

    tree_->Branch("flavorsIndex"  , &flavorsIndex_,"flavorsIndex/i" );

    tree_->Branch("nb"  , &nb_  ,"nb/i"  );
    tree_->Branch("nc"  , &nc_  ,"nc/i"  );
//    tree_->Branch("nbc"  , &nbc_  ,"nbc/i"  );
//    tree_->Branch("ncb"  , &ncb_  ,"ncb/i"  );

//    tree_->Branch("bIDs", &bIDs_,"bIDs/l");
//    tree_->Branch("cIDs", &cIDs_,"cIDs/l");
//    tree_->Branch("gnjet20bf", &gnjet20bf_,"gnjet20bf/I");

    tree_->Branch("bJets" , "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &bJets_ );
    tree_->Branch("cJets" , "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &cJets_ );
//    tree_->Branch("bJets" ,  bJets_);
//    tree_->Branch("cJets" ,  cJets_);
    tree_->Branch("bIDs"  , bIDs_  );
    tree_->Branch("cIDs"  , cIDs_  );

    tree_->Branch("bpIDs"  , bpIDs_  );
    tree_->Branch("cpIDs"  , cpIDs_  );
    tree_->Branch("bsDRs"  , bsDRs_  );
    tree_->Branch("csDRs"  , csDRs_  );

    tree_->Branch("bprogenitor_pdgId"  , bprogenitor_pdgId_  );
    tree_->Branch("cprogenitor_pdgId"  , cprogenitor_pdgId_  );

    tree_->Branch("bprogenitor_pt"  , bprogenitor_pt_  );
    tree_->Branch("bprogenitor_eta"  , bprogenitor_eta_  );
    tree_->Branch("bprogenitor_phi"  , bprogenitor_phi_  );

    tree_->Branch("cprogenitor_pt"  , cprogenitor_pt_  );
    tree_->Branch("cprogenitor_eta"  , cprogenitor_eta_  );
    tree_->Branch("cprogenitor_phi"  , cprogenitor_phi_  );

    tree_->Branch("genParticles_pt" , genParticles_pt_ );
    tree_->Branch("genParticles_eta", genParticles_eta_);
    tree_->Branch("genParticles_phi", genParticles_phi_);
    tree_->Branch("genParticles_m"  , genParticles_m_  );
    tree_->Branch("genParticles_pdgId", genParticles_pdgId_);
    tree_->Branch("genParticles_mother", genParticles_mother_);

    tree_->Branch("genttbarM"      , &genttbarM_      , "genttbarM/D"  );
    tree_->Branch("nGenJet20"      , &nGenJet20_      , "nGenJet20/I"  );
    tree_->Branch("nGenbJet20"     , &nGenbJet20_     , "nGenbJet20/I"  );
    tree_->Branch("nGenaddbJet20"  , &nGenaddbJet20_  , "nGenaddbJet20/I"  );
    tree_->Branch("nGencJet20"     , &nGencJet20_     , "nGencJet20/I"  );
    tree_->Branch("nGenJet40"      , &nGenJet40_      , "nGenJet40/I"  );
    tree_->Branch("nGenbJet40"     , &nGenbJet40_     , "nGenbJet40/I"  );
    tree_->Branch("nGenaddbJet40"  , &nGenaddbJet40_  , "nGenaddbJet40/I"  );
    tree_->Branch("nGencJet40"     , &nGencJet40_     , "nGencJet40/I"  );
    tree_->Branch("ttbarGen_dileptonic"  , &ttbarGen_dileptonic_  , "ttbarGen_dileptonic/I"  );
    tree_->Branch("genLep1_pt"      , &genLep1_pt_      , "genLep1_pt/D"  );
    tree_->Branch("genLep2_pt"      , &genLep2_pt_      , "genLep2_pt/D"  );
    tree_->Branch("genLep1_eta"      , &genLep1_eta_      , "genLep1_eta/D"  );
    tree_->Branch("genLep2_eta"      , &genLep2_eta_      , "genLep2_eta/D"  );
    tree_->Branch("genLep1_phi"      , &genLep1_phi_      , "genLep1_phi/D"  );
    tree_->Branch("genLep2_phi"      , &genLep2_phi_      , "genLep2_phi/D"  );


  }
}

void GenericEventTTBB::clear()
{
  // Clear up
  electrons_pt_->clear();
  electrons_eta_->clear();
  electrons_phi_->clear();
  electrons_m_->clear();
  electrons_Q_->clear();
  electrons_type_->clear();
  electrons_relIso_->clear();

  electrons_mva_->clear();
  electrons_scEta_->clear();
  electrons_qConsistent_->clear();

  muons_pt_->clear();
  muons_eta_->clear();
  muons_phi_->clear();
  muons_m_->clear();
  muons_Q_->clear();
  muons_type_->clear();
  muons_relIso_->clear();

  jets_pt_->clear();
  jets_eta_->clear();
  jets_phi_->clear();
  jets_m_->clear();
  jets_bTag_->clear();
  jets_partonflavor_->clear();

  jets_JESUp_->clear();
  jets_JESDn_->clear();

  jpsis_pt_->clear();
  jpsis_eta_->clear();
  jpsis_phi_->clear();
  jpsis_m_->clear();
  jpsis_lxy_->clear();

  jpsis_pt1_ ->clear();
  jpsis_eta1_->clear();
  jpsis_phi1_->clear();
  jpsis_pt2_ ->clear();
  jpsis_eta2_->clear();
  jpsis_phi2_->clear();

  jpsis_nPixHits1_->clear();
  jpsis_nPixHits2_->clear();

  if ( isMC_ )
  {
    pdfWeights_->clear();

    jets_JER_->clear();
    jets_JERUp_->clear();
    jets_JERDn_->clear();

    genJets_pt_ ->clear();
    genJets_eta_->clear();
    genJets_phi_->clear();
    genJets_m_  ->clear();
    genJets_decayFromBHadron_->clear();
    genJets_decayFromCHadron_->clear();

    genParticles_pt_ ->clear();
    genParticles_eta_->clear();
    genParticles_phi_->clear();
    genParticles_m_  ->clear();
    genParticles_pdgId_->clear();
    genParticles_mother_->clear();

    bJets_->clear();
    cJets_->clear();
    bIDs_->clear();
    cIDs_->clear();

    bpIDs_->clear();
    cpIDs_->clear();
    bsDRs_->clear();
    csDRs_->clear();

    bprogenitor_pdgId_->clear();
    cprogenitor_pdgId_->clear();

    bprogenitor_pt_->clear();
    bprogenitor_eta_->clear();
    bprogenitor_phi_->clear();

    cprogenitor_pt_->clear();
    cprogenitor_eta_->clear();
    cprogenitor_phi_->clear();
  }
}

void GenericEventTTBB::setBranch(TTree* tree)
{
  tree_ = tree;

  tree_->SetBranchAddress("run", &run_);
  tree_->SetBranchAddress("lumi", &lumi_);
  tree_->SetBranchAddress("event", &event_);

  tree_->SetBranchAddress("puWeight", &puWeight_);
  tree_->SetBranchAddress("puWeightUp", &puWeightUp_);
  tree_->SetBranchAddress("puWeightDn", &puWeightDn_);
  tree_->SetBranchAddress("nVertex", &nVertex_);
  tree_->SetBranchAddress("nPileup", &nPileup_);

  tree_->SetBranchAddress("muons_pt"  , &muons_pt_  );
  tree_->SetBranchAddress("muons_eta" , &muons_eta_ );
  tree_->SetBranchAddress("muons_phi" , &muons_phi_ );
  tree_->SetBranchAddress("muons_m"   , &muons_m_   );
  tree_->SetBranchAddress("muons_Q"   , &muons_Q_   );
  tree_->SetBranchAddress("muons_type", &muons_type_);
  tree_->SetBranchAddress("muons_relIso", &muons_relIso_);

  tree_->SetBranchAddress("electrons_pt"  , &electrons_pt_  );
  tree_->SetBranchAddress("electrons_eta" , &electrons_eta_ );
  tree_->SetBranchAddress("electrons_phi" , &electrons_phi_ );
  tree_->SetBranchAddress("electrons_m"   , &electrons_m_   );
  tree_->SetBranchAddress("electrons_Q"   , &electrons_Q_   );
  tree_->SetBranchAddress("electrons_type", &electrons_type_);
  tree_->SetBranchAddress("electrons_relIso", &electrons_relIso_);

  tree_->SetBranchAddress("electrons_mva", &electrons_mva_);
  tree_->SetBranchAddress("electrons_scEta", &electrons_scEta_);
  tree_->SetBranchAddress("electrons_qConsistent", &electrons_qConsistent_);

  tree_->SetBranchAddress("jets_pt"  , &jets_pt_  );
  tree_->SetBranchAddress("jets_eta" , &jets_eta_ );
  tree_->SetBranchAddress("jets_phi" , &jets_phi_ );
  tree_->SetBranchAddress("jets_m"   , &jets_m_   );
  tree_->SetBranchAddress("jets_bTag", &jets_bTag_);
  tree_->SetBranchAddress("jets_partonflavor", &jets_partonflavor_);

  tree_->SetBranchAddress("jets_JESUp", &jets_JESUp_);
  tree_->SetBranchAddress("jets_JESDn", &jets_JESDn_);

  tree_->SetBranchAddress("met_pt"     , &met_pt_     );
  tree_->SetBranchAddress("metJESUp_pt", &metJESUp_pt_);
  tree_->SetBranchAddress("metJESDn_pt", &metJESDn_pt_);

  tree_->SetBranchAddress("met_phi"     , &met_phi_     );
  tree_->SetBranchAddress("metJESUp_phi", &metJESUp_phi_);
  tree_->SetBranchAddress("metJESDn_phi", &metJESDn_phi_);

  tree_->SetBranchAddress("jpsis_pt" , &jpsis_pt_ );
  tree_->SetBranchAddress("jpsis_eta", &jpsis_eta_);
  tree_->SetBranchAddress("jpsis_phi", &jpsis_phi_);
  tree_->SetBranchAddress("jpsis_m"  , &jpsis_m_  );
  tree_->SetBranchAddress("jpsis_lxy", &jpsis_lxy_);

  tree_->SetBranchAddress("jpsis_pt1" , &jpsis_pt1_ );
  tree_->SetBranchAddress("jpsis_eta1", &jpsis_eta1_);
  tree_->SetBranchAddress("jpsis_phi1", &jpsis_phi1_);
  tree_->SetBranchAddress("jpsis_pt2" , &jpsis_pt2_ );
  tree_->SetBranchAddress("jpsis_eta2", &jpsis_eta2_);
  tree_->SetBranchAddress("jpsis_phi2", &jpsis_phi2_);

  tree_->SetBranchAddress("jpsis_nPixHits1", &jpsis_nPixHits1_);
  tree_->SetBranchAddress("jpsis_nPixHits2", &jpsis_nPixHits2_);

  if ( isMC_ )
  {
    tree_->SetBranchAddress("pdfWeights", &pdfWeights_);
    tree_->SetBranchAddress("gluon2t_N"     , &gluon2t_N_     );

    tree_->SetBranchAddress("jets_JER"  , &jets_JER_  );
    tree_->SetBranchAddress("jets_JERUp", &jets_JERDn_);
    tree_->SetBranchAddress("jets_JERDn", &jets_JERDn_);

    tree_->SetBranchAddress("metJER_pt"  , &metJER_pt_  );
    tree_->SetBranchAddress("metJERUp_pt", &metJERUp_pt_);
    tree_->SetBranchAddress("metJERDn_pt", &metJERDn_pt_);

    tree_->SetBranchAddress("metJER_phi"  , &metJER_phi_  );
    tree_->SetBranchAddress("metJERUp_phi", &metJERUp_phi_);
    tree_->SetBranchAddress("metJERDn_phi", &metJERDn_phi_);

    tree_->SetBranchAddress("genWeight", &genWeight_);

    tree_->SetBranchAddress("pdf_id1", &pdf_id1_);
    tree_->SetBranchAddress("pdf_id2", &pdf_id2_);
    tree_->SetBranchAddress("pdf_x1" , &pdf_x1_ );
    tree_->SetBranchAddress("pdf_x2" , &pdf_x2_ );
    tree_->SetBranchAddress("pdf_q"  , &pdf_q_  );

    tree_->SetBranchAddress("genJets_pt" , &genJets_pt_ );
    tree_->SetBranchAddress("genJets_eta", &genJets_eta_);
    tree_->SetBranchAddress("genJets_phi", &genJets_phi_);
    tree_->SetBranchAddress("genJets_m"  , &genJets_m_  );
    tree_->SetBranchAddress("genJets_decayFromBHadron", &genJets_decayFromBHadron_);
    tree_->SetBranchAddress("genJets_decayFromCHadron", &genJets_decayFromCHadron_);

    tree_->SetBranchAddress("flavorsIndex", &flavorsIndex_);
    tree_->SetBranchAddress("nb", &nb_);
    tree_->SetBranchAddress("nc", &nc_);

//    tree_->SetBranchAddress("nbc", &nbc_);
//    tree_->SetBranchAddress("ncb", &ncb_);

//    tree_->SetBranchAddress("bIDs", &bIDs_);
//    tree_->SetBranchAddress("cIDs", &cIDs_);
//    tree_->SetBranchAddress("gnjet20bf", &gnjet20bf_);

    tree_->SetBranchAddress("bJets", &bJets_);
    tree_->SetBranchAddress("cJets", &cJets_);
    tree_->SetBranchAddress("bIDs" , &bIDs_ );
    tree_->SetBranchAddress("cIDs" , &cIDs_ );

    tree_->SetBranchAddress("bpIDs" , &bpIDs_ );
    tree_->SetBranchAddress("cpIDs" , &cpIDs_ );
    tree_->SetBranchAddress("bsDRs" , &bsDRs_ );
    tree_->SetBranchAddress("csDRs" , &csDRs_ );

    tree_->SetBranchAddress("bprogenitor_pdgId", &bprogenitor_pdgId_);
    tree_->SetBranchAddress("cprogenitor_pdgId", &cprogenitor_pdgId_);

    tree_->SetBranchAddress("bprogenitor_pt", &bprogenitor_pt_);
    tree_->SetBranchAddress("bprogenitor_eta", &bprogenitor_eta_);
    tree_->SetBranchAddress("bprogenitor_phi", &bprogenitor_phi_);

    tree_->SetBranchAddress("cprogenitor_pt", &cprogenitor_pt_);
    tree_->SetBranchAddress("cprogenitor_eta", &cprogenitor_eta_);
    tree_->SetBranchAddress("cprogenitor_phi", &cprogenitor_phi_);

    tree_->SetBranchAddress("genParticles_pt" , &genParticles_pt_ );
    tree_->SetBranchAddress("genParticles_eta", &genParticles_eta_);
    tree_->SetBranchAddress("genParticles_phi", &genParticles_phi_);
    tree_->SetBranchAddress("genParticles_m"  , &genParticles_m_  );
    tree_->SetBranchAddress("genParticles_pdgId", &genParticles_pdgId_);
    tree_->SetBranchAddress("genParticles_mother", &genParticles_mother_);

    tree_->SetBranchAddress("genttbarM"      , &genttbarM_  );
    tree_->SetBranchAddress("nGenbJet20"     , &nGenbJet20_  );
    tree_->SetBranchAddress("nGenaddbJet20"  , &nGenaddbJet20_  );
    tree_->SetBranchAddress("nGencJet20"     , &nGencJet20_  );
    tree_->SetBranchAddress("nGenbJet40"     , &nGenbJet40_  );
    tree_->SetBranchAddress("nGenaddbJet40"  , &nGenaddbJet40_  );
    tree_->SetBranchAddress("nGencJet40"     , &nGencJet40_  );
    tree_->SetBranchAddress("ttbarGen_dileptonic"  , &ttbarGen_dileptonic_  );


    tree_->SetBranchAddress("genLep1_pt"     , &genLep1_pt_  );
    tree_->SetBranchAddress("genLep2_pt"     , &genLep2_pt_  );
    tree_->SetBranchAddress("genLep1_eta"     , &genLep1_eta_  );
    tree_->SetBranchAddress("genLep2_eta"     , &genLep2_eta_  );
    tree_->SetBranchAddress("genLep1_phi"     , &genLep1_phi_  );
    tree_->SetBranchAddress("genLep2_phi"     , &genLep2_phi_  );
 
  }
}

GenericEventTTBB::~GenericEventTTBB()
{
  delete muons_pt_ ;
  delete muons_eta_;
  delete muons_phi_;
  delete muons_m_  ;
  delete muons_Q_   ;
  delete muons_type_;
  delete muons_relIso_;

  delete electrons_pt_ ;
  delete electrons_eta_;
  delete electrons_phi_;
  delete electrons_m_  ;
  delete electrons_Q_   ;
  delete electrons_type_;
  delete electrons_relIso_;

  delete electrons_mva_;
  delete electrons_scEta_;
  delete electrons_qConsistent_;

  delete jets_pt_ ;
  delete jets_eta_;
  delete jets_phi_;
  delete jets_m_  ;
  delete jets_bTag_;
  delete jets_partonflavor_;
  delete jets_JESUp_;
  delete jets_JESDn_;

  delete jpsis_pt_ ;
  delete jpsis_eta_;
  delete jpsis_phi_;
  delete jpsis_m_  ;
  delete jpsis_lxy_;

  delete jpsis_pt1_ ;
  delete jpsis_eta1_;
  delete jpsis_phi1_;
  delete jpsis_pt2_ ;
  delete jpsis_eta2_;
  delete jpsis_phi2_;

  delete jpsis_nPixHits1_;
  delete jpsis_nPixHits2_;

  if ( isMC_ )
  {
    delete pdfWeights_;

    // JER
    delete jets_JER_  ;
    delete jets_JERUp_;
    delete jets_JERDn_;

    // GenJets
    delete genJets_pt_ ;
    delete genJets_eta_;
    delete genJets_phi_;
    delete genJets_m_  ;
    delete genJets_decayFromBHadron_;
    delete genJets_decayFromCHadron_;

    // Generator information
    delete genParticles_pt_ ;
    delete genParticles_eta_;
    delete genParticles_phi_;
    delete genParticles_m_  ;
    delete genParticles_pdgId_ ;
    delete genParticles_mother_;

    delete bJets_;
    delete cJets_;
    delete bIDs_;
    delete cIDs_;

    delete bpIDs_;
    delete cpIDs_;
    delete bsDRs_;
    delete csDRs_;

    delete bprogenitor_pdgId_;
    delete cprogenitor_pdgId_;

    delete bprogenitor_pt_;
    delete bprogenitor_eta_;
    delete bprogenitor_phi_;

    delete cprogenitor_pt_;
    delete cprogenitor_eta_;
    delete cprogenitor_phi_;
  }
}
