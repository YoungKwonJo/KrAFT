#ifndef KrAFT_GenericNtuple_KLeptonJetTreeReducer_H
#define KrAFT_GenericNtuple_KLeptonJetTreeReducer_H

#include "KrAFT/GenericNtuple/interface/KFlatTreeReducerBase.h"

class KLeptonJetTreeReducer : public KFlatTreeReducerBase
{
public:
  KLeptonJetTreeReducer(const std::string modeName,
                       const std::string inputFileName,
                       const std::string outputFileName);
  void run() { KFlatTreeReducerBase::run(); }

private:
  virtual bool analyze();

private:
  // Cache for Leptons, depend on decay modes
  typedef doubles* doublesP;
  typedef ints* intsP;

private:
  // Branch items for output tree
  double lepton_pt_, lepton_eta_, lepton_phi_, lepton_iso_;

  doublesP jets_pt_, jets_bTag_;
  doublesP jetsUp_pt_, jetsUp_bTag_;
  doublesP jetsDn_pt_, jetsDn_bTag_;
  doublesP jetsResUp_pt_, jetsResUp_bTag_;
  doublesP jetsResDn_pt_, jetsResDn_bTag_;

};

#endif
