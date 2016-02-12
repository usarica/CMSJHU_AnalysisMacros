#ifndef OPTIONPARSER_H
#define OPTIONPARSER_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include "TString.h"
#include "PDGHelpers.h"
#include "melaHelpers.h"
#include "ParticleComparators.h"
#include "HiggsComparators.h"

using namespace std;

class OptionParser{
public:
  OptionParser(int argc, char** argv);
  ~OptionParser(){ deconfigureMela(); };

  void analyze();
  void splitOption(string rawoption, string& wish, string& value, char delimiter='=');
  void splitOptionRecursive(string rawoption, vector<string>& splitoptions, char delimiter=',');
  void interpretOption(string wish, string value);
  void printOptionsHelp();
  void printOptionSummary();

  Bool_t processGenInfo(){ bool doProcess=true; if (includeGenInfo==0) doProcess=false; return doProcess; }
  Bool_t processRecoInfo(){ bool doProcess=true; if (includeRecoInfo==0) doProcess=false; return doProcess; }
  Bool_t isAnExcludedBranch(string branchname);
  Int_t analysisLevel(){ return fileLevel; }
  Int_t doGenHZZdecay(){ return isGenHZZ; }
  Int_t doRecoHZZdecay(){ return isRecoHZZ; }
  Int_t genDecayProducts(){ return genDecayMode; }
  Int_t recoDecayProducts(){ return recoDecayMode; }
  Int_t recoSelectionMode(){ return recoSelBehaviour; }
  Int_t recoSmearingMode(){ return recoSmearBehaviour; }
  HiggsComparators::CandidateSelection getHiggsCandidateSelectionScheme(bool isGen=false){ if (isGen) return genHiggsCandidateSelectionScheme; else return recoHiggsCandidateSelectionScheme; }
  string inputDir(){ return indir+"/"; }
  string outputDir(){ return outdir+"/"; }
  string getTempDir(){ return tmpDir; }
  string outputFilename(){ return coutput; }
  vector<string> inputfiles(){ return filename; }
  Int_t maxEventsToProcess(){ return maxEvents; };
  vector < pair<Int_t, Int_t> > getSkippedEvents(){ return eventSkipRanges; };

  // MELA-related options
  Double_t mH(){ return mPOLE; }
  Double_t GammaH(){ return wPOLE; }
  Int_t sqrts(){ return erg_tev; }
  Bool_t initializeMELA(){ return (includeGenDecayProb.size()>0 || includeRecoDecayProb.size()>0 || includeGenProdProb.size()>0 || includeRecoProdProb.size()>0); }
  Bool_t doRemoveLepMasses(){ bool doProcess=true; if (removeDaughterMasses==0) doProcess=false; return doProcess; }
  Bool_t doComputeDecayAngles(){ bool doProcess=true; if (computeDecayAngles==0) doProcess=false; return doProcess; }
  Bool_t doComputeVBFAngles(){ bool doProcess=true; if (computeVBFAngles==0) doProcess=false; return doProcess; }
  Bool_t doComputeVHAngles(){ bool doProcess=true; if (computeVHAngles==0) doProcess=false; return doProcess; }
  Int_t computeVHAnglesOrder(){ return computeVHAngles; }
  Bool_t hasGenDecayME(string str);
  Bool_t hasRecoDecayME(string str);
  Bool_t hasGenProdME(string str);
  Bool_t hasRecoProdME(string str);
  pair<TVar::Production, TVar::MatrixElement> getSampleProductionId(){ return sampleProductionId; }


protected:
  void extractSkippedEvents(string rawoption);

  void configureMela();
  void deconfigureMela();
  void extractMelaGenProdId(string rawoption);

  Bool_t checkListVariable(vector<string>& list, string var);

  vector<string> rawOptions;

  Double_t mPOLE;
  Double_t wPOLE;
  Double_t wPOLEStandard;
  Int_t erg_tev;
  Int_t includeGenInfo;
  Int_t includeRecoInfo;
  Int_t removeDaughterMasses;
  Int_t computeDecayAngles;
  Int_t computeVBFAngles;
  Int_t computeVHAngles;
  Int_t fileLevel;
  Int_t isGenHZZ;
  Int_t isRecoHZZ;
  Int_t genDecayMode;
  Int_t recoDecayMode;
  Int_t recoSelBehaviour;
  Int_t recoSmearBehaviour;
  HiggsComparators::CandidateSelection genHiggsCandidateSelectionScheme;
  HiggsComparators::CandidateSelection recoHiggsCandidateSelectionScheme;

  string indir;
  string outdir;
  string coutput;
  string tmpDir;
  vector<string> filename;
  vector<string> excludedBranch;
  Int_t maxEvents;

  // Mela probabilities to include, has to be in abbreviated form (eg. "All", "None", "p0plus", "g1", "g1_prime2" etc.)
  pair<TVar::Production, TVar::MatrixElement> sampleProductionId;
  vector<string> includeGenDecayProb;
  vector<string> includeRecoDecayProb;
  vector<string> includeGenProdProb;
  vector<string> includeRecoProdProb;

  vector < pair<Int_t, Int_t> > eventSkipRanges;
};

#endif
