#include "TSystem.h"
#include "TInterpreter.h"
#include "TList.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "../interface/convertPythia.h"

using namespace PDGHelpers;

convertPythia::convertPythia(OptionParser* options_) : converter(options_){
  configure();
  run();
}

void convertPythia::configure(){
  string tmpdir = options->getTempDir();
  string strCmd = "mkdir -p ";
  strCmd.append(tmpdir);
  gSystem->Exec(strCmd.c_str());
}
void convertPythia::finalizeRun(){
  string tmpdir = options->getTempDir();
  string strCmd = "rm -rf ";
  strCmd.append(tmpdir);
  gSystem->Exec(strCmd.c_str());
}
void convertPythia::run(){
  Float_t MC_weight=0;
  Int_t isSelected=0;

  int globalNEvents = 0;
  int maxProcEvents = options->maxEventsToProcess();
  vector < pair<Int_t, Int_t> > eventSkipList = options->getSkippedEvents();

  tree->bookAllBranches(false);

  for (int f=0; f<filename.size(); f++){
    string cinput = filename.at(f);
    cout << "Processing " << cinput << "..." << endl;
    TFile* fin = 0;
    fin = getIntermediateFile(cinput);
    if (fin!=0 && fin->IsZombie()){
      if (fin->IsOpen()) fin->Close();
      delete fin;
      fin=0;
    }
    else if (fin!=0 && !fin->IsOpen()){ // Separate if-statement on purpose
      delete fin;
      fin=0;
    }
    else if (fin!=0){
      TTree* tin = (TTree*)fin->Get("tmpTree");
      foutput->cd();

      int nInputEvents = tin->GetEntries();
      int nProcessed = 0;
      for (int ev=0; ev<nInputEvents; ev++){
        if (globalNEvents>=maxProcEvents && maxProcEvents>=0) break;
        bool doSkipEvent = false;
        for (int es=0; es<eventSkipList.size(); es++){
          if (
            (eventSkipList.at(es).first<=globalNEvents && eventSkipList.at(es).second>=globalNEvents)
            ||
            (eventSkipList.at(es).first<=globalNEvents && eventSkipList.at(es).second<0)
            )doSkipEvent=true;
        }
        if (doSkipEvent){ globalNEvents++; continue; }

        double weight;
        bool genSuccess=false, smearedSuccess=false;

        vector<Particle*> genParticleList;
        vector<Particle*> smearedParticleList;
        vector<ZZCandidate*> genCandList; // Bookkeeping
        vector<ZZCandidate*> smearedCandList; // Bookkeeping

        tree->initializeBranches();

        readEvent(tin, ev, genParticleList, genSuccess, smearedParticleList, smearedSuccess, weight);

        if (weight!=0){
          MC_weight = (float)weight;

          Event genEvent;
          if (genSuccess){
            genEvent.setWeight(weight);
            vectorInt hasGenHiggs;
            for (int p=0; p<genParticleList.size(); p++){
              Particle* genPart = genParticleList.at(p); // Has mother info from Pythia reading
              if (isAHiggs(genPart->id)){
                hasGenHiggs.push_back(p);
                if (options->doGenHZZdecay()==-1 && (genPart->genStatus==1 || genPart->genStatus==2)) genEvent.addIntermediate(genPart);
              }
              if (genPart->genStatus==1){
                if (isALepton(genPart->id)) genEvent.addLepton(genPart);
                else if (isANeutrino(genPart->id)) genEvent.addNeutrino(genPart);
                else if (isAGluon(genPart->id) || isAQuark(genPart->id)) genEvent.addJet(genPart);
              }
              else if (genPart->genStatus==-1) tree->fillMotherInfo(genPart);
            }

            genEvent.constructVVCandidates(options->doGenHZZdecay(), options->genDecayProducts());
            for (int p=0; p<genParticleList.size(); p++){
              Particle* genPart = genParticleList.at(p);
              if (genPart->genStatus==-1) genEvent.addVVCandidateMother(genPart);
            }
            genEvent.addVVCandidateAppendages();
            ZZCandidate* genCand=0;
            if (hasGenHiggs.size()>0){
              for (int gk=0; gk<hasGenHiggs.size(); gk++){
                ZZCandidate* tmpCand = HiggsComparators::matchAHiggsToParticle(genEvent, genParticleList.at(hasGenHiggs.at(gk)));
                if (tmpCand!=0){
                  if (genCand==0) genCand=tmpCand;
                  else genCand = HiggsComparators::candComparator(genCand, tmpCand, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
                }
              }
            }
            else genCand = HiggsComparators::candidateSelector(genEvent, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
            if (genCand!=0) tree->fillCandidate(genCand, true);
            else{
              cout << cinput << " (" << ev << "): No gen. level Higgs candidate was found!" << endl;
              genSuccess=false;
            }
          }

          Event smearedEvent;
          if (smearedSuccess){
            smearedEvent.setWeight(weight);
            for (int p=0; p<smearedParticleList.size(); p++){
              Particle* smearedPart = smearedParticleList.at(p);
              if (isALepton(smearedPart->id)) smearedEvent.addLepton(smearedPart);
              else if (isANeutrino(smearedPart->id)) smearedEvent.addNeutrino(smearedPart);
              else if (smearedPart->id==0) smearedEvent.addJet(smearedPart);
              else smearedEvent.addParticle(smearedPart);
            }
            smearedEvent.constructVVCandidates(options->doRecoHZZdecay(), options->recoDecayProducts());
            if (options->recoSelectionMode()==0) smearedEvent.applyParticleSelection();
            smearedEvent.addVVCandidateAppendages();
            ZZCandidate* rCand = HiggsComparators::candidateSelector(smearedEvent, options->getHiggsCandidateSelectionScheme(false), options->doRecoHZZdecay());
            if (rCand!=0){
              isSelected=1;
              tree->fillCandidate(rCand, false);
            }
            else{
              isSelected=0;
              smearedSuccess=false;
            }
          }

          tree->fillEventVariables(MC_weight, isSelected);
          if ((smearedSuccess && options->processRecoInfo()) || (genSuccess && options->processGenInfo())){
            tree->record();
            nProcessed++;
          }
        }

        for (int p=0; p<smearedCandList.size(); p++){ // Bookkeeping
          ZZCandidate* tmpCand = (ZZCandidate*)smearedCandList.at(p);
          if (tmpCand!=0) delete tmpCand;
        }
        for (int p=0; p<smearedParticleList.size(); p++){ // Bookkeeping
          Particle* tmpPart = (Particle*)smearedParticleList.at(p);
          if (tmpPart!=0) delete tmpPart;
        }

        for (int p=0; p<genCandList.size(); p++){ // Bookkeeping
          ZZCandidate* tmpCand = (ZZCandidate*)genCandList.at(p);
          if (tmpCand!=0) delete tmpCand;
        }
        for (int p=0; p<genParticleList.size(); p++){ // Bookkeeping
          Particle* tmpPart = (Particle*)genParticleList.at(p);
          if (tmpPart!=0) delete tmpPart;
        }

        // Bookkeeping
        smearedCandList.clear();
        smearedParticleList.clear();
        genCandList.clear();
        genParticleList.clear();

        globalNEvents++;
      }
      fin->Close();
      cout << "Processed number of events from the input file (recorded events / sample size observed / cumulative traversed): " << nProcessed << " / " << nInputEvents << " / " << globalNEvents << endl;
    }
  }
  finalizeRun();
}


TFile* convertPythia::getIntermediateFile(string cinput){
  string coutput = options->getTempDir();
  string strCmd = "root -b -l -q loadLib.cc 'trimPythia.cc+(\"";
  strCmd.append(cinput);
  strCmd.append("\", \"");
  strCmd.append(coutput);
  strCmd.append("\", ");
  strCmd.append(to_string(options->analysisLevel()));
  strCmd.append(")'");
  gSystem->Exec(strCmd.c_str());
  string strtmp=coutput;
  strtmp.append("pythiaTemp.root");
  TFile* ftmp = new TFile(strtmp.c_str(), "read");
  return ftmp;
}


void convertPythia::readEvent(TTree* tin, int ev, vector<Particle*>& genCollection, bool& genSuccess, vector<Particle*>& recoCollection, bool& smearedSuccess, double& eventWeight){
  int nEvents = tin->GetEntries();
  vectorDouble weights;
  if (ev>=nEvents){
    double weight = 0;
    genSuccess=false;
    smearedSuccess=false;
    weights.push_back(0);
  }
  else{
    vectorDouble* geneventinfoweights=0;
    TBranch* b_geneventinfoweights=0;

    vectorDouble* reco_GenParticle_FV[4]={ 0 };
    vectorInt* reco_GenParticle_id=0;
    vectorInt* reco_GenParticle_status=0;
    TBranch* b_reco_GenParticle_FV[4]={ 0 };
    TBranch* b_reco_GenParticle_id=0;
    TBranch* b_reco_GenParticle_status=0;

    vectorDouble* reco_GenJet_FV[4]={ 0 };
    vectorInt* reco_GenJet_id=0;
    vectorInt* reco_GenJet_status=0;
    TBranch* b_reco_GenJet_FV[4]={ 0 };
    TBranch* b_reco_GenJet_id=0;
    TBranch* b_reco_GenJet_status=0;

    tin->SetBranchAddress("genWeights", &geneventinfoweights, &b_geneventinfoweights);

    tin->SetBranchAddress("reco_GenParticle_X", reco_GenParticle_FV, b_reco_GenParticle_FV);
    tin->SetBranchAddress("reco_GenParticle_Y", reco_GenParticle_FV+1, b_reco_GenParticle_FV+1);
    tin->SetBranchAddress("reco_GenParticle_Z", reco_GenParticle_FV+2, b_reco_GenParticle_FV+2);
    tin->SetBranchAddress("reco_GenParticle_E", reco_GenParticle_FV+3, b_reco_GenParticle_FV+3);
    tin->SetBranchAddress("reco_GenParticle_id", &reco_GenParticle_id, &b_reco_GenParticle_id);
    tin->SetBranchAddress("reco_GenParticle_status", &reco_GenParticle_status, &b_reco_GenParticle_status);

    tin->SetBranchAddress("reco_GenJet_X", reco_GenJet_FV, b_reco_GenJet_FV);
    tin->SetBranchAddress("reco_GenJet_Y", reco_GenJet_FV+1, b_reco_GenJet_FV+1);
    tin->SetBranchAddress("reco_GenJet_Z", reco_GenJet_FV+2, b_reco_GenJet_FV+2);
    tin->SetBranchAddress("reco_GenJet_E", reco_GenJet_FV+3, b_reco_GenJet_FV+3);
    tin->SetBranchAddress("reco_GenJet_id", &reco_GenJet_id, &b_reco_GenJet_id);
    tin->SetBranchAddress("reco_GenJet_status", &reco_GenJet_status, &b_reco_GenJet_status);

    tin->GetEntry(ev);

    // Gen. particles
    int motherID[2];
    int mctr=0;
    for (int a = 0; a < reco_GenParticle_id->size(); a++){
      int istup = reco_GenParticle_status->at(a);
      if (istup==21 && mctr<2){
        motherID[mctr] = a;
        mctr++;
      }
      int idup = reco_GenParticle_id->at(a);
      TLorentzVector partFourVec(reco_GenParticle_FV[0]->at(a), reco_GenParticle_FV[1]->at(a), reco_GenParticle_FV[2]->at(a), reco_GenParticle_FV[3]->at(a));

      Particle* onePart = new Particle(idup, partFourVec);
      onePart->setGenStatus(PDGHelpers::convertPythiaStatus(istup));
      onePart->setLifetime(0);
      genCollection.push_back(onePart);
    }
    // Assign the mothers
    for (int a = 0; a < genCollection.size(); a++){
      for(int m=0;m<2;m++) genCollection.at(a)->addMother(genCollection.at(motherID[m]));
    }
    genSuccess=(genCollection.size()>0);
    // Reco. particles
    for (int a = 0; a < reco_GenJet_id->size(); a++){
      int istup = reco_GenJet_status->at(a);
      int idup = reco_GenJet_id->at(a);
      TLorentzVector partFourVec(reco_GenJet_FV[0]->at(a), reco_GenJet_FV[1]->at(a), reco_GenJet_FV[2]->at(a), reco_GenJet_FV[3]->at(a));

      Particle* onePart = new Particle(idup, partFourVec);
      if (options->recoSmearingMode()>0){ // Reverse of LHE mode
        Particle* smearedPart = LHEParticleSmear::smearParticle(onePart);
        delete onePart;
        onePart = smearedPart;
      }
      onePart->setGenStatus(PDGHelpers::convertPythiaStatus(istup));
      onePart->setLifetime(0);
      recoCollection.push_back(onePart);
    }
    smearedSuccess=(recoCollection.size()>0);

    tin->ResetBranchAddresses();

    if (geneventinfoweights!=0 && geneventinfoweights->size()>0){
      for (int w=0; w<geneventinfoweights->size(); w++) weights.push_back(geneventinfoweights->at(w));
    }
    else if (!smearedSuccess && !genSuccess) weights.push_back(0);
    else weights.push_back(1.);
  }
  eventWeight = weights.at(0);
}

