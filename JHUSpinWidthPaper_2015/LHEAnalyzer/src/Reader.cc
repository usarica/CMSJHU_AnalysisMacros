#include "TSystem.h"
#include "TInterpreter.h"
#include "TList.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "../interface/Reader.h"

using namespace PDGHelpers;

Reader::Reader(OptionParser* options_) : converter(options_){
  configure();
  run();
}

void Reader::configure(){

}
void Reader::finalizeRun(){}

template<typename returnType> bool Reader::setVariable(const Event* ev, string& branchname, returnType(*evalVar)(const Event*, string&)){
  returnType result = evalVar(ev, branchname);
  returnType* resultPtr = &result;
  if (dynamic_cast<vectorInt*>(resultPtr)!=0 || dynamic_cast<vectorDouble*>(resultPtr)!=0){
    for (int el=0; el<result.size(); el++) tree->setVal(branchname, result.at(el));
    return true;
  }
  else if (dynamic_cast<Int_t*>(resultPtr)!=0 || dynamic_cast<Float_t*>(resultPtr)!=0){
    tree->setVal(branchname, result);
    return true;
  }
  else return false;
}

void Reader::bindInputBranches(HVVTree* tin){
  vector<string> inputBranches = tin->getBranchList();
  vector<pair<string, BaseTree::BranchTypes> > unreservedBranches;
  for (int br=0; br<inputBranches.size(); br++){
    int pos=-1;
    string branchname = inputBranches.at(br);

    BaseTree::BranchTypes iBT = tin->searchArray(branchname, pos);
    pos=-1;
    BaseTree::BranchTypes oBT = tree->searchArray(branchname, pos);
    bool outfound = !(pos==-1 || oBT==BaseTree::nBranchTypes || iBT!=oBT);
    if (!outfound && tree->getBranchList().size()==0) unreservedBranches.push_back(pair<string, BaseTree::BranchTypes>(branchname, iBT));
  }
  if (unreservedBranches.size()>0){
    for (int b=0; b<unreservedBranches.size(); b++) tree->reserveBranch(unreservedBranches.at(b).first, unreservedBranches.at(b).second, false);
    tree->actuateBranches(false);
  }

  for (int br=0; br<inputBranches.size(); br++){
    int pos=-1;
    string branchname = inputBranches.at(br);

    BaseTree::BranchTypes iBT = tin->searchArray(branchname, pos);
    pos=-1;
    BaseTree::BranchTypes oBT = tree->searchArray(branchname, pos);
    bool outfound = !(pos==-1 || oBT==BaseTree::nBranchTypes || iBT!=oBT);
    if (!outfound) continue; // If still not found, forget it!
    else{
      void* iBVRef = tin->getBranchHandleRef(branchname);
      void* oBVRef = tree->getBranchHandleRef(branchname);
      if (iBVRef==0 || oBVRef==0) continue;
      if (iBT == BaseTree::bInt){
        Int_t* iBRef = (Int_t*)iBVRef;
        Int_t* oBRef = (Int_t*)oBVRef;
        pair<Int_t*, Int_t*> refpair(iBRef, oBRef);
        intBranchMap.push_back(refpair);
      }
      else if (iBT == BaseTree::bFloat){
        Float_t* iBRef = (Float_t*)iBVRef;
        Float_t* oBRef = (Float_t*)oBVRef;
        pair<Float_t*, Float_t*> refpair(iBRef, oBRef);
        floatBranchMap.push_back(refpair);
      }
      else if (iBT == BaseTree::bVectorInt){
        vectorInt** iBRef = (vectorInt**)iBVRef;
        vectorInt** oBRef = (vectorInt**)oBVRef;
        pair<vectorInt**, vectorInt**> refpair(iBRef, oBRef);
        vectorIntBranchMap.push_back(refpair);
      }
      else if (iBT == BaseTree::bVectorDouble){
        vectorDouble** iBRef = (vectorDouble**)iBVRef;
        vectorDouble** oBRef = (vectorDouble**)oBVRef;
        pair<vectorDouble**, vectorDouble**> refpair(iBRef, oBRef);
        vectorDoubleBranchMap.push_back(refpair);
      }
    }
  }
}
void Reader::resetBranchBinding(){
  intBranchMap.clear();
  floatBranchMap.clear();
  vectorIntBranchMap.clear();
  vectorDoubleBranchMap.clear();
}
void Reader::synchMappedBranches(){
  for (int b=0; b<intBranchMap.size(); b++) *(intBranchMap.at(b).second) = *(intBranchMap.at(b).first);
  for (int b=0; b<floatBranchMap.size(); b++) *(floatBranchMap.at(b).second) = *(floatBranchMap.at(b).first);
  for (int b=0; b<vectorIntBranchMap.size(); b++){
    vectorInt* inhandle = *(vectorIntBranchMap.at(b).first);
    vectorInt* outhandle = *(vectorIntBranchMap.at(b).second);
    for (int el=0; el<inhandle->size(); el++){
      outhandle->push_back(inhandle->at(el));
    }
  }
  for (int b=0; b<vectorDoubleBranchMap.size(); b++){
    vectorDouble* inhandle = *(vectorDoubleBranchMap.at(b).first);
    vectorDouble* outhandle = *(vectorDoubleBranchMap.at(b).second);
    for (int el=0; el<inhandle->size(); el++){
      outhandle->push_back(inhandle->at(el));
    }
  }
}



void Reader::run(){
  Float_t MC_weight=0;
  Int_t isSelected=0;

  int globalNEvents = 0;
  int maxProcEvents = options->maxEventsToProcess();
  vector < pair<Int_t, Int_t> > eventSkipList = options->getSkippedEvents();

  bool firstFile = true;
  for (int f=0; f<filename.size(); f++){
    string cinput = filename.at(f);
    cout << "Processing " << cinput << "..." << endl;
    TFile* fin = new TFile(cinput.c_str(), "read");
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
      int nProcessed = 0;

      HVVTree* tin = new HVVTree("SelectedTree", fin);
      if (tin->getTree()!=0){
        tin->setOptions(options);
        tin->bookAllBranches(true);
        cout << "Input tree branches booked...";
        bindInputBranches(tin);
        cout << "bound..." << endl;
        foutput->cd();
        if (firstFile){
          tree->bookAllBranches(false);
          firstFile=false;
        }

        int nInputEvents = tin->getTree()->GetEntries();
        cout << "Number of input events to process: " << nInputEvents << endl;
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

          vector<Particle*> genParticleList;
          vector<Particle*> recoParticleList;
          vector<ZZCandidate*> genCandList; // Bookkeeping
          vector<ZZCandidate*> recoCandList; // Bookkeeping

          tree->initializeBranches();

          tin->getTree()->GetEntry(ev);
          synchMappedBranches();

          Event genEvent, recoEvent;
          ZZCandidate* genCand=0;
          ZZCandidate* recoCand=0;
          if (options->processGenInfo()){
            readEvent(genEvent, genParticleList, true);
            if (genEvent.getNZZCandidates()<=1) genCand = genEvent.getZZCandidate(0);
            else genCand = HiggsComparators::candidateSelector(genEvent, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
          }
          if (options->processRecoInfo()){
            readEvent(recoEvent, recoParticleList, false);
            if (recoEvent.getNZZCandidates()<=1) recoCand = recoEvent.getZZCandidate(0);
            else recoCand = HiggsComparators::candidateSelector(recoEvent, options->getHiggsCandidateSelectionScheme(false), options->doRecoHZZdecay());
          }

          if (genCand!=0){
            if (options->doComputeDecayAngles()) tree->fillDecayAngles(genCand, true);
            if (options->doComputeVBFAngles()) tree->fillVBFProductionAngles(genCand, true);
            if (options->doComputeVHAngles()) tree->fillVHProductionAngles(genCand, true);
            if (options->initializeMELA()) tree->fillMELAProbabilities(genCand, true);
          }
          if (recoCand!=0){
            if (options->recoSelectionMode()!=0){
              tree->fillCandidate(recoCand, false);
              tree->fillEventVariables(*((Float_t*)tree->getBranchHandleRef("MC_weight")), 0 /*isSelected*/);
            }
            else{
              if (options->doComputeDecayAngles()) tree->fillDecayAngles(recoCand, false);
              if (options->doComputeVBFAngles()) tree->fillVBFProductionAngles(recoCand, false);
              if (options->doComputeVHAngles()) tree->fillVHProductionAngles(recoCand, false);
              if (options->initializeMELA()) tree->fillMELAProbabilities(recoCand, false);
            }
          }
          else if (options->recoSelectionMode()!=0) tree->fillEventVariables(*((Float_t*)tree->getBranchHandleRef("MC_weight")), 0 /*isSelected*/);

          if ((recoCand!=0 && options->processRecoInfo()) || (genCand!=0 && options->processGenInfo())){
            tree->record();
            nProcessed++;
          }

          for (int p=0; p<recoCandList.size(); p++){ // Bookkeeping
            ZZCandidate* tmpCand = (ZZCandidate*)recoCandList.at(p);
            if (tmpCand!=0) delete tmpCand;
          }
          for (int p=0; p<recoParticleList.size(); p++){ // Bookkeeping
            Particle* tmpPart = (Particle*)recoParticleList.at(p);
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
          recoCandList.clear();
          recoParticleList.clear();
          genCandList.clear();
          genParticleList.clear();

          globalNEvents++;
        }
        resetBranchBinding();
        cout << "Processed number of events from the input file (recorded events / sample size observed / cumulative traversed): " << nProcessed << " / " << nInputEvents << " / " << globalNEvents << endl;
      }
      delete tin;
      fin->Close();
    }
  }
  finalizeRun();
}



void Reader::readEvent(Event& outEvent, vector<Particle*>& particles, bool isGen){
  string varname;

  string isSelected = "isSelected";

  string strId = "Id";
  vector<string> strMassPtEtaPhi; strMassPtEtaPhi.push_back("Pt"); strMassPtEtaPhi.push_back("Eta"); strMassPtEtaPhi.push_back("Phi"); strMassPtEtaPhi.push_back("Mass");
  vector<string> strMassPtPzPhi; strMassPtPzPhi.push_back("Pt"); strMassPtPzPhi.push_back("Pz"); strMassPtPzPhi.push_back("Phi"); strMassPtPzPhi.push_back("Mass");

  string strLepCore = "Lep";
  if (isGen) strLepCore.insert(0, "Gen");

  // Search the tree for the Higgs daughters
  vector<Particle*> candFinalDaughters;
  for (int v=0; v<2; v++){
    for (int d=0; d<2; d++){
      int iPart = 2*v+d+1;
      char cIPart[2];
      sprintf(cIPart, "%i", iPart);
      string strIPart = string(cIPart);

      bool success=true;
      Float_t* ref_partFV[4]={ 0 };
      // Get Lep1-4 PtEtaPhiM
      for (int fv=0; fv<4; fv++){
        varname = strLepCore + strIPart + strMassPtEtaPhi.at(fv);
        ref_partFV[fv] = (Float_t*)tree->getBranchHandleRef(varname);
        if (ref_partFV[fv]==0){ success=false; break; }
      }
      varname = strLepCore + strIPart + strId;
      Int_t* ref_partId = (Int_t*)tree->getBranchHandleRef(varname);
      if (ref_partId==0) success=false;
      if (!success) continue;

      Int_t partId = *ref_partId;
      TLorentzVector partFV; partFV.SetPtEtaPhiM(*(ref_partFV[0]), *(ref_partFV[1]), *(ref_partFV[2]), *(ref_partFV[3]));
      Particle* part = new Particle((int)partId, partFV);
      if(isGen) part->setGenStatus(1);
      particles.push_back(part);
      candFinalDaughters.push_back(part);
    }
  }

  // Search the tree for the associated particles
  vector<Particle*> associatedParticles;
  string strAPCore = "AssociatedParticle";
  if (isGen) strAPCore.insert(0, "Gen");
  vectorDouble** ref_apartFV[4];
  bool apart_success=true;
  for (int fv=0; fv<4; fv++){
    varname = strAPCore + strMassPtEtaPhi.at(fv);
    ref_apartFV[fv] = (vectorDouble**)tree->getBranchHandleRef(varname);
    if (ref_apartFV[fv]==0) { apart_success=false; break; }
    else if ((*ref_apartFV[fv])==0) { apart_success=false; break; }
  }
  varname = strAPCore + strId;
  vectorInt** ref_apartId = (vectorInt**)tree->getBranchHandleRef(varname);
  if (ref_apartId==0) apart_success=false;
  if ((*ref_apartId)==0) apart_success=false;
  if (apart_success){
    vectorInt* idhandle = *ref_apartId;
    vectorDouble* fvhandle[4];

    for (int fv=0; fv<4; fv++){
      fvhandle[fv]=*(ref_apartFV[fv]);
      if (fvhandle[fv]->size()!=idhandle->size()) { apart_success=false; break; }
    }
    for (int el=0; el<idhandle->size(); el++){
      Int_t partId = idhandle->at(el);
      TLorentzVector partFV; partFV.SetPtEtaPhiM(fvhandle[0]->at(el), fvhandle[1]->at(el), fvhandle[2]->at(el), fvhandle[3]->at(el));
      Particle* part = new Particle((int)partId, partFV);
      if (isGen) part->setGenStatus(1);
      particles.push_back(part);
      associatedParticles.push_back(part);
    }
  }

  // Reconstruct the gen and reco candidates
  if (isGen){ // Do not find the best combination, leave it to the input tree

    // First search the tree for the mothers
    vector<Particle*> motherParticles;
    string strMPCore = "Mother";
    strMPCore.insert(0, "Gen");
    vectorDouble** ref_mpartFV[4];
    bool mpart_success=true;
    for (int fv=0; fv<4; fv++){
      varname = strMPCore + strMassPtPzPhi.at(fv);
      ref_mpartFV[fv] = (vectorDouble**)tree->getBranchHandleRef(varname);
      if (ref_mpartFV[fv]==0) { mpart_success=false; break; }
      else if ((*ref_mpartFV[fv])==0) { mpart_success=false; break; }
    }
    varname = strMPCore + strId;
    vectorInt** ref_mpartId = (vectorInt**)tree->getBranchHandleRef(varname);
    if (ref_mpartId==0) mpart_success=false;
    if ((*ref_mpartId)==0) mpart_success=false;
    if (mpart_success){
      vectorInt* idhandle = *ref_mpartId;
      vectorDouble* fvhandle[4];

      for (int fv=0; fv<4; fv++){
        fvhandle[fv]=*(ref_mpartFV[fv]);
        if (fvhandle[fv]->size()!=idhandle->size()) { mpart_success=false; break; }
      }
      for (int el=0; el<idhandle->size(); el++){
        Int_t partId = idhandle->at(el);
        TLorentzVector partFV; partFV.SetXYZM(fvhandle[0]->at(el)*cos(fvhandle[2]->at(el)), fvhandle[0]->at(el)*sin(fvhandle[2]->at(el)), fvhandle[1]->at(el), fvhandle[3]->at(el));
        Particle* part = new Particle((int)partId, partFV);
        part->setGenStatus(-1);
        particles.push_back(part);
        motherParticles.push_back(part);
      }
    }

    for (int d=0; d<candFinalDaughters.size(); d++){
      Particle* part = candFinalDaughters.at(d);
      if (isALepton(part->id)) outEvent.addLepton(part);
      else if (isANeutrino(part->id)) outEvent.addNeutrino(part);
      else if (isAGluon(part->id) || isAQuark(part->id)) outEvent.addJet(part);
      else if (part->id==0) outEvent.addJet(part);
      else outEvent.addParticle(part);
    }
    outEvent.constructVVCandidates(options->doGenHZZdecay(), options->genDecayProducts());
    for (int p=0; p<motherParticles.size(); p++){
      Particle* part = motherParticles.at(p);
      outEvent.addVVCandidateMother(part);
    }
    for (int d=0; d<associatedParticles.size(); d++){
      Particle* part = associatedParticles.at(d);
      if (isALepton(part->id)) outEvent.addLepton(part);
      else if (isANeutrino(part->id)) outEvent.addNeutrino(part);
      else if (isAGluon(part->id) || isAQuark(part->id)) outEvent.addJet(part);
      else if (part->id==0) outEvent.addJet(part);
      else outEvent.addParticle(part);
    }
  }
  else{
    for (int d=0; d<candFinalDaughters.size(); d++){
      Particle* part = candFinalDaughters.at(d);
      if (isALepton(part->id)) outEvent.addLepton(part);
      else if (isANeutrino(part->id)) outEvent.addNeutrino(part);
      else if (isAGluon(part->id) || isAQuark(part->id)) outEvent.addJet(part);
      else if (part->id==0) outEvent.addJet(part);
      else outEvent.addParticle(part);
    }

    if (options->recoSelectionMode()==0){ // Do not find the best combination, leave it to the input tree
      outEvent.constructVVCandidates(options->doRecoHZZdecay(), options->recoDecayProducts());
      for (int d=0; d<associatedParticles.size(); d++){
        Particle* part = associatedParticles.at(d);
        if (isALepton(part->id)) outEvent.addLepton(part);
        else if (isANeutrino(part->id)) outEvent.addNeutrino(part);
        else if (isAGluon(part->id) || isAQuark(part->id)) outEvent.addJet(part);
        else if (part->id==0) outEvent.addJet(part);
        else outEvent.addParticle(part);
      }
    }
    else{ // Find the best combination
      for (int d=0; d<associatedParticles.size(); d++){
        Particle* part = associatedParticles.at(d);
        if (isALepton(part->id)) outEvent.addLepton(part);
        else if (isANeutrino(part->id)) outEvent.addNeutrino(part);
        else if (isAGluon(part->id) || isAQuark(part->id)) outEvent.addJet(part);
        else if (part->id==0) outEvent.addJet(part);
        else outEvent.addParticle(part);
      }
      outEvent.constructVVCandidates(options->doRecoHZZdecay(), options->recoDecayProducts());
      outEvent.applyParticleSelection();
    }
  }
  outEvent.addVVCandidateAppendages();

/*
  if (isGen) cout << "NGenCandidates: " << outEvent.getNZZCandidates() << endl;
  else cout << "NRecoCandidates: " << outEvent.getNZZCandidates() << endl;
  for (int cc=0; cc<outEvent.getNZZCandidates(); cc++){
  cout << outEvent.getZZCandidate(cc)->m()
  << " ("
  << outEvent.getZZCandidate(cc)->getSortedV(0)->m()
  << " -> "
  << "pT = " << outEvent.getZZCandidate(cc)->getSortedV(0)->getDaughter(0)->pt() << " + " << outEvent.getZZCandidate(cc)->getSortedV(0)->getDaughter(1)->pt()
  << ", "
  << outEvent.getZZCandidate(cc)->getSortedV(1)->m()
  << " -> "
  << "pT = " << outEvent.getZZCandidate(cc)->getSortedV(1)->getDaughter(0)->pt() << " + " << outEvent.getZZCandidate(cc)->getSortedV(1)->getDaughter(1)->pt()
  << ")\n";
  }
  if (isGen) cout << "Recorded (m1, m2): (" << *((Float_t*)tree->getBranchHandleRef("GenZ1Mass")) << ", " << *((Float_t*)tree->getBranchHandleRef("GenZ2Mass")) << ")\n";
  else cout << "Recorded (m1, m2): (" << *((Float_t*)tree->getBranchHandleRef("Z1Mass")) << ", " << *((Float_t*)tree->getBranchHandleRef("Z2Mass")) << ")\n";
  cout << endl;
*/
}

