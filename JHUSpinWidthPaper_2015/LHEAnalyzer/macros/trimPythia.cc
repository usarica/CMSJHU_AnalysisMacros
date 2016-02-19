//std includes
#include <iostream>
#include <utility>
#include <vector>
//root includes
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
//CMSSW includes
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "FWCore/Utilities/interface/TypeWithDict.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace std;

vector<pair<int, int>> findDuplicates(const vector<double>* fourvector, vector<int> id, vector<int> status);

void trimPythia(TString cinput, TString outdir="./", int fileLevel=1){
  TString coutput = outdir;
  coutput.Append("pythiaTemp.root");

  TFile* ftemp = new TFile(coutput, "recreate");
  TTree* tmpTree = new TTree("tmpTree", "");

  vector<double> geneventinfoweights;

  vector<double> reco_GenJet_FV[4];
  vector<int> reco_GenJet_id;
  vector<int> reco_GenJet_status;

  vector<double> reco_GenParticle_FV[4];
  vector<int> reco_GenParticle_id;
  vector<int> reco_GenParticle_status;

  tmpTree->Branch("genWeights", &geneventinfoweights);

  tmpTree->Branch("reco_GenJet_X", reco_GenJet_FV);
  tmpTree->Branch("reco_GenJet_Y", reco_GenJet_FV+1);
  tmpTree->Branch("reco_GenJet_Z", reco_GenJet_FV+2);
  tmpTree->Branch("reco_GenJet_E", reco_GenJet_FV+3);
  tmpTree->Branch("reco_GenJet_id", &reco_GenJet_id);
  tmpTree->Branch("reco_GenJet_status", &reco_GenJet_status);

  tmpTree->Branch("reco_GenParticle_X", reco_GenParticle_FV);
  tmpTree->Branch("reco_GenParticle_Y", reco_GenParticle_FV+1);
  tmpTree->Branch("reco_GenParticle_Z", reco_GenParticle_FV+2);
  tmpTree->Branch("reco_GenParticle_E", reco_GenParticle_FV+3);
  tmpTree->Branch("reco_GenParticle_id", &reco_GenParticle_id);
  tmpTree->Branch("reco_GenParticle_status", &reco_GenParticle_status);



  TFile f(cinput, "read");
  if (!f.IsZombie()){
    TTree* events = (TTree*)f.Get("Events");
    events->SetBranchStatus("*", 0);

    edm::Wrapper< vector<reco::GenJet> >* jetWrapper;
    edm::Wrapper< vector<reco::GenParticle> >* genparticleWrapper;
    edm::Wrapper< GenEventInfoProduct >* geneventinfoWrapper;

    TString suffix;
    if (fileLevel == 1) suffix = "SIM";
    else if (fileLevel == 2) suffix = "GEN";
    else {
      cout << "trimPythia should not be called with fileLevel=" << fileLevel << endl;
      assert(0);
    }
    events->SetBranchStatus("recoGenJets_ak5GenJets__"+suffix+"*", 1);
    events->SetBranchAddress("recoGenJets_ak5GenJets__"+suffix+".", &jetWrapper);
    events->SetBranchStatus("recoGenParticles_genParticles__"+suffix+"*", 1);
    events->SetBranchAddress("recoGenParticles_genParticles__"+suffix+".", &genparticleWrapper);
    events->SetBranchStatus("GenEventInfoProduct_generator__"+suffix+"*", 1);
    events->SetBranchAddress("GenEventInfoProduct_generator__"+suffix+".", &geneventinfoWrapper);

    for (int ev=0; ev<events->GetEntries(); ev++){
      events->GetEntry(ev);

      geneventinfoweights.clear();

      for (int v=0; v<4; v++){
        reco_GenJet_FV[v].clear();
        reco_GenParticle_FV[v].clear();
      }
      reco_GenJet_id.clear();
      reco_GenParticle_id.clear();
      reco_GenJet_status.clear();
      reco_GenParticle_status.clear();

      const gen::PdfInfo* genpdf = geneventinfoWrapper->product()->pdf();
      const vector<double> genweights = geneventinfoWrapper->product()->weights();
      for (unsigned int g=0; g<genweights.size(); g++) geneventinfoweights.push_back(genweights.at(g));

      const vector<reco::GenJet>* reco_GenJets = jetWrapper->product();
      const vector<reco::GenParticle>* reco_GenParticles = genparticleWrapper->product();

      for (unsigned int p=0; p<reco_GenParticles->size(); p++){
        const reco::GenParticle& part = reco_GenParticles->at(p);
        if (
          ((part.pdgId()==25 || part.pdgId()==32) && part.status()==22) // Generated Higgs
          ||
          ((std::abs(part.pdgId())>=11 && std::abs(part.pdgId())<=16) && (part.status()==23 || part.status()==1)) // Generated leptons
          ||
          ((std::abs(part.pdgId())<=6 || std::abs(part.pdgId())==21) && (part.status()==21 || part.status()==23 || part.status()==22)) // Generated partons
          ){
          reco_GenParticle_FV[0].push_back(part.px());
          reco_GenParticle_FV[1].push_back(part.py());
          reco_GenParticle_FV[2].push_back(part.pz());
          reco_GenParticle_FV[3].push_back(part.energy());
          reco_GenParticle_id.push_back(part.pdgId());
          reco_GenParticle_status.push_back(part.status());
        }
      }

      vector<pair<int, int>> reco_GenParticle_duplicates = findDuplicates(reco_GenParticle_FV, reco_GenParticle_id, reco_GenParticle_status);

      // Add status==1 or status==23 (with precedence for status==1) leptons to genJet collection
      for (unsigned int p=0; p<reco_GenParticle_id.size(); p++){
        if ((std::abs(reco_GenParticle_id.at(p))>=11 && std::abs(reco_GenParticle_id.at(p))<=16) && std::abs(reco_GenParticle_id.at(p)) % 2 == 1){
          bool match=false;
          for (unsigned int dd=0; dd<reco_GenParticle_duplicates.size(); dd++){
            pair<int, int> duplicate = reco_GenParticle_duplicates.at(dd);
            int iOriginal = duplicate.first; // Status==23 particle
            if (iOriginal==int(p)){ match=true; break; }
          }
          if (match) continue;
          for(int fv=0;fv<4;fv++) reco_GenJet_FV[fv].push_back(reco_GenParticle_FV[fv].at(p));
          reco_GenJet_id.push_back(reco_GenParticle_id.at(p));
          reco_GenJet_status.push_back(reco_GenParticle_status.at(p));
        }
      }

      vector<int> removalArray;
      for (unsigned int dd=0; dd<reco_GenParticle_duplicates.size(); dd++){
        pair<int, int> duplicate = reco_GenParticle_duplicates.at(dd);
        int iTransfer = duplicate.second; // Status==1 particle
        bool inserted=false;
        for (unsigned int it = 0; it<removalArray.size(); it++){
          int iIndex = removalArray.at(it);
          if (iTransfer > iIndex){
            removalArray.insert(removalArray.begin()+it, iTransfer);
            inserted=true;
            break;
          }
        }
        if (!inserted) removalArray.push_back(iTransfer);
      }
      // Remove status==1 duplicates from genParticles
      for (unsigned int dd=0; dd<removalArray.size(); dd++){
        int iTransfer = removalArray.at(dd);
        for (int fv=0; fv<4; fv++){
          reco_GenParticle_FV[fv].erase(reco_GenParticle_FV[fv].begin()+iTransfer);
        }
        reco_GenParticle_id.erase(reco_GenParticle_id.begin()+iTransfer);
        reco_GenParticle_status.erase(reco_GenParticle_status.begin()+iTransfer);
      }

      for (unsigned int p=0; p<reco_GenJets->size(); p++){
        const reco::GenJet& jet = reco_GenJets->at(p);
        reco_GenJet_FV[0].push_back(jet.px());
        reco_GenJet_FV[1].push_back(jet.py());
        reco_GenJet_FV[2].push_back(jet.pz());
        reco_GenJet_FV[3].push_back(jet.energy());
        reco_GenJet_id.push_back(jet.pdgId());
        reco_GenJet_status.push_back(jet.status());
      }

      tmpTree->Fill();
    }
    f.Close();
  }
  else if (f.IsOpen()) f.Close();
  ftemp->WriteTObject(tmpTree);
  delete tmpTree;
  ftemp->Close();
}

vector<pair<int, int>> findDuplicates(const vector<double>* fourvector, vector<int> id, vector<int> status){
  vector<pair<int, int>> duplicates;

  for (unsigned int xx=0; xx<id.size(); xx++){
    TLorentzVector p_tm(fourvector[0].at(xx), fourvector[1].at(xx), fourvector[2].at(xx), fourvector[3].at(xx));
    int id_tm = id.at(xx);
    int st_tm = status.at(xx);
    if (st_tm==23 || st_tm==1){
      for (unsigned int yy=xx+1; yy<id.size(); yy++){
        TLorentzVector p_tbm(fourvector[0].at(yy), fourvector[1].at(yy), fourvector[2].at(yy), fourvector[3].at(yy));
        int id_tbm = id.at(yy);
        int st_tbm = status.at(yy);

        if (id_tbm==id_tm && st_tbm!=st_tm && (st_tbm==1 || st_tbm==23)){
          double dot_tm = p_tbm.Dot(p_tm);
          double diff_sqmass = dot_tm - p_tm.M2();
          diff_sqmass = fabs((double)diff_sqmass);
          if (diff_sqmass<0.005*fabs(p_tm.M2())){
            int iFirst = (st_tm==23 ? xx : yy); // Should be order-independent
            int iSecond = (st_tbm==1 ? yy : xx);
            pair<int, int> dupinst(iFirst, iSecond);
            duplicates.push_back(dupinst);
          }
        }
      }
    }
  }

  return duplicates;
}

