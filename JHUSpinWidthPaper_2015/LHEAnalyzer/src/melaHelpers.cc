#include "../interface/melaHelpers.h"


namespace melaHelpers{
  Mela* melaHandle=0;
  Float_t genPOLEWidth=4.07e-3;
  Float_t standardPOLEWidth=4.07e-3;
}

using namespace PDGHelpers;

void melaHelpers::setSamplePoleWidth(Float_t width_){ genPOLEWidth = width_; }
void melaHelpers::setStandardPoleWidth(Float_t width_){ standardPOLEWidth = width_; }


Float_t melaHelpers::melaBranchMEInterpreter(const ZZCandidate* cand, string& branchname){
  Float_t result=0;

  // Initialize all variables that could be used
  TVar::MatrixElement myME = TVar::JHUGen;
  TVar::Production myProduction = TVar::ZZGG;
  TVar::Process myProcess = TVar::SelfDefine_spin0;
  TVar::SuperMelaSyst superMELASyst = TVar::SMSyst_None;
  TVar::EventScaleScheme rScaleScheme = TVar::DefaultScaleScheme;
  TVar::LeptonInterference leptonInterfScheme = TVar::DefaultLeptonInterf;
  bool hasSuperMELA=false;

  double selfDHggcoupl[SIZE_HGG][2] ={ { 0 } };
  double selfDHvvcoupl[SIZE_HVV_VBF][2] ={ { 0 } };
  double selfDHwwcoupl[SIZE_HWW_VBF][2] ={ { 0 } };
  double selfDTotalHvvcoupl[SIZE_HVV][2] ={ { 0 } };
  double selfDTotalHwwcoupl[SIZE_HVV][2] ={ { 0 } };
  double ggvvcoupl[2]={ 0, 0 };
  double mePoleScale = 1.;

  // Here we go...
  if (branchname.find("VAJHU")!=string::npos || branchname.find("m4l")!=string::npos) myME = TVar::JHUGen;
  else if (branchname.find("VAMCFM")!=string::npos) { myME = TVar::MCFM; myProcess = TVar::HSMHiggs; }

  if (branchname.find("vbf")!=string::npos) myProduction = TVar::JJVBF;
  else if (branchname.find("hjj")!=string::npos) myProduction = TVar::JJGG;
  else if (branchname.find("hj")!=string::npos){ myProduction = TVar::JH; myProcess = TVar::HSMHiggs; }
  else if (branchname.find("wh")!=string::npos) myProduction = TVar::WH;
  else if (branchname.find("zh")!=string::npos) myProduction = TVar::ZH;
  else if (branchname.find("tth")!=string::npos) myProduction = TVar::ttH;
  else if (branchname.find("bbh")!=string::npos) myProduction = TVar::bbH;
  else if (branchname.find("m4l")!=string::npos) { myProduction = TVar::ZZGG; hasSuperMELA=true; } // Important to quote it here before bkg_VAMCFM
  else if (branchname.find("bkg_s")!=string::npos) myProduction = TVar::ZZQQB_S;
  else if (branchname.find("bkg_tu")!=string::npos) myProduction = TVar::ZZQQB_TU;
  else if (branchname.find("bkg")!=string::npos) myProduction = TVar::ZZQQB;

  if (branchname.find("p0plus_m4l")!=string::npos) myProcess = TVar::HSMHiggs;
  else if (branchname.find("bkg")!=string::npos) myProcess = TVar::bkgZZ;
  else if (branchname.find("ggzz_VAMCFM")!=string::npos || branchname.find("VVzz_VAMCFM")!=string::npos) myProcess = TVar::bkgZZ;
  else if ((branchname.find("ggzz")!=string::npos || branchname.find("VVzz")!=string::npos) && branchname.find("VAMCFM")!=string::npos) myProcess = TVar::bkgZZ_SMHiggs;

  if (branchname.find("ScaleUp")!=string::npos) superMELASyst = TVar::SMSyst_ScaleUp;
  else if (branchname.find("ScaleDown")!=string::npos) superMELASyst = TVar::SMSyst_ScaleDown;
  else if (branchname.find("ResUp")!=string::npos) superMELASyst = TVar::SMSyst_ResUp;
  else if (branchname.find("ResDown")!=string::npos) superMELASyst = TVar::SMSyst_ResDown;
  melaHelpers::melaHandle->setProcess(myProcess, myME, myProduction);
  // ...phew!
  //cout << "melaHelpers::melaBranchMEInterpreter: " << branchname << " received. ";
  //cout << "Production: " << myProduction << ", " << "ME: " << myME << ", " << "Process: " << myProcess << ", " << "superMELASyst: " << superMELASyst << endl;

  // Set MELA flags
  if (branchname.find("Gen")!=string::npos){
    leptonInterfScheme = TVar::InterfOn;
    melaHelpers::melaHandle->setMelaHiggsWidth(melaHelpers::standardPOLEWidth);
    if (branchname.find("VAMCFM")!=string::npos && branchname.find("GHscaled")!=string::npos){
      melaHelpers::melaHandle->setMelaHiggsWidth(melaHelpers::genPOLEWidth);
      mePoleScale = melaHelpers::genPOLEWidth/melaHelpers::standardPOLEWidth; // pGen for MCFM needs to be scaled consistently.
    }
  }
  melaHelpers::melaHandle->setMelaLeptonInterference(leptonInterfScheme);

  TLorentzVector nullFourVector(0,0,0,0);
  TLorentzVector pHiggs = cand->p4;
  int Higgs_id = cand->id; // == 25

  vector<TLorentzVector> Higgs_daughters;
  vector<int> Higgs_daughter_ids;
  for (int v=0; v<2; v++){
    Particle* pV = cand->getSortedV(v);
    if (pV==0) continue;
    for (int d=0; d<2; d++){
      Particle* pVD = pV->getDaughter(d);
      if (pVD==0) continue;
      int iVD = 2*v+d;
      TLorentzVector mom(pVD->x(), pVD->y(), pVD->z(), pVD->t());
      Higgs_daughters.push_back(mom);
      Higgs_daughter_ids.push_back(pVD->id);
    }
  }

  if (Higgs_daughter_ids.size()<4 && (myME==TVar::MCFM || myProduction==TVar::ZZGG)) return result;

  if (hasSuperMELA){ // Exit route for SuperMELA
    TVar::LeptonFlavor superFlavor=TVar::Flavor_Dummy;
    if (abs(Higgs_daughter_ids.at(0))==abs(Higgs_daughter_ids.at(1)) &&
      abs(Higgs_daughter_ids.at(0))==abs(Higgs_daughter_ids.at(2)) &&
      abs(Higgs_daughter_ids.at(0))==abs(Higgs_daughter_ids.at(3))){
      if (abs(Higgs_daughter_ids.at(0))==11) superFlavor=TVar::Flavor_4e;
      else superFlavor=TVar::Flavor_4mu;
    }
    else superFlavor=TVar::Flavor_2e2mu;
    if (superFlavor!=TVar::Flavor_Dummy)
      melaHelpers::melaHandle->computePM4l(
      cand->m(),
      superFlavor,
      superMELASyst,
      result
      );
    return result;
  }

  // Find channel flavor
  int flavor=-1;
  if (
    abs(Higgs_daughter_ids.at(1))==abs(Higgs_daughter_ids.at(3))
    &&
    abs(Higgs_daughter_ids.at(0))==abs(Higgs_daughter_ids.at(2))
    ) flavor = 1;
  else flavor=3;

  // Special treatment for VH involves calling the routine with a single Higgs
  if ((myProduction==TVar::WH || myProduction==TVar::ZH) && myME==TVar::JHUGen){ // Avoid factor of 4
    Higgs_daughters.clear();
    Higgs_daughter_ids.clear();
    for (int d=0; d<3; d++) { Higgs_daughters.push_back(nullFourVector); Higgs_daughter_ids.push_back(0); }
    Higgs_daughters.push_back(pHiggs); Higgs_daughter_ids.push_back(Higgs_id);
  }

  vector<TLorentzVector> V_daughters;
  vector<int> V_daughter_ids;
  // Find the associated production leptons/jets
  if ((myProduction==TVar::WH || myProduction==TVar::ZH)){ // VH
    if (cand->getNSortedVs()<3) return result;
    else{
      for (int av=2; av<cand->getNSortedVs(); av++){ // Vs are already constructed after the jets, leptons etc. are ordered. Looking at their daughters therefore also orders the jets, leptons etc. by default.
        Particle* pV = cand->getSortedV(av);
        if (pV==0) continue;
        Particle* pVD1 = pV->getDaughter(0);
        Particle* pVD2 = pV->getDaughter(1);
        if (
          (
          (
          ((myProduction==TVar::WH && isAWBoson(pV->id)) || (myProduction==TVar::ZH && isAZBoson(pV->id))) && (isALepton(pVD1->id) || isANeutrino(pVD1->id))
          ) && branchname.find("leptonic")!=string::npos // Leptonic VH
          )
          ||
          (
          (
          (myProduction==TVar::WH && ((isAWBoson(pV->id) && isAQuark(pVD1->id)) || (pV->id==0 && pVD1->id==0)))
          ||
          (myProduction==TVar::ZH && ((isAZBoson(pV->id) && isAQuark(pVD1->id)) || (pV->id==0 && pVD1->id==0)))
          ) && branchname.find("hadronic")!=string::npos // Hadronic VH
          )
          ){
          TLorentzVector mom1(pVD1->x(), pVD1->y(), pVD1->z(), pVD1->t());
          V_daughters.push_back(mom1);
          V_daughter_ids.push_back(pVD1->id);
          TLorentzVector mom2(pVD2->x(), pVD2->y(), pVD2->z(), pVD2->t());
          V_daughters.push_back(mom2);
          V_daughter_ids.push_back(pVD2->id);
          break;
        }
      }
      if (V_daughter_ids.size()!=2) return result;
    }
  }
  else{ // Hadronic 2-jet MEs
    if (cand->getNAssociatedJets()>0){ // Uniform treatment for 1-jet cases
      if (
        (cand->getNAssociatedJets()<2 && myProduction==TVar::JJGG)
        ||
        (cand->getNAssociatedJets()>1 && myProduction==TVar::JH)
        ) return result;

      int nrealjets = min(2, cand->getNAssociatedJets());
      for (int av=0; av<nrealjets; av++){ // Jets are already ordered by pT.
        Particle* pJ = cand->getAssociatedJet(av);
        TLorentzVector mom(pJ->x(), pJ->y(), pJ->z(), pJ->t());
        V_daughters.push_back(mom);
        V_daughter_ids.push_back(pJ->id);
      }
      if (V_daughter_ids.size()==1){
        TLorentzVector mom(0, 0, 0, 0);
        V_daughters.push_back(mom);
        V_daughter_ids.push_back(0);
      }
    }
  }


  vector<string> gList[2];
  vector<pair<int, double> > gCoef;
  if (myProcess!=TVar::bkgZZ && !hasSuperMELA){
    if (myProduction==TVar::ttH || myProduction==TVar::bbH || myProduction==TVar::JH || myProduction==TVar::JJGG || (myProduction==TVar::ZZGG && branchname.find("prod")!=string::npos)){
      gList[0].push_back("g2"); gList[1].push_back("0plus");
      gList[0].push_back("g4"); gList[1].push_back("0minus");
      if (myProduction==TVar::ttH){ // ttH
        gCoef.push_back(pair<int, double>(0, 1));
        gCoef.push_back(pair<int, double>(2, 1.593));
      }
      else if (myProduction==TVar::bbH){ // bbH, the same as ttH for now
        gCoef.push_back(pair<int, double>(0, 1));
        gCoef.push_back(pair<int, double>(2, 1.593));
      }
      else if (myProduction==TVar::JH){ // 1-jet ggH, process == HSMHiggs, so selfDefine is not enabled.
        gCoef.push_back(pair<int, double>(0, 0));
        gCoef.push_back(pair<int, double>(2, 0));
      }
      else if (myProduction==TVar::JJGG){ // 2-jets ggH
        gCoef.push_back(pair<int, double>(0, 1));
        gCoef.push_back(pair<int, double>(2, sqrt(1.0017)));
      }
      else if (myProduction==TVar::ZZGG && branchname.find("prod")!=string::npos){ // 0-jet ggH prod. ME, the same as 2-jets ggH for now in absolute scale
        gCoef.push_back(pair<int, double>(0, 1));
        gCoef.push_back(pair<int, double>(2, sqrt(1.0017)));
      }
    }
    else{
      double coeffScale = 1.;
      if (myME==TVar::MCFM && branchname.find("Gen")!=string::npos && branchname.find("GHscaled")!=string::npos) coeffScale = pow(mePoleScale, 0.25);

      gList[0].push_back("g1"); gList[1].push_back("0plus"); gCoef.push_back(pair<int, double>(0, 1.*coeffScale));
      gList[0].push_back("g2"); gList[1].push_back("0hplus");
      gList[0].push_back("g4"); gList[1].push_back("0minus");
      gList[0].push_back("g1_prime2"); gList[1].push_back("0_g1prime2");
      if (myProduction==TVar::WH){ // WH
        gCoef.push_back(pair<int, double>(1, 0.0999064*coeffScale));
        gCoef.push_back(pair<int, double>(3, 0.123645*coeffScale));
        if (branchname.find("Gen")==string::npos){
          if (myME==TVar::MCFM) gCoef.push_back(pair<int, double>(11, -525.556*coeffScale));
          else gCoef.push_back(pair<int, double>(5, -525.556*coeffScale));
        }
        else{ // For a mostly positive discriminant
          if (myME==TVar::MCFM) gCoef.push_back(pair<int, double>(11, 525.556*coeffScale));
          else gCoef.push_back(pair<int, double>(5, 525.556*coeffScale));
        }
      }
      else if (myProduction==TVar::ZH){ // ZH
        gCoef.push_back(pair<int, double>(1, 0.112507*coeffScale));
        gCoef.push_back(pair<int, double>(3, 0.144082*coeffScale));
        if (branchname.find("Gen")==string::npos){
          if (myME==TVar::MCFM) gCoef.push_back(pair<int, double>(11, -517.761*coeffScale));
          else gCoef.push_back(pair<int, double>(5, -517.761*coeffScale));
        }
        else{ // For a mostly positive discriminant
          if (myME==TVar::MCFM) gCoef.push_back(pair<int, double>(11, 517.761*coeffScale));
          else gCoef.push_back(pair<int, double>(5, 517.761*coeffScale));
        }
      }
      else if (myProduction==TVar::JJVBF){ // VBF, to be revised
        gCoef.push_back(pair<int, double>(1, 0.270955*coeffScale));
        gCoef.push_back(pair<int, double>(3, 0.29724129*coeffScale));
        if (branchname.find("Gen")==string::npos){
          if (myME==TVar::MCFM) gCoef.push_back(pair<int, double>(11, -2132.143*coeffScale));
          else gCoef.push_back(pair<int, double>(5, -2132.143*coeffScale));
        }
        else{ // For a mostly positive discriminant
          if (myME==TVar::MCFM) gCoef.push_back(pair<int, double>(11, 2132.143*coeffScale));
          else gCoef.push_back(pair<int, double>(5, 2132.143*coeffScale));
        }
      }
      else if (myProduction==TVar::ZZGG){ // 0-jet ggH dec. ME
        gCoef.push_back(pair<int, double>(1, 1.65684*coeffScale));
        gCoef.push_back(pair<int, double>(3, 2.55052*coeffScale));
        if (branchname.find("Gen")==string::npos) gCoef.push_back(pair<int, double>(11, -12100.42*coeffScale));
        else gCoef.push_back(pair<int, double>(11, 12100.42*coeffScale)); // For a mostly positive discriminant
      }
    }

    int sgList = gList[0].size();
    bool** gFind;
    if (sgList>0){
      string trimmedBranchname = branchname;
      gFind = new bool*[sgList];
      for (int gg=0; gg<sgList; gg++){
        gFind[gg] = new bool[2];
        for (int im=0; im<2; im++) gFind[gg][im] = false;
      }
      for (int gg=0; gg<sgList; gg++){
        for (int im=1; im>=0; im--){
          for (int tt=0; tt<2; tt++){
            string chvar = gList[tt].at(gg);
            if (im==1) chvar.append("_pi2");
            size_t posFound = trimmedBranchname.find(chvar);
            if (posFound!=string::npos){
              gFind[gg][im]=true; // Does not care if it also finds g1g1_pi2, for example.
              trimmedBranchname.erase(trimmedBranchname.begin()+posFound, trimmedBranchname.begin()+posFound+chvar.length());
            }
          }
        }
      }
    }

    if (myProduction==TVar::ttH || myProduction==TVar::bbH || myProduction==TVar::JJGG || (myProduction==TVar::ZZGG && branchname.find("prod")!=string::npos)){
      for (int gg=0; gg<sgList; gg++){
        for (int im=0; im<2; im++){
          if (gFind[gg][im]) selfDHggcoupl[gCoef.at(gg).first][im] = gCoef.at(gg).second;
        }
      }
    }
    else if (myProduction==TVar::JH) {} // Do nothing
    else if (myProduction==TVar::JJVBF || myProduction==TVar::WH || myProduction==TVar::ZH || myProduction==TVar::ZZGG){
      for (int gg=0; gg<sgList; gg++){
        for (int im=0; im<2; im++){
          if (gFind[gg][im]){
            // Account for array size differences between MCFM and JHUGen associated production
            if (myME!=TVar::MCFM && myProduction!=TVar::ZZGG){
              selfDHvvcoupl[gCoef.at(gg).first][im] = gCoef.at(gg).second;
              selfDHwwcoupl[gCoef.at(gg).first][im] = gCoef.at(gg).second;
            }
            else{
              selfDTotalHvvcoupl[gCoef.at(gg).first][im] = gCoef.at(gg).second;
              selfDTotalHwwcoupl[gCoef.at(gg).first][im] = gCoef.at(gg).second;
            }
          }
        }
      }
    }
    // The arrays's task is complete.
    if (sgList>0){
      for (int gg=0; gg<sgList; gg++) delete[] gFind[gg];
      delete[] gFind;
    }
  }

  // Note: No implementation of tt/bbH yet!
  Float_t tmpME = 0;
  if ((myProduction==TVar::JJVBF || myProduction==TVar::JJGG || myProduction==TVar::JH) && myME==TVar::JHUGen){
    if (V_daughters.size()>1){
      Float_t auxME = 1;
      melaHelpers::melaHandle->computeProdP(
        V_daughters.at(0), V_daughter_ids.at(0),
        V_daughters.at(1), V_daughter_ids.at(1),
        pHiggs, Higgs_id,
        nullFourVector, 0,
        selfDHggcoupl,
        selfDHvvcoupl,
        selfDHwwcoupl,
        tmpME
        );
      melaHelpers::melaHandle->get_PAux(auxME);
      result = tmpME*auxME;
    }
  }
  else if ((myProduction==TVar::WH || myProduction==TVar::ZH) && myME==TVar::JHUGen){
    if (V_daughters.size()>1){
      // Unfortunately, cannot use vector::data
      TLorentzVector myjets[2] ={ V_daughters.at(0), V_daughters.at(1) };
      TLorentzVector daughters[4] ={ Higgs_daughters.at(0), Higgs_daughters.at(1), Higgs_daughters.at(2), Higgs_daughters.at(3) };
      int daughterids[4] ={ Higgs_daughter_ids.at(0), Higgs_daughter_ids.at(1), Higgs_daughter_ids.at(2), Higgs_daughter_ids.at(3) };
      int jetids[2] ={ V_daughter_ids.at(0), V_daughter_ids.at(1) };

      melaHelpers::melaHandle->computeProdP(
        myjets,
        daughters,
        jetids,
        daughterids,
        false,
        selfDHvvcoupl,
        tmpME);
      result = tmpME;
    }
  }
  else if ((myProduction==TVar::ZZGG || myProduction==TVar::ZZQQB) && !hasSuperMELA){
    Float_t helcosthetaZ1=0, helcosthetaZ2=0, helphi=0, costhetastar=0, phistarZ1=0;
    if (cand!=0) mela::computeAngles(
      cand->getSortedV(0)->getDaughter(0)->p4, cand->getSortedV(0)->getDaughter(0)->id,
      cand->getSortedV(0)->getDaughter(1)->p4, cand->getSortedV(0)->getDaughter(1)->id,
      cand->getSortedV(1)->getDaughter(0)->p4, cand->getSortedV(1)->getDaughter(0)->id,
      cand->getSortedV(1)->getDaughter(1)->p4, cand->getSortedV(1)->getDaughter(1)->id,
      costhetastar, helcosthetaZ1, helcosthetaZ2, helphi, phistarZ1
      );
    // Protect against NaN
    if (!(costhetastar==costhetastar)) costhetastar=0;
    if (!(helcosthetaZ1==helcosthetaZ1)) helcosthetaZ1=0;
    if (!(helcosthetaZ2==helcosthetaZ2)) helcosthetaZ2=0;
    if (!(helphi==helphi)) helphi=0;
    if (!(phistarZ1==phistarZ1)) phistarZ1=0;

    if (myProduction==TVar::ZZGG && myProcess!=TVar::bkgZZ){
      melaHelpers::melaHandle->computeP(
        (float)cand->m(), (float)cand->getSortedV(0)->m(), (float)cand->getSortedV(1)->m(),
        costhetastar, helcosthetaZ1, helcosthetaZ2, helphi, phistarZ1,
        flavor,
        selfDTotalHvvcoupl,
        tmpME
        );
    }
    else{ // qqZZ or ggZZ bkg.
      melaHelpers::melaHandle->computeP(
        (float)cand->m(), (float)cand->getSortedV(0)->m(), (float)cand->getSortedV(1)->m(),
        costhetastar, helcosthetaZ1, helcosthetaZ2, helphi, phistarZ1,
        flavor,
        tmpME,
        branchname.find("wconst")!=string::npos
        );
    }
    result = tmpME;
  }
  return result;
}


// Updated constrainedRemoveLeptonMass implementation
void melaHelpers::constrainedRemoveLeptonMass(TLorentzVector& p1, TLorentzVector& p2){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (p1==nullFourVector || p2==nullFourVector) return; // Need to add a smarter if-statement to take H->bb into account

  TLorentzVector pZ_old = p1 + p2;
  TLorentzVector p1_new = p1;
  TLorentzVector p2_new = p2;

  double energy = p1.T();
  double mom = p1.P();
  double mom_new = (energy + mom) / 2.;
  if (mom != 0){
    double x_new = p1.X()*mom_new / mom;
    double y_new = p1.Y()*mom_new / mom;
    double z_new = p1.Z()*mom_new / mom;
    p1_new.SetXYZT(x_new, y_new, z_new, mom_new);
  }
  energy = p2.T();
  mom = p2.P();
  mom_new = (energy + mom) / 2.;
  if (mom != 0){
    double x_new = p2.X()*mom_new / mom;
    double y_new = p2.Y()*mom_new / mom;
    double z_new = p2.Z()*mom_new / mom;
    p2_new.SetXYZT(x_new, y_new, z_new, mom_new);
  }
  TLorentzVector pZ_new = p1_new + p2_new;
  double delta_mZ=pZ_old.M()-pZ_new.M();

  p1_new.Boost(-(pZ_new.BoostVector()));
  mom = p1_new.T();
  mom_new = mom + delta_mZ / 2.;
  if (mom!=0) p1_new *= mom_new/mom;
  p1_new.Boost((pZ_old.BoostVector()));

  p2_new.Boost(-(pZ_new.BoostVector()));
  mom = p2_new.T();
  mom_new = mom + delta_mZ / 2.;
  if (mom!=0) p2_new *= mom_new/mom;
  p2_new.Boost((pZ_old.BoostVector()));

  p1=p1_new; p2=p2_new;
}

// Updated computeAngles implementation
void melaHelpers::computeAngles(
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1){
  /*
  cout << "dec begin" << endl;
  cout << "Entry:\n";
  cout << p4M11.X() << "\t" << p4M11.Y() << "\t" << p4M11.Z() << "\t" << p4M11.T() << endl;
  cout << p4M12.X() << "\t" << p4M12.Y() << "\t" << p4M12.Z() << "\t" << p4M12.T() << endl;
  cout << p4M21.X() << "\t" << p4M21.Y() << "\t" << p4M21.Z() << "\t" << p4M21.T() << endl;
  cout << p4M22.X() << "\t" << p4M22.Y() << "\t" << p4M22.Z() << "\t" << p4M22.T() << endl;
  */
  if (mela::forbidMassiveLeptons){
    if (!(fabs(Z1_lept1Id)==23 || fabs(Z1_lept1Id)==24 || fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==25 || fabs(Z1_lept2Id)==23 || fabs(Z1_lept2Id)==24 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==25)) melaHelpers::constrainedRemoveLeptonMass(p4M11, p4M12);
    if (!(fabs(Z2_lept1Id)==23 || fabs(Z2_lept1Id)==24 || fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==25 || fabs(Z2_lept2Id)==23 || fabs(Z2_lept2Id)==24 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==25)) melaHelpers::constrainedRemoveLeptonMass(p4M21, p4M22);
  }
  /*
  cout << "a.b.:\n";
  cout << p4M11.X() << "\t" << p4M11.Y() << "\t" << p4M11.Z() << "\t" << p4M11.T() << endl;
  cout << p4M12.X() << "\t" << p4M12.Y() << "\t" << p4M12.Z() << "\t" << p4M12.T() << endl;
  cout << p4M21.X() << "\t" << p4M21.Y() << "\t" << p4M21.Z() << "\t" << p4M21.T() << endl;
  cout << p4M22.X() << "\t" << p4M22.Y() << "\t" << p4M22.Z() << "\t" << p4M22.T() << endl;
  */

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;

  // Sort Z1 leptons so that:
  if ((Z1_lept1Id*Z1_lept2Id<0 && Z1_lept1Id<0) || // for OS pairs: lep1 must be the negative one
    (Z1_lept1Id*Z1_lept2Id>0 && p4M11.Phi()<=p4M12.Phi()) //for SS pairs: use random deterministic convention
    ) {
    swap(p4M11, p4M12);
  }

  // Same for Z2 leptons
  if ((Z2_lept1Id*Z2_lept2Id<0 && Z2_lept1Id<0) ||
    (Z2_lept1Id*Z2_lept2Id>0 && p4M21.Phi()<=p4M22.Phi())
    ) {
    swap(p4M21, p4M22);
  }

  // BEGIN THE CALCULATION

  // build H 4-vectors
  TLorentzVector p4H = p4Z1 + p4Z2;

  // -----------------------------------

  //// costhetastar
  TVector3 boostX = -(p4H.BoostVector());
  TLorentzVector thep4Z1inXFrame(p4Z1);
  TLorentzVector thep4Z2inXFrame(p4Z2);
  thep4Z1inXFrame.Boost(boostX);
  thep4Z2inXFrame.Boost(boostX);
  TVector3 theZ1X_p3 = TVector3(thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z());
  TVector3 theZ2X_p3 = TVector3(thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z());
  costhetastar = theZ1X_p3.CosTheta();

  //// --------------------------- costheta1
  TVector3 boostV1 = -(p4Z1.BoostVector());
  if (boostV1.Mag()>=1.) {
    cout << "Warning: Mela::computeAngles: Z1 boost with beta=1, scaling down" << endl;
    boostV1*=0.9999/boostV1.Mag();
  }
  TLorentzVector p4M11_BV1(p4M11);
  TLorentzVector p4M12_BV1(p4M12);
  TLorentzVector p4M21_BV1(p4M21);
  TLorentzVector p4M22_BV1(p4M22);
  p4M11_BV1.Boost(boostV1);
  p4M12_BV1.Boost(boostV1);
  p4M21_BV1.Boost(boostV1);
  p4M22_BV1.Boost(boostV1);

  TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
  //// costheta1
  costheta1 = -p4V2_BV1.Vect().Unit().Dot(p4M11_BV1.Vect().Unit());

  //// --------------------------- costheta2
  TVector3 boostV2 = -(p4Z2.BoostVector());
  if (boostV2.Mag()>=1.) {
    cout << "Warning: Mela::computeAngles: Z2 boost with beta=1, scaling down" << endl;
    boostV2*=0.9999/boostV2.Mag();
  }
  TLorentzVector p4M11_BV2(p4M11);
  TLorentzVector p4M12_BV2(p4M12);
  TLorentzVector p4M21_BV2(p4M21);
  TLorentzVector p4M22_BV2(p4M22);
  p4M11_BV2.Boost(boostV2);
  p4M12_BV2.Boost(boostV2);
  p4M21_BV2.Boost(boostV2);
  p4M22_BV2.Boost(boostV2);

  TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
  //// costheta2
  costheta2 = -p4V1_BV2.Vect().Unit().Dot(p4M21_BV2.Vect().Unit());

  //// --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
  TLorentzVector p4M11_BX(p4M11);
  TLorentzVector p4M12_BX(p4M12);
  TLorentzVector p4M21_BX(p4M21);
  TLorentzVector p4M22_BX(p4M22);

  p4M11_BX.Boost(boostX);
  p4M12_BX.Boost(boostX);
  p4M21_BX.Boost(boostX);
  p4M22_BX.Boost(boostX);
  TLorentzVector p4V1_BX = p4M11_BX + p4M12_BX;

  TVector3 beamAxis(0, 0, 1);
  TVector3 p3V1_BX = p4V1_BX.Vect().Unit();
  TVector3 normal1_BX = (p4M11_BX.Vect().Cross(p4M12_BX.Vect())).Unit();
  TVector3 normal2_BX = (p4M21_BX.Vect().Cross(p4M22_BX.Vect())).Unit();
  TVector3 normalSC_BX = (beamAxis.Cross(p3V1_BX)).Unit();


  //// Phi
  float tmpSgnPhi = p3V1_BX.Dot(normal1_BX.Cross(normal2_BX));
  float sgnPhi = 0;
  if (fabs(tmpSgnPhi)>0) sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  float dot_BX12 = normal1_BX.Dot(normal2_BX);
  if (fabs(dot_BX12)>=1) dot_BX12 *= 1./fabs(dot_BX12);
  Phi = sgnPhi * acos(-1.*dot_BX12);


  //// Phi1
  float tmpSgnPhi1 = p3V1_BX.Dot(normal1_BX.Cross(normalSC_BX));
  float sgnPhi1 = 0;
  if (fabs(tmpSgnPhi1)>0) sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);
  float dot_BX1SC = normal1_BX.Dot(normalSC_BX);
  if (fabs(dot_BX1SC)>=1) dot_BX1SC *= 1./fabs(dot_BX1SC);
  Phi1 = sgnPhi1 * acos(dot_BX1SC);

  if (isnan(costhetastar) || isnan(costheta1) || isnan(costheta2) || isnan(Phi) || isnan(Phi1)){
    cout << "WARNING: NaN in computeAngles: "
      << costhetastar << " "
      << costheta1  << " "
      << costheta2  << " "
      << Phi  << " "
      << Phi1  << " " << endl;
    cout << "   boostV1: " <<boostV1.Pt() << " " << boostV1.Eta() << " " << boostV1.Phi() << " " << boostV1.Mag() << endl;
    cout << "   boostV2: " <<boostV2.Pt() << " " << boostV2.Eta() << " " << boostV2.Phi() << " " << boostV2.Mag() << endl;
  }
}


void melaHelpers::computeVBFangles(
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1,
  float& Q2V1,
  float& Q2V2,
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  TLorentzVector jet1, int jet1Id,
  TLorentzVector jet2, int jet2Id,
  TLorentzVector* injet1, int injet1Id, // Gen. partons in lab frame
  TLorentzVector* injet2, int injet2Id
  ){

  if (mela::forbidMassiveLeptons){
    if (!(fabs(Z1_lept1Id)==23 || fabs(Z1_lept1Id)==24 || fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==25 || fabs(Z1_lept2Id)==23 || fabs(Z1_lept2Id)==24 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==25)) melaHelpers::constrainedRemoveLeptonMass(p4M11, p4M12);
    if (!(fabs(Z2_lept1Id)==23 || fabs(Z2_lept1Id)==24 || fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==25 || fabs(Z2_lept2Id)==23 || fabs(Z2_lept2Id)==24 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==25)) melaHelpers::constrainedRemoveLeptonMass(p4M21, p4M22);
  }

  TLorentzVector jet1massless, jet2massless;
  mela::computeJetMassless(jet1, jet1massless);
  mela::computeJetMassless(jet2, jet2massless);

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;
  TLorentzVector pH = p4Z1+p4Z2;
  //jet1 is defined as going forwards (bigger pz), jet2 going backwards (smaller pz)
  if (jet1massless.Z() < jet2massless.Z()) { swap(jet1massless, jet2massless); swap(jet1Id, jet2Id); }

  //Find the incoming partons - first boost so the pT(HJJ) = 0, then boost away the pz.
  //This preserves the z direction.  Then assume that the partons come in this z direction.
  //This is exactly correct at LO (since pT=0 anyway).
  //Then associate the one going forwards with jet1 and the one going backwards with jet2
  TLorentzRotation movingframe;
  TLorentzVector pHJJ = pH+jet1massless+jet2massless;
  TLorentzVector pHJJ_perp(pHJJ.X(), pHJJ.Y(), 0, pHJJ.T());
  movingframe.Boost(-pHJJ_perp.BoostVector());
  pHJJ.Boost(-pHJJ_perp.BoostVector());
  movingframe.Boost(-pHJJ.BoostVector());
  pHJJ.Boost(-pHJJ.BoostVector());   //make sure to boost HJJ AFTER boosting movingframe

  TLorentzVector P1(0, 0, pHJJ.T()/2, pHJJ.T()/2);
  TLorentzVector P2(0, 0, -pHJJ.T()/2, pHJJ.T()/2);
  // Transform incoming partons back to original frame
  P1.Transform(movingframe.Inverse());
  P2.Transform(movingframe.Inverse());
  //movingframe, HJJ, and HJJ_T will not be used anymore
  if (injet1!=0 && injet2!=0){ // Handle gen. partons if they are available
    if (fabs((*injet1+*injet2).P()-pHJJ.P())<pHJJ.P()*1e-4){
      P1=*injet1;
      P2=*injet2;
      if (P1.Z() < P2.Z()){
        swap(P1, P2);
        swap(injet1Id, injet2Id);
      }
      // In the case of gen. partons, check if the intermediates are a Z or a W.
      int diff1Id = jet1Id-injet1Id;
      int diff2Id = jet2Id-injet2Id;
      if (
        !(
        (diff1Id==0 && diff2Id==0 && !(injet1Id==21 || injet2Id==21)) // Two Z bosons
        ||
        ((fabs(diff1Id)==1 || fabs(diff1Id)==3 || fabs(diff1Id)==5) && (fabs(diff2Id)==1 || fabs(diff2Id)==3 || fabs(diff2Id)==5)) // Two W bosons, do not check W+ vs W-
        )
        ){
        int diff12Id = jet1Id-injet2Id;
        int diff21Id = jet2Id-injet1Id;
        if (
          ((diff12Id==0 || diff21Id==0) && !(injet1Id==21 || injet2Id==21)) // At least one Z boson
          ||
          ((fabs(diff12Id)==1 || fabs(diff12Id)==3 || fabs(diff12Id)==5) || (fabs(diff21Id)==1 || fabs(diff21Id)==3 || fabs(diff21Id)==5)) // At least one W boson
          ){
          swap(P1, P2);
          swap(injet1Id, injet2Id);
        }
      }
    }
  }

  TLorentzRotation ZZframe;
  ZZframe.Boost(-pH.BoostVector());
  P1.Transform(ZZframe);
  P2.Transform(ZZframe);
  p4Z1.Transform(ZZframe);
  p4Z2.Transform(ZZframe);
  jet1massless.Transform(ZZframe);
  jet2massless.Transform(ZZframe);

  TLorentzVector V1 = P1-jet1massless; // V1 = (-p12) - p11 = -Z1
  TLorentzVector V2 = P2-jet2massless; // V2 = (-p22) - p21 = -Z2
  Q2V1 = -V1.M2();
  Q2V2 = -V2.M2();

  costhetastar = -V1.Vect().Unit().Dot(p4Z2.Vect().Unit());
  costheta1 = -V1.Vect().Unit().Dot(jet1massless.Vect().Unit());
  costheta2 = -V2.Vect().Unit().Dot(jet2massless.Vect().Unit());

  TVector3 normvec1 = P1.Vect().Cross(jet1massless.Vect()).Unit(); // p11 x p12 = (-p12) x p11
  TVector3 normvec2 = P2.Vect().Cross(jet2massless.Vect()).Unit(); // p21 x p22 = (-p22) x p21
  TVector3 normvec3 = V1.Vect().Cross(p4Z2.Vect()).Unit(); // == z x Z1

  double cosPhi = normvec1.Dot(normvec2);
  double sgnPhi = normvec1.Cross(normvec2).Dot(-V1.Vect());
  double cosPhi1 = normvec1.Dot(normvec3);
  double sgnPhi1 = normvec1.Cross(normvec3).Dot(-V1.Vect());
  if (fabs(cosPhi)>1) cosPhi *= 1./fabs(cosPhi);
  if (fabs(cosPhi1)>1) cosPhi1 *= 1./fabs(cosPhi1);
  Phi = TMath::Sign(acos(-cosPhi), sgnPhi);            //TMath::Sign(a,b) = |a|*(b/|b|)
  Phi1 = TMath::Sign(acos(cosPhi1), sgnPhi1);
}


void melaHelpers::computeVHangles(
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1,
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  TLorentzVector jet1, int jet1Id,
  TLorentzVector jet2, int jet2Id,
  TLorentzVector* injet1, int injet1Id, // Gen. partons in lab frame
  TLorentzVector* injet2, int injet2Id
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);

  //cout << p4M11.X() << '\t' << p4M11.Y() << '\t' << p4M11.Z() << '\t' << p4M11.T() << endl;
  //cout << p4M12.X() << '\t' << p4M12.Y() << '\t' << p4M12.Z() << '\t' << p4M12.T() << endl;
  //cout << p4M21.X() << '\t' << p4M21.Y() << '\t' << p4M21.Z() << '\t' << p4M21.T() << endl;
  //cout << p4M22.X() << '\t' << p4M22.Y() << '\t' << p4M22.Z() << '\t' << p4M22.T() << endl;
  if (mela::forbidMassiveLeptons){
    if (!(fabs(Z1_lept1Id)==23 || fabs(Z1_lept1Id)==24 || fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==25 || fabs(Z1_lept2Id)==23 || fabs(Z1_lept2Id)==24 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==25)) melaHelpers::constrainedRemoveLeptonMass(p4M11, p4M12);
    if (!(fabs(Z2_lept1Id)==23 || fabs(Z2_lept1Id)==24 || fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==25 || fabs(Z2_lept2Id)==23 || fabs(Z2_lept2Id)==24 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==25)) melaHelpers::constrainedRemoveLeptonMass(p4M21, p4M22);
  }
  //cout << endl;
  //cout << p4M11.X() << '\t' << p4M11.Y() << '\t' << p4M11.Z() << '\t' << p4M11.T() << endl;
  //cout << p4M12.X() << '\t' << p4M12.Y() << '\t' << p4M12.Z() << '\t' << p4M12.T() << endl;
  //cout << p4M21.X() << '\t' << p4M21.Y() << '\t' << p4M21.Z() << '\t' << p4M21.T() << endl;
  //cout << p4M22.X() << '\t' << p4M22.Y() << '\t' << p4M22.Z() << '\t' << p4M22.T() << endl;

  // Build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;
  TLorentzVector pH = p4Z1 + p4Z2;

  TLorentzVector jet1massless, jet2massless;
  mela::computeJetMassless(jet1, jet1massless);
  mela::computeJetMassless(jet2, jet2massless);

  // Apply convention for outgoing particles
  if (
    (jet1Id*jet2Id<0 && jet1Id<0) // for OS pairs: jet1 must be the particle
    ||
    (jet1Id*jet2Id>0 && jet1massless.Phi()<=jet2massless.Phi()) // for SS pairs: use random deterministic convention
    ){
    swap(jet1massless, jet2massless);
    swap(jet1Id, jet2Id);
  }

  //Find the incoming partons - first boost so the pT(HJJ) = 0, then boost away the pz.
  //This preserves the z direction.  Then assume that the partons come in this z direction.
  //This is exactly correct at LO (since pT=0 anyway).
  //Then associate the one going forwards with jet1 and the one going backwards with jet2
  TLorentzRotation movingframe;
  TLorentzVector pHJJ = pH+jet1massless+jet2massless;
  TLorentzVector pHJJ_perp(pHJJ.X(), pHJJ.Y(), 0, pHJJ.T());
  movingframe.Boost(-pHJJ_perp.BoostVector());
  pHJJ.Boost(-pHJJ_perp.BoostVector());
  movingframe.Boost(-pHJJ.BoostVector());
  pHJJ.Boost(-pHJJ.BoostVector());   //make sure to boost HJJ AFTER boosting movingframe

  TLorentzVector P1(0, 0, -pHJJ.T()/2, pHJJ.T()/2);
  TLorentzVector P2(0, 0, pHJJ.T()/2, pHJJ.T()/2);
  // Transform incoming partons back to the original frame
  P1.Transform(movingframe.Inverse());
  P2.Transform(movingframe.Inverse());
  //movingframe, HJJ, and HJJ_T will not be used anymore
  if (injet1!=0 && injet2!=0){ // Handle gen. partons if they are available
    if (fabs((*injet1+*injet2).P()-pHJJ.P())<=pHJJ.P()*1e-4){
      P1=*injet1;
      P2=*injet2;
      // Apply convention for incoming (!) particles
      if (
        (injet1Id*injet2Id<0 && injet1Id>0) // for OS pairs: parton 2 must be the particle
        ||
        (injet1Id*injet2Id>0 && P1.Z()>=P2.Z()) //for SS pairs: use random deterministic convention
        ){
        swap(P1, P2);
        swap(injet1Id, injet2Id);
      }
    }
  }

  // Rotate every vector such that Z1 - Z2 axis is the "beam axis" analogue of decay
  TLorentzRotation ZZframe;
  TVector3 beamAxis(0, 0, 1);
  if (p4Z1==nullFourVector || p4Z2==nullFourVector){
    TVector3 pNewAxis = (p4Z2-p4Z1).Vect().Unit(); // Let Z2 be in the z direction so that once the direction of H is reversed, Z1 is in the z direction
    if (pNewAxis != nullFourVector.Vect()){
      TVector3 pNewAxisPerp = pNewAxis.Cross(beamAxis);
      ZZframe.Rotate(acos(pNewAxis.Dot(beamAxis)), pNewAxisPerp);

      P1.Transform(ZZframe);
      P2.Transform(ZZframe);
      jet1massless = -jet1massless;
      jet2massless = -jet2massless;
      jet1massless.Transform(ZZframe);
      jet2massless.Transform(ZZframe);
      jet1massless = -jet1massless;
      jet2massless = -jet2massless;
    }
  }
  else{
    ZZframe.Boost(-pH.BoostVector());
    p4Z1.Boost(-pH.BoostVector());
    p4Z2.Boost(-pH.BoostVector());
    TVector3 pNewAxis = (p4Z2-p4Z1).Vect().Unit(); // Let Z2 be in the z direction so that once the direction of H is reversed, Z1 is in the z direction
    TVector3 pNewAxisPerp = pNewAxis.Cross(beamAxis);
    ZZframe.Rotate(acos(pNewAxis.Dot(beamAxis)), pNewAxisPerp);
    P1.Transform(ZZframe);
    P2.Transform(ZZframe);
    jet1massless = -jet1massless;
    jet2massless = -jet2massless;
    jet1massless.Transform(ZZframe);
    jet2massless.Transform(ZZframe);
    jet1massless = -jet1massless;
    jet2massless = -jet2massless;
  }

  melaHelpers::computeAngles(
    -P1, 23, // Id is 23 to avoid an attempt to remove quark mass
    -P2, 0, // Id is 0 to avoid swapping
    jet1massless, 23,
    jet2massless, 0,
    costhetastar,
    costheta1,
    costheta2,
    Phi,
    Phi1);
}


