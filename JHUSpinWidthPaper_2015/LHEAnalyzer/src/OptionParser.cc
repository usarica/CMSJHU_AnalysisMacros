#include "../interface/OptionParser.h"


OptionParser::OptionParser(int argc, char** argv) :
mPOLE(125.), // mH
wPOLE(4.07e-3), // GammaH
wPOLEStandard(4.07e-3), // GHSM
erg_tev(13), // C.o.M. energy in TeV
includeGenInfo(1), // Record gen. level quantities
includeRecoInfo(1), // Record reco. level quantities
removeDaughterMasses(1), // Force lepton masses to 0 in MELA
computeDecayAngles(1), // Decay angles
computeVBFAngles(0), // VBF production angles
computeVHAngles(0), // VH production angles
sampleProductionId(TVar::ZZGG, TVar::JHUGen), // Sample gen. production mode
fileLevel(0), // -1: ReadMode, 0: LHE, 1: Pythia,
isGenHZZ(1), // H->ZZ or H->WW
isRecoHZZ(1), // H->ZZ or H->WW
genDecayMode(0), // 4l with HZZ, 2l2nu with HWW, see Event::constructVVCandidates(bool isZZ, int fstype)
recoDecayMode(0), // 4l with HZZ, 2l2nu with HWW, see Event::constructVVCandidates(bool isZZ, int fstype)
recoSelBehaviour(0),
recoSmearBehaviour(0),
genHiggsCandidateSelectionScheme(HiggsComparators::BestZ1ThenZ2ScSumPt),
recoHiggsCandidateSelectionScheme(HiggsComparators::BestZ1ThenZ2ScSumPt),
indir("./"),
outdir("./"),
tmpDir("./tmpStore/"),
coutput("tmp.root"),
maxEvents(-1)
{
  for (int a=0; a<argc; a++){
    string tmpArg(argv[a]);
    rawOptions.push_back(tmpArg);
  }

  analyze();
}
void OptionParser::analyze(){
  bool hasInvalidOption=false;
  bool redefinedOutputFile=false;
  bool hasDecayAngles=false;
  bool hasGenProdProb=false;
  char rawdelimiter = '=';
  for (int opt=1; opt<rawOptions.size(); opt++){
    string wish, value;
    splitOption(rawOptions.at(opt), wish, value, rawdelimiter);
    interpretOption(wish, value);
    if (wish=="outfile") redefinedOutputFile=true;
    if (wish=="computeDecayAngles") hasDecayAngles=true;
    if (wish=="includeGenProdProb") hasGenProdProb=true;
  }

  if (filename.size()==0){ cerr << "You have to specify the input files." << endl; if(!hasInvalidOption) hasInvalidOption=true; }
  else{
    for (int f=0; f<filename.size(); f++){
      if ((filename.at(f).find(".lhe")!=string::npos && fileLevel!=0) || (filename.at(f).find(".root")!=string::npos && fileLevel==0)){
        cerr << "Inconsistent file name " << filename.at(f) << " and fileLevel option " << fileLevel << "!" << endl;
        if(!hasInvalidOption) hasInvalidOption=true;
      }
    }
  }

  if (maxEvents>=0){
    for (int es=0; es<eventSkipRanges.size(); es++) maxEvents += (eventSkipRanges.at(es).second-eventSkipRanges.at(es).first+1);
  }

  // Check for any invalid options
  if (isGenHZZ==-1){
    cout << "Gen. Higgs decay is disabled. Disabling MEs, decay angles and anything reco. as well." << endl;
    includeGenDecayProb.clear();
    includeRecoDecayProb.clear();
    includeGenProdProb.clear();
    includeRecoProdProb.clear();
    includeRecoInfo=0;
    hasDecayAngles=false; computeDecayAngles=0;
  }
  if (isRecoHZZ==-1){ cerr << "Reco. Higgs decay cannot be disabled." << endl; if (!hasInvalidOption) hasInvalidOption=true; }
  if (includeGenInfo==0 && includeRecoInfo==0){ cerr << "Cannot omit both reco. and gen. level info." << endl; if (!hasInvalidOption) hasInvalidOption=true; }
  if (mPOLE==0 || wPOLE==0 || erg_tev==0){ cerr << "Cannot have mH, GammaH or sqrts == 0" << endl; if(!hasInvalidOption) hasInvalidOption=true; }
  if (genHiggsCandidateSelectionScheme>=HiggsComparators::nCandidateSelections){ cerr << "Gen. H selection scheme is invalid!" << endl; if(!hasInvalidOption) hasInvalidOption=true; }
  if (recoHiggsCandidateSelectionScheme>=HiggsComparators::nCandidateSelections){ cerr << "Reco. H selection scheme is invalid!" << endl; if(!hasInvalidOption) hasInvalidOption=true; }
  if (hasGenProdProb && sampleProductionId.first==TVar::ZZGG){ cerr << "sampleProductionId==ZZGG is not a valid option (ME is not implemented). Use decay MEs instead for ZZGG or specify another production." << endl; if (!hasInvalidOption) hasInvalidOption=true; }

  if (!redefinedOutputFile) cout << "WARNING: No output file specified. Defaulting to " << coutput << "." << endl;
  if (!hasDecayAngles && fileLevel<0 && recoSelBehaviour==0) { cout << "Disabling the re-calculation of decay angles in ReadMode by default since no relevant option is specified." << endl; computeDecayAngles=0; }
  if (fileLevel<0 && recoSelBehaviour!=0) { cout << "Enabling the (re-)calculation of all angles in ReadMode since the re-selection option is specified." << endl; computeVBFAngles=1; computeVHAngles=1; computeDecayAngles=1; }

  // Print help if needed and abort at this point, nowhere later
  if (hasInvalidOption) printOptionsHelp();

  // Initialize the global Mela if needed
  configureMela();
}
void OptionParser::splitOption(string rawoption, string& wish, string& value, char delimiter){
  size_t posEq = rawoption.find(delimiter);
  if (posEq!=string::npos){
    wish=rawoption;
    value=rawoption.substr(posEq+1);
    wish.erase(wish.begin()+posEq,wish.end());
  }
  else{
    wish="";
    value=rawoption;
  }
}
void OptionParser::splitOptionRecursive(string rawoption, vector<string>& splitoptions, char delimiter){
  string suboption=rawoption, result=rawoption;
  string remnant;
  while (result!=""){
    splitOption(suboption, result, remnant, delimiter);
    if (result!="") splitoptions.push_back(result);
    suboption = remnant;
  }
  if (remnant!="") splitoptions.push_back(remnant);
}
Bool_t OptionParser::isAnExcludedBranch(string branchname){
  bool isExcluded=false;
  for (int eb=0; eb<excludedBranch.size(); eb++){
    if (branchname.find(excludedBranch.at(eb))!=string::npos && !(branchname.find("Gen")!=string::npos && excludedBranch.at(eb).find("Gen")==string::npos)){
      isExcluded=true;
      break;
    }
  }
  return isExcluded;
}
void OptionParser::extractSkippedEvents(string rawoption){
  vector<string> skipPair;
  splitOptionRecursive(rawoption, skipPair, ',');
  for (int p=0; p<skipPair.size(); p++){
    string strlow, strhigh;
    splitOption(skipPair.at(p), strlow, strhigh, '.');

    bool firstInclusive = true;
    size_t posFirstInc=strlow.find("[");
    size_t posFirstExc=strlow.find("(");
    if (posFirstInc!=string::npos || posFirstExc!=string::npos){
      if (posFirstInc!=string::npos && posFirstExc!=string::npos){
        cerr << "Invalid skipEvents range. Ignoring..." << endl;
        continue;
      }
      if (posFirstExc!=string::npos){
        firstInclusive=false;
        strlow = strlow.substr(posFirstExc+1);
      }
      else strlow = strlow.substr(posFirstInc+1);
    }
    Int_t firstId=(Int_t)atoi(strlow.c_str());
    if (!firstInclusive) firstId++;

    bool lastInclusive = false;
    size_t posLastInc=strhigh.find("]");
    size_t posLastExc=strhigh.find(")");
    if (posLastInc!=string::npos || posLastExc!=string::npos){
      if (posLastInc!=string::npos && posLastExc!=string::npos){
        cerr << "Invalid skipEvents range. Ignoring..." << endl;
        continue;
      }
      if (posLastExc!=string::npos) strhigh.erase(strhigh.begin()+posLastExc, strhigh.end());
      else{
        lastInclusive=true;
        strhigh.erase(strhigh.begin()+posLastInc, strhigh.end());
      }
    }
    Int_t lastId=(Int_t)atoi(strhigh.c_str());
    if (!lastInclusive) lastId--;

    if ((lastId>=0 && lastId<firstId) || (firstId<0 && lastId<0)){
      cerr << "Invalid skipEvents range. Ignoring..." << endl;
      continue;
    }
    pair<Int_t, Int_t> tmpPair(firstId, lastId);
    cout << "OptionParser: Will skip events " << tmpPair.first << " - " << tmpPair.second << endl;
    eventSkipRanges.push_back(tmpPair);
  }
}

void OptionParser::configureMela(){
  Int_t needMela = includeGenDecayProb.size()+includeRecoDecayProb.size()+includeGenProdProb.size()+includeRecoProdProb.size();
  if (needMela>0){
    melaHelpers::melaHandle = new Mela((int)erg_tev, (float)mPOLE);
  }
  melaHelpers::setSamplePoleWidth(wPOLE);
  melaHelpers::setStandardPoleWidth(wPOLEStandard);
  mela::applyLeptonMassCorrection(doRemoveLepMasses()); // Remains fixed, so nota problem to set it here
}
void OptionParser::deconfigureMela(){
  if (melaHelpers::melaHandle!=0) delete melaHelpers::melaHandle;
}
void OptionParser::extractMelaGenProdId(string rawoption){
  vector<string> prod_me_pair;
  splitOptionRecursive(rawoption, prod_me_pair, ',');
  if (prod_me_pair.size()!=2){
    cerr << "Incorrect specification for sampleProductionId. Has to follow the format (TVar::Production, TVar::MatrixElement) (no parentheses)." << endl;
    printOptionsHelp(); // The process ends here.
  }
  else{
    TVar::Production tmpProd; TVar::MatrixElement tmpME;

    if (prod_me_pair.at(0) == "JJGG") tmpProd = TVar::JJGG;
    else if (prod_me_pair.at(0) == "JJVBF") tmpProd = TVar::JJVBF;
    else if (prod_me_pair.at(0) == "JH") tmpProd = TVar::JH;
    else if (prod_me_pair.at(0) == "ZH") tmpProd = TVar::ZH;
    else if (prod_me_pair.at(0) == "WH") tmpProd = TVar::WH;
    else if (prod_me_pair.at(0) == "ttH") tmpProd = TVar::ttH;
    else if (prod_me_pair.at(0) == "bbH") tmpProd = TVar::bbH;
    else tmpProd = TVar::ZZGG;

    if (prod_me_pair.at(1) == "MCFM") tmpME = TVar::MCFM;
    else tmpME = TVar::JHUGen;

    pair<TVar::Production, TVar::MatrixElement> tmpPair(tmpProd, tmpME);
    sampleProductionId = tmpPair;
  }
}
Bool_t OptionParser::checkListVariable(vector<string>& list, string var){
  for (int v=0; v<list.size(); v++){
    if (list.at(v)==var) return true; // Look for exact match
  }
  return false;
}
Bool_t OptionParser::hasGenDecayME(string str){
  if (str=="" || str=="*"){
    return (includeGenDecayProb.size()>0 && processGenInfo());
  }
  return (checkListVariable(includeGenDecayProb, str) && processGenInfo());
}
Bool_t OptionParser::hasRecoDecayME(string str){
  if (str=="" || str=="*"){
    return (includeRecoDecayProb.size()>0 && processGenInfo());
  }
  return (checkListVariable(includeRecoDecayProb, str) && processRecoInfo());
}
Bool_t OptionParser::hasRecoProdME(string str){
  if (str=="" || str=="*"){
    return (includeRecoProdProb.size()>0 && processGenInfo());
  }
  return (checkListVariable(includeRecoProdProb, str) && processRecoInfo());
}
Bool_t OptionParser::hasGenProdME(string str){ // This one is a little bit trickier to avoid unneeded gen. prod. MEs
  if (str=="" || str=="*"){
    return (includeGenProdProb.size()>0 && processGenInfo());
  }
  return (checkListVariable(includeGenProdProb, str) && processGenInfo());
}

void OptionParser::interpretOption(string wish, string value){
  if (wish.empty()){
    if (value.find(".lhe")!=string::npos || value.find(".root")!=string::npos) filename.push_back(value);
    else if (value.find("help")!= string::npos) printOptionsHelp();
    else cerr << "Unknown unspecified argument: " << value << endl;
  }
  else if (wish=="includeGenInfo") includeGenInfo = (int)atoi(value.c_str());
  else if (wish=="includeRecoInfo") includeRecoInfo = (int)atoi(value.c_str());
  else if (wish=="fileLevel") fileLevel = (int)atoi(value.c_str());
  else if (wish=="isGenHZZ"){
    isGenHZZ = (int)atoi(value.c_str());
    if (isGenHZZ==0) PDGHelpers::setHVVmass(PDGHelpers::Wmass);
    else PDGHelpers::setHVVmass(PDGHelpers::Zmass);
  }
  else if (wish=="isRecoHZZ") isRecoHZZ = (int)atoi(value.c_str());
  else if (wish=="genDecayMode") genDecayMode = (int)atoi(value.c_str());
  else if (wish=="recoDecayMode") recoDecayMode = (int)atoi(value.c_str());
  else if (wish=="recoSelBehaviour" || wish=="recoSelBehavior") recoSelBehaviour = (int)atoi(value.c_str());
  else if (wish=="recoSmearBehaviour" || wish=="recoSmearBehavior") recoSmearBehaviour = (int)atoi(value.c_str());
  else if (wish=="genCandidateSelection" || wish=="genCandSel"){
    if (value=="BestZ1ThenZ2" || value=="BestZ1ThenZ2ScSumPt") genHiggsCandidateSelectionScheme = HiggsComparators::BestZ1ThenZ2ScSumPt;
  }
  else if (wish=="recoCandidateSelection" || wish=="recoCandSel"){
    if (value=="BestZ1ThenZ2" || value=="BestZ1ThenZ2ScSumPt") recoHiggsCandidateSelectionScheme = HiggsComparators::BestZ1ThenZ2ScSumPt;
  }
  else if (wish=="indir") indir = value;
  else if (wish=="outdir") outdir = value;
  else if (wish=="outfile") coutput = value;
  else if (wish=="tmpDir" || wish=="tempDir") tmpDir = value;
  else if (wish=="excludeBranch") splitOptionRecursive(value, excludedBranch, ',');
  else if (wish=="maxevents" || wish=="maxEvents") maxEvents = (int)atoi(value.c_str());
  else if (wish=="skipevents" || wish=="skipEvents") extractSkippedEvents(value);

  else if (wish=="mH" || wish=="MH" || wish=="mPOLE") mPOLE = (double)atof(value.c_str());
  else if (wish=="GH" || wish=="GaH" || wish=="GammaH" || wish=="wPOLE") wPOLE = (double)atof(value.c_str());
  else if (wish=="GHSM" || wish=="GaHSM" || wish=="GammaHSM" || wish=="wPOLEStandard") wPOLEStandard = (double)atof(value.c_str());
  else if (wish=="sqrts") erg_tev = (int)atoi(value.c_str());
  else if (wish=="removeDaughterMasses") removeDaughterMasses = (int)atoi(value.c_str());
  else if (wish=="computeDecayAngles") computeDecayAngles = (int)atoi(value.c_str());
  else if (wish=="computeVBFProdAngles") computeVBFAngles = (int)atoi(value.c_str());
  else if (wish=="computeVHProdAngles") computeVHAngles = (int)atoi(value.c_str());

  else if (wish=="includeRecoDecayProb") splitOptionRecursive(value, includeRecoDecayProb, ',');
  else if (wish=="includeRecoProdProb") splitOptionRecursive(value, includeRecoProdProb, ',');
  else if (wish=="includeGenDecayProb") splitOptionRecursive(value, includeGenDecayProb, ',');
  else if (wish=="includeGenProdProb") splitOptionRecursive(value, includeGenProdProb, ',');
  else if (wish=="sampleProductionId") extractMelaGenProdId(value);

  else cerr << "Unknown specified argument: " << value << " with specifier " << wish << endl;
}

void OptionParser::printOptionsHelp(){
  cout << endl;
  cout << "The options implemented for the LHEAnalyzer (format: specifier=value):\n\n";

  cout << "- No option specifier: Input files with extension .lhe or .root. Multiple input files can be passed as different arguments.\n\n";
  cout << "- indir: Location of input files. Default=\"./\"\n\n";
  cout << "- fileLevel: -1==ReadMode, 0==LHE, 1==Pythia8 (GEN-SIM), 2==Pythia8 (GEN). \".lhe\" extension only allowed for 0, and \".root\" is the only format for the others. Default=0\n\n";
  cout << "- outfile: Output file name. Default=\"tmp.root\"\n\n";
  cout << "- outdir: Location of the output file. Default=\"./\"\n\n";
  cout << "- tmpDir: Location of temporary files. Default=\"./tmpStore/\"\n\n";
  cout << "- maxevents / maxEvents: Maximum number of events to process. Default=-1 (all events)\n\n";
  cout << "- skipevents / skipEvents: Events to skip at the beginning. Default=none.\n\tAssignment is made in the form [ev1.ev2],[ev3.ev4] to allow multiple ranges. Use [ or ] for inclusive, ( or ) for exclusive ranges. Counting the events begins from 0.\n\n";

  cout << "- sqrts: pp collision c.o.m. energy. Default=13 (TeV)\n\n";
  cout << "- removeDaughterMasses: Switch to control the removal of lepton masses in the angle computation. Default=1\n\n";
  cout << "- computeDecayAngles: Switch to control the decay angles computation. Default=1 in LHE or Pythia modes, 0 in ReadMode.\n\n";
  cout << "- computeVBFProdAngles: Switch to control the VBF production angles computation. Default=0\n\n";
  cout << "- computeVHProdAngles: Switch to control the VH production angles computation. Default=0.\n\tPossible values are 0 (==disable), 1 (==compute from jets only), 2 (==compute from leptons only), 3 (==compute from jets and leptons, with separate variable suffixes \"*_VHhadronic\" and \"*_VHleptonic\".).\n\n";
  cout << "- includeGenDecayProb, includeGenProdProb, includeRecoDecayProb, includeRecoProdProb: Comma-separated list of spin-0 gen. or reco. decay or production MEs. Default=Empty (==None).\n\tThe explicit values tested are g1, g2, g4, g1_prime2; g1_pi2, g2_pi2, g4_pi2, g1_prime2_pi2, m4l, None, All.\n\t Gen. MEs are present for the purpose of reweighting. The appropriate target and origin combinations are left to the user at an analysis step.\n\n";
  cout << "- sampleProductionId: Production mechanism used in the computation of includeGenProdProb. Follows the format (TVar::Production, TVar::MatrixElement) with no parentheses. Default=(ZZGG, JHUGen)\n\tProduction==ZZGG is an invalid option since one can use decay MEs for this purpose, and the ghg_i couplings are not enabled in the MEs.\n\n";

  cout << "- mH / MH / mPOLE: Mass of the Higgs. Used in common for generator and reco. objects. Default=125 (GeV)\n\n";
  cout << "- GH / GaH / GammaH / wPOLE: Width of the generated Higgs. Used in generator objects. Default=4.07 (MeV)\n\n";
  cout << "- GHSM / GaHSM / GammaHSM / wPOLEStandard: Standard SM width. Used in scaling Mela probabilities properly. Default=4.07 (MeV).\n\n";
  cout << "- includeGenInfo, includeRecoInfo: Flags to control the writing of gen. and reco. info., respectively. Cannot be both false (0). Default=(1, 1)\n\n";
  cout << "- isGenHZZ, isRecoHZZ: Gen. or reco. H->VV decay. -1==H undecayed (gen.-only), 0==H->WW decay, 1==H->ZZ decay. isGenHZZ also (re)sets the default V mass in H->VV decay. Defaults=(1, 1)\n\n";
  cout << "- genDecayMode, recoDecayMode: Gen. or reco. H->VV->final states. Defaults=(0, 0)\n\tIf H->ZZ decay is specified, -1 - 5==Any, 4l, 4q, 2l2q, 2l2nu, 2q2nu, 4nu.\n\tIf H->WW decay is specified, -1 - 2==Any, 2l2nu, 4nu, lnu2q.\n\n";
  cout << "- recoSelBehavior / recoSelBehaviour: Selection behaviour on all reco. final states. Default=0.\n\t0==Apply selection in LHE and Pythia modes, apply no re-selection in ReadMode.\n\t1==Opposite of 0. Also enables the computation of all angles in ReadMode, overriding the relevant command line options.\n\n";
  cout << "- recoSmearBehavior / recoSmearBehaviour: Smearing behaviour on all reco. final states. Does not apply to ReadMode. Default=0.\n\t0==Apply smearing in LHE mode, no smearing in Pythia mode\n\t1==Opposite of 0\n\n";
  cout << "- genCandidateSelection, recoCandidateSelection: Higgs candidate selection algorithm. Values accepted are\n\t->BestZ1ThenZ2 (=BestZ1ThenZ2ScSumPt).\n\tDefaults==(BestZ1ThenZ2, BestZ1ThenZ2)\n\n";

  cout << "- excludeBranch: Comma-separated list of excluded branches. Default is to include all branches called via HVVTree::bookAllBranches.\n\n";


  cout << endl;
  assert(0);
}

void OptionParser::printOptionSummary(){

}

