#include "../interface/BaseTree.h"

BaseTree::BaseTree(string treename){ initTree(treename, ""); }
BaseTree::BaseTree(string treename, string treetitle){ initTree(treename, treetitle); }
BaseTree::BaseTree(string treename, TFile* fin){ getTreeFromFile(treename, fin); }

void BaseTree::getTreeFromFile(string treename, TFile* fin){
  if (fin!=0 && !fin->IsZombie() && fin->IsOpen()) hvvtree = (TTree*)fin->Get(treename.c_str());
  else hvvtree=0;
  if (hvvtree==0) cout << "Failed to extract the tree named " << treename << "!" << endl;
}

bool BaseTree::bookBranch(string branchname, BaseTree::BranchTypes bType, bool doSetAddress){
  bool success=true;
  if (hvvtree!=0){
    if (bType==BaseTree::bInt){
      Int_t* container = 0;
      if (!doSetAddress || hvvtree->GetBranchStatus(branchname.c_str())){
        container = new Int_t; // Cannot be a null pointer even if setting branch address
        pair<string, Int_t*> varPair(branchname, container);
        intBranches.push_back(varPair);
      }
      else success=false;
    }
    else if (bType==BaseTree::bFloat){
      Float_t* container = 0;
      if (!doSetAddress || hvvtree->GetBranchStatus(branchname.c_str())){
        container = new Float_t; // Cannot be a null pointer even if setting branch address
        pair<string, Float_t*> varPair(branchname, container);
        floatBranches.push_back(varPair);
      }
      else success=false;
    }
    else if (bType==BaseTree::bVectorInt){
      vectorInt* container = 0;
      if (!doSetAddress) container = new vectorInt;
      if (!doSetAddress || hvvtree->GetBranchStatus(branchname.c_str())){
        pair<string, vectorInt*> varPair(branchname, container);
        vectorIntBranches.push_back(varPair);
      }
      else success=false;
    }
    else if (bType==BaseTree::bVectorDouble){
      vectorDouble* container = 0;
      if (!doSetAddress) container = new vectorDouble;
      if (!doSetAddress || hvvtree->GetBranchStatus(branchname.c_str())){
        pair<string, vectorDouble*> varPair(branchname, container);
        vectorDoubleBranches.push_back(varPair);
      }
      else success=false;
    }
    success=true;
  }
  else success=false;
  if (!success) cerr << "BaseTree::bookBranch: Failed to book branch named " << branchname << "!" << endl;
  return success;
}
bool BaseTree::actuateBranches(bool doSetAddress){
  bool success=true;
  if (hvvtree!=0){
    if (!doSetAddress){
      for (int el=0; el<intBranches.size(); el++){
        if (!hvvtree->GetBranchStatus(intBranches.at(el).first.c_str()))
          hvvtree->Branch(intBranches.at(el).first.c_str(), intBranches.at(el).second);
      }
      for (int el=0; el<floatBranches.size(); el++){
        if (!hvvtree->GetBranchStatus(floatBranches.at(el).first.c_str()))
          hvvtree->Branch(floatBranches.at(el).first.c_str(), floatBranches.at(el).second);
      }
      for (int el=0; el<vectorIntBranches.size(); el++){
        if (!hvvtree->GetBranchStatus(vectorIntBranches.at(el).first.c_str()))
          hvvtree->Branch(vectorIntBranches.at(el).first.c_str(), vectorIntBranches.at(el).second);
      }
      for (int el=0; el<vectorDoubleBranches.size(); el++){
        if (!hvvtree->GetBranchStatus(vectorDoubleBranches.at(el).first.c_str()))
          hvvtree->Branch(vectorDoubleBranches.at(el).first.c_str(), vectorDoubleBranches.at(el).second);
      }
    }
    else{
      for (int el=0; el<intBranches.size(); el++) hvvtree->SetBranchAddress(intBranches.at(el).first.c_str(), intBranches.at(el).second); // Already a pointer
      for (int el=0; el<floatBranches.size(); el++) hvvtree->SetBranchAddress(floatBranches.at(el).first.c_str(), floatBranches.at(el).second); // Already a pointer
      for (int el=0; el<vectorIntBranches.size(); el++) hvvtree->SetBranchAddress(vectorIntBranches.at(el).first.c_str(), &(vectorIntBranches.at(el).second)); // Need to pass the address of the pointer to std::vector
      for (int el=0; el<vectorDoubleBranches.size(); el++) hvvtree->SetBranchAddress(vectorDoubleBranches.at(el).first.c_str(), &(vectorDoubleBranches.at(el).second)); // Need to pass the address of the pointer to std::vector
    }
  }
  else success=false;
  if (!success) cerr << "BaseTree::actuateBranch: Failed to actuate the branches!" << endl;
  return success;
}
vector<string> BaseTree::getBranchList(){
  vector<string> branchlist;
  for (int el=0; el<intBranches.size(); el++) branchlist.push_back(intBranches.at(el).first);
  for (int el=0; el<floatBranches.size(); el++) branchlist.push_back(floatBranches.at(el).first);
  for (int el=0; el<vectorIntBranches.size(); el++) branchlist.push_back(vectorIntBranches.at(el).first);
  for (int el=0; el<vectorDoubleBranches.size(); el++) branchlist.push_back(vectorDoubleBranches.at(el).first);
  return branchlist;
}

BaseTree::BranchTypes BaseTree::searchArray(string branchname, int& position){
  for (int el=0; el<intBranches.size(); el++){
    if (branchname==intBranches.at(el).first){
      position = el;
      return bInt;
    }
  }
  for (int el=0; el<floatBranches.size(); el++){
    if (branchname==floatBranches.at(el).first){
      position = el;
      return bFloat;
    }
  }
  for (int el=0; el<vectorIntBranches.size(); el++){
    if (branchname==vectorIntBranches.at(el).first){
      position = el;
      return bVectorInt;
    }
  }
  for (int el=0; el<vectorDoubleBranches.size(); el++){
    if (branchname==vectorDoubleBranches.at(el).first){
      position = el;
      return bVectorDouble;
    }
  }
  return BaseTree::nBranchTypes;
}
void BaseTree::cleanBranches(){
  for (int el=0; el<intBranches.size(); el++){
    if (intBranches.at(el).second!=0) delete intBranches.at(el).second;
    intBranches.at(el).second=0;
  }
  intBranches.clear();
  for (int el=0; el<floatBranches.size(); el++){
    if (floatBranches.at(el).second!=0) delete floatBranches.at(el).second;
    floatBranches.at(el).second=0;
  }
  floatBranches.clear();
  for (int el=0; el<vectorIntBranches.size(); el++){
    if (vectorIntBranches.at(el).second!=0){
      vectorIntBranches.at(el).second->clear();
      delete vectorIntBranches.at(el).second;
    }
    vectorIntBranches.at(el).second=0;
  }
  vectorIntBranches.clear();
  for (int el=0; el<vectorDoubleBranches.size(); el++){
    if (vectorDoubleBranches.at(el).second!=0){
      vectorDoubleBranches.at(el).second->clear();
      delete vectorDoubleBranches.at(el).second;
    }
    vectorDoubleBranches.at(el).second=0;
  }
  vectorDoubleBranches.clear();
}
void BaseTree::initializeBranches(){
  for (int el=0; el<intBranches.size(); el++){
    if (intBranches.at(el).second!=0) *(intBranches.at(el).second)=0;
  }
  for (int el=0; el<floatBranches.size(); el++){
    if (floatBranches.at(el).second!=0) *(floatBranches.at(el).second)=0;
  }
  for (int el=0; el<vectorIntBranches.size(); el++){
    if (vectorIntBranches.at(el).second!=0) vectorIntBranches.at(el).second->clear();
  }
  for (int el=0; el<vectorDoubleBranches.size(); el++){
    if (vectorDoubleBranches.at(el).second!=0) vectorDoubleBranches.at(el).second->clear();
  }
}
void BaseTree::printEntry(int ev){
  initializeBranches();
  hvvtree->GetEntry(ev);
  for (int el=0; el<intBranches.size(); el++){
    if (intBranches.at(el).second!=0) cout << intBranches.at(el).first << ":\t" << *(intBranches.at(el).second) << endl;
  }
  for (int el=0; el<floatBranches.size(); el++){
    if (floatBranches.at(el).second!=0) cout << floatBranches.at(el).first << ":\t" << *(floatBranches.at(el).second) << endl;
  }
  for (int el=0; el<vectorIntBranches.size(); el++){
    if (vectorIntBranches.at(el).second!=0){
      cout << vectorIntBranches.at(el).first << ":\t";
      for (int v=0; v<vectorIntBranches.at(el).second->size(); v++) cout << vectorIntBranches.at(el).second->at(v) << '\t';
      cout << endl;
    }
  }
  for (int el=0; el<vectorDoubleBranches.size(); el++){
    if (vectorDoubleBranches.at(el).second!=0){
      cout << vectorDoubleBranches.at(el).first << ":\t";
      for (int v=0; v<vectorDoubleBranches.at(el).second->size(); v++) cout << vectorDoubleBranches.at(el).second->at(v) << '\t';
      cout << endl;
    }
  }
  initializeBranches();
}




