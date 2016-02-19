import loadLib
import os
import ROOT
import sys

def trimPythia(cinput, outdir="./", pythiaLevel=1, jetAlgorithm="ak5"):
  pythiaLevel = int(pythiaLevel)
  coutput = os.path.join(outdir, "pythiaTemp.root")

  ftemp = ROOT.TFile(coutput, "recreate")
  tmpTree = ROOT.TTree("tmpTree", "")

  geneventinfoweights = ROOT.vector("double")()

  reco_GenJet_FV = [ROOT.vector("double")() for i in range(4)]
  reco_GenJet_id = ROOT.vector("int")()
  reco_GenJet_status = ROOT.vector("int")()

  reco_GenParticle_FV = [ROOT.vector("double")() for i in range(4)]
  reco_GenParticle_id = ROOT.vector("int")()
  reco_GenParticle_status = ROOT.vector("int")()

  tmpTree.Branch("genWeights", "vector<double>", ROOT.AddressOf(geneventinfoweights))

  tmpTree.Branch("reco_GenJet_X", "vector<double>", ROOT.AddressOf(reco_GenJet_FV[0]))
  tmpTree.Branch("reco_GenJet_Y", "vector<double>", ROOT.AddressOf(reco_GenJet_FV[1]))
  tmpTree.Branch("reco_GenJet_Z", "vector<double>", ROOT.AddressOf(reco_GenJet_FV[2]))
  tmpTree.Branch("reco_GenJet_E", "vector<double>", ROOT.AddressOf(reco_GenJet_FV[3]))
  tmpTree.Branch("reco_GenJet_id", "vector<int>", ROOT.AddressOf(reco_GenJet_id))
  tmpTree.Branch("reco_GenJet_status", "vector<int>", ROOT.AddressOf(reco_GenJet_status))

  tmpTree.Branch("reco_GenParticle_X", "vector<double>", ROOT.AddressOf(reco_GenParticle_FV[0]))
  tmpTree.Branch("reco_GenParticle_Y", "vector<double>", ROOT.AddressOf(reco_GenParticle_FV[1]))
  tmpTree.Branch("reco_GenParticle_Z", "vector<double>", ROOT.AddressOf(reco_GenParticle_FV[2]))
  tmpTree.Branch("reco_GenParticle_E", "vector<double>", ROOT.AddressOf(reco_GenParticle_FV[3]))
  tmpTree.Branch("reco_GenParticle_id", "vector<int>", ROOT.AddressOf(reco_GenParticle_id))
  tmpTree.Branch("reco_GenParticle_status", "vector<int>", ROOT.AddressOf(reco_GenParticle_status))



  f = ROOT.TFile(cinput, "read")
  if not f.IsZombie():
    events = f.Get("Events")
    events.SetBranchStatus("*", 0)

    if pythiaLevel == 1: suffix = "SIM"
    elif pythiaLevel == 0: suffix = "GEN"
    else:
      print "trimPythia should not be called with pythiaLevel=%i" % pythiaLevel
      assert False
    events.SetBranchStatus("recoGenJets_"+jetAlgorithm+"GenJets__"+suffix+"*", 1)
    events.SetBranchStatus("recoGenParticles_genParticles__"+suffix+"*", 1)
    events.SetBranchStatus("GenEventInfoProduct_generator__"+suffix+"*", 1)

    for ev in events:
      geneventinfoweights.clear()

      for v in range(4):
        reco_GenJet_FV[v].clear()
        reco_GenParticle_FV[v].clear()
      reco_GenJet_id.clear()
      reco_GenParticle_id.clear()
      reco_GenJet_status.clear()
      reco_GenParticle_status.clear()

      geneventinfoWrapper = getattr(ev, "GenEventInfoProduct_generator__"+suffix+".")
      genpdf = geneventinfoWrapper.product().pdf()
      genweights = geneventinfoWrapper.product().weights()
      for genweight in genweights: geneventinfoweights.push_back(genweight)

      jetWrapper = getattr(ev, "recoGenJets_"+jetAlgorithm+"GenJets__"+suffix+".")
      genparticleWrapper = getattr(ev, "recoGenParticles_genParticles__"+suffix+".")
      reco_GenJets = jetWrapper.product()
      reco_GenParticles = genparticleWrapper.product()

      for part in reco_GenParticles:
        if (
          ((part.pdgId()==25 or part.pdgId()==32) and part.status()==22) # Generated Higgs
          or
          ((abs(part.pdgId())>=11 and abs(part.pdgId())<=16) and (part.status()==23 or part.status()==1)) # Generated leptons
          or
          ((abs(part.pdgId())<=6 or abs(part.pdgId())==21) and (part.status()==21 or part.status()==23 or part.status()==22)) # Generated partons
          ):
          reco_GenParticle_FV[0].push_back(part.px())
          reco_GenParticle_FV[1].push_back(part.py())
          reco_GenParticle_FV[2].push_back(part.pz())
          reco_GenParticle_FV[3].push_back(part.energy())
          reco_GenParticle_id.push_back(part.pdgId())
          reco_GenParticle_status.push_back(part.status())

      reco_GenParticle_duplicates = findDuplicates(reco_GenParticle_FV, reco_GenParticle_id, reco_GenParticle_status)

      # Add status==1 or status==23 (with precedence for status==1) leptons to genJet collection
      for p, id in enumerate(reco_GenParticle_id):
        if (abs(id)>=11 and abs(id) <= 16 and abs(id) % 2 == 1):
          match = False
          for iOriginal, iDuplicate in reco_GenParticle_duplicates:
            if iOriginal==id:
              match = True
              break
          if match: continue
          for fv in range(4): reco_GenJet_FV[fv].push_back(reco_GenParticle_FV[fv].at(p))
          reco_GenJet_id.push_back(id)
          reco_GenJet_status.push_back(reco_GenParticle_status.at(p))

      removalArray = []
      for iFirst, iTransfer in reco_GenParticle_duplicates:
        # iTransfer is the status==1 particle
        inserted = False
        for it, iIndex in enumerate(removalArray):
          if iTransfer > iIndex:
            removalArray.insert(it, iTransfer)
            inserted = True
            break
        if not inserted: removalArray.append(iTransfer)
      # Remove status==1 duplicates from genParticles
      for iTransfer in removalArray:
        for fv in reco_GenParticle_FV:
          fv.erase(fv.begin()+iTransfer)
        reco_GenParticle_id.erase(reco_GenParticle_id.begin()+iTransfer)
        reco_GenParticle_status.erase(reco_GenParticle_status.begin()+iTransfer)

      for jet in reco_GenJets:
        reco_GenJet_FV[0].push_back(jet.px())
        reco_GenJet_FV[1].push_back(jet.py())
        reco_GenJet_FV[2].push_back(jet.pz())
        reco_GenJet_FV[3].push_back(jet.energy())
        reco_GenJet_id.push_back(jet.pdgId())
        reco_GenJet_status.push_back(jet.status())

      tmpTree.Fill()
    f.Close()

  elif f.IsOpen(): f.Close()
  ftemp.WriteTObject(tmpTree)
  del tmpTree
  ftemp.Close()

def findDuplicates(fourvectors, ids, statuses):
  #fourvectors = [pxs, pys, pzs, Es], convert to list of 4 momenta
  tlorentzvectors = [ROOT.TLorentzVector(*pxpypzE) for pxpypzE in zip(*fourvectors)]
  assert len(tlorentzvectors) == len(ids) == len(statuses)
  duplicates = []

  for xx, (p_tm, id_tm, st_tm) in enumerate(zip(tlorentzvectors, ids, statuses)):
    if st_tm==23 or st_tm==1:
      for yy, (p_tbm, id_tbm, st_tbm) in enumerate(zip(tlorentzvectors[xx+1:], ids[xx+1:], statuses[xx+1:]), start=xx+1):
        if (id_tbm==id_tm and st_tbm!=st_tm and (st_tbm==1 or st_tbm==23)):
          dot_tm = p_tbm.Dot(p_tm)
          diff_sqmass = abs(dot_tm - p_tm.M2())
          if (diff_sqmass<0.005*abs(p_tm.M2())):
            iFirst = (xx if st_tm==23 else yy) #Should be order-independent
            iSecond = (yy if st_tbm==1 else xx)
            duplicates.append((iFirst, iSecond))

  return duplicates

if __name__ == "__main__":
    trimPythia(*sys.argv[1:5])
