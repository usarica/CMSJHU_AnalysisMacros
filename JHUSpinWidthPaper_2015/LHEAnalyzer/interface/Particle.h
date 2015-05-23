#include <vector>
#include "TLorentzVector.h"
//#include <string>
//#include "TString.h"

class Particle{
public:

  // Constructors

  Particle();
  Particle(int id_, TLorentzVector p4_);
  Particle(const Particle& particle_);
  Particle& operator=(const Particle& particle_);
  ~Particle(){};

  // Data

  int id;
  TLorentzVector p4;
  bool passSelection;


  // Member functions
  void setSelected(bool isSelected=true){ passSelection = isSelected; }

  void addMother(Particle* myParticle){ mothers.push_back(myParticle); };
  void addDaughter(Particle* myParticle){ daughters.push_back(myParticle); };

  int getNMothers() const{ return mothers.size(); };
  int getNDaughters() const{ return daughters.size(); };

  Particle* getMother(int index) const;
  Particle* getDaughter(int index) const;

  double charge()const;
  double m()const{ return p4.M(); }
  double x()const{ return p4.X(); }
  double y()const{ return p4.Y(); }
  double z()const{ return p4.Z(); }
  double t()const{ return p4.T(); }
  double pt()const{ return p4.Pt(); }
  double eta()const{ return p4.Eta(); }
  double phi()const{ return p4.Phi(); }
  double rapidity()const{ return p4.Rapidity(); }
  double deltaR(const TLorentzVector& v)const{ return p4.DeltaR(v); }


protected:
  std::vector<Particle*> mothers;
  std::vector<Particle*> daughters;
};

