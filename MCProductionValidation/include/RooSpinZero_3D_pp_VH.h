#ifndef ROOSPINZERO_3D_PP_VH
#define ROOSPINZERO_3D_PP_VH

#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/RooSpinZero.h>
#else
#include "RooSpinZero.h"
#endif
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

class RooSpinZero_3D_pp_VH : public RooSpinZero {
public:

  Double_t sqrts;

  RooSpinZero_3D_pp_VH(){}
  RooSpinZero_3D_pp_VH(
    const char *name, const char *title,
    modelMeasurables _measurables,
    modelParameters _parameters,
    Double_t _sqrts,
    int _Vdecay1=1, int _Vdecay2=1
    );

  RooSpinZero_3D_pp_VH(const RooSpinZero_3D_pp_VH& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooSpinZero_3D_pp_VH(*this, newname); }
  inline virtual ~RooSpinZero_3D_pp_VH(){}

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:

  Double_t evaluate() const;

  void evaluatePolarizationTerms(Double_t& A00term, Double_t& Appterm, Double_t& Ammterm, Double_t& A00ppterm, Double_t& A00mmterm, Double_t& Appmmterm, const Int_t code, bool isGammaV1=false, bool isGammaV2=false) const;

  Double_t partonicLuminosity(Double_t mVal, Double_t YVal, Double_t sqrts) const;

};

#endif
