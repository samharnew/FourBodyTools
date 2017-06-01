
#ifndef PHILIPPE_MODEL_WRAPPER_HH
#define PHILIPPE_MODEL_WRAPPER_HH

#include "TRandom.h"
#include "TString.h"
#include "Mint/FitParameter.h"
#include "Mint/NamedParameter.h"
#include "Mint/Minimiser.h"
#include "Mint/Neg2LL.h"
#include "Mint/Neg2LLSum.h"
#include "Mint/DalitzEventList.h"
#include "Mint/NamedDecayTreeList.h"
#include "Mint/DecayTree.h"
#include "Mint/DiskResidentEventList.h"
#include "Mint/CLHEPPhysicalConstants.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "Mint/PdfBase.h"
#include "Mint/DalitzPdfBase.h"
#include "Mint/DalitzPdfBaseFastInteg.h"
#include "Mint/FitAmplitude.h"
#include "Mint/FitAmpSum.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/DalitzEvent.h"
#include "Mint/AmpRatios.h"
#include "Mint/IEventGenerator.h"
#include "Mint/DalitzBWBoxSet.h"
#include "Mint/DalitzBoxSet.h"
#include "Mint/SignalGenerator.h"
#include "Mint/FromFileGenerator.h"
#include "Mint/DalitzSumPdf.h"
#include "Mint/cexp.h"
#include "Mint/DalitzPdfNormChecker.h"
#include "Mint/IFastAmplitudeIntegrable.h"
#include "Mint/DalitzPdfBaseFlexiFastInteg.h"
#include "Mint/DalitzPdfSaveInteg.h"
#include "Mint/Chi2Binning.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/FitAmpList.h"
#include "Mint/DalitzPdfBaseMCInteg.h"

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "AmpRatio.h"

#include "TRandom3.h"

class PhilippeModelWrapper{
  
  DalitzEventList eventListMC;

  DalitzEventPattern pat;

  MINT::counted_ptr<FitAmpSum> fas;
  MINT::counted_ptr<FitAmpSum> fascc;

  MINT::FitParameter a1_re ;
  MINT::FitParameter a1_im ;
  MINT::FitParameter pi_re ;
  MINT::FitParameter pi_im ;
  MINT::FitParameter pi2_re;
  MINT::FitParameter pi2_im;
  MINT::FitParameter a1p_re;
  MINT::FitParameter a1p_im;

  MINT::counted_ptr<FitAmpList> copyList;

  MINT::counted_ptr<FitAmpList> a1List ;
  MINT::counted_ptr<FitAmpList> piList ;
  MINT::counted_ptr<FitAmpList> pi2List;
  MINT::counted_ptr<FitAmpList> a1pList;

  MINT::counted_ptr<MINT::IReturnComplex> ar_a1 ;
  MINT::counted_ptr<MINT::IReturnComplex> ar_pi ;
  MINT::counted_ptr<MINT::IReturnComplex> ar_pi2;
  MINT::counted_ptr<MINT::IReturnComplex> ar_a1p;

  FitAmpSum* fas_a1 ;
  FitAmpSum* fas_pi ;
  FitAmpSum* fas_pi2;
  FitAmpSum* fas_a1p;

  MINT::NamedParameter<std::string> ModelName;
  
  static PhilippeModelWrapper* s_wrapper; 

  public:
  
  PhilippeModelWrapper();

  static PhilippeModelWrapper* getStaticWrapper();

  MINT::counted_ptr<FitAmpSum> getFitAmpSum();
  MINT::counted_ptr<FitAmpSum> getFitAmpSumCC();

  
  TString getModelName();

  void radomiseModel(int seed);
  void resetModel();


  ~PhilippeModelWrapper();
  

};







#endif 

