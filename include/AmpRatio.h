/**
 * <B>D4piCleoAnalysis</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Feb 2016
 *
 *  
 **/

 
#ifndef MINT_AMP_RATIO
#define MINT_AMP_RATIO

#include "Mint/IReturnComplex.h"
#include "Mint/FitParameter.h"
#include <complex.h>

class AmpRatio : virtual public MINT::IReturnComplex{

  MINT::FitParameter& _re;
  MINT::FitParameter& _im;

  public:

  AmpRatio(MINT::FitParameter& re, MINT::FitParameter& im);
    
  std::complex<double> ComplexVal();
  
};


#endif

