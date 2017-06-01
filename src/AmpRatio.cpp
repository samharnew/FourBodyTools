#include "AmpRatio.h"

AmpRatio::AmpRatio(MINT::FitParameter& re, MINT::FitParameter& im) : 
  _re(re),
  _im(im)
{


}
    
std::complex<double> AmpRatio::ComplexVal(){

  std::complex<double> result(_re,_im); 
  return result;
  
}
