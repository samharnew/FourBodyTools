 
#include "ModInspired2DBinning.h"





ModInspired2DBinning::ModInspired2DBinning(MINT::counted_ptr<FitAmpSum> fasA, MINT::counted_ptr<FitAmpSum> fasB, DalitzEventPattern& pat):
  _pat (pat ),
  _fasA(fasA),
  _fasB(fasB)
{
  
}
  

double ModInspired2DBinning::getPhaseDiff( DalitzEvent& evt ) const{
  std::complex<double> valA = _fasA->getVal( evt );
  std::complex<double> valB = _fasB->getVal( evt );
  return arg(valA*conj(valB));
}

double ModInspired2DBinning::getAmpRatio ( DalitzEvent& evt ) const{
  double valA = _fasA->Prob( evt );
  double valB = _fasB->Prob( evt );
  return sqrt(valB/valA);
}

int ModInspired2DBinning::getBinNumX( DalitzEvent& evt ) const{
  return 0;
}
int ModInspired2DBinning::getBinNumY( DalitzEvent& evt ) const{
  return 0;	
}

int ModInspired2DBinning::getBinNum( DalitzEvent& evt ) const{
  return 0;
}

ModInspired2DBinning::~ModInspired2DBinning(){


}

