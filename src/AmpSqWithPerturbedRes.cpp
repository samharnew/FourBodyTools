#include "AmpSqWithPerturbedRes.h"



AmpSqWithPerturbedRes::AmpSqWithPerturbedRes(){



}





  
double AmpSqWithPerturbedRes::getVal( 
  	                double norm   , // this is |A(x)|^2
  	                double phi    , // this is arg(A(x)/B(x))
  	                double ratio  , // this is |B(x)| / |A(x)|
  	                double realDel, // this is the real part of the purtubation
  	                double imagDel  // this is the imaginary part of the purtubation
  	            ) const{

  
  double magDelSq = realDel*realDel + imagDel*imagDel;

  double unpertubedVal = 1.0 + ratio*ratio + 2.0*ratio*cos(phi);
  double pertubation   = (magDelSq + 2.0*realDel)*ratio*ratio 
                       + 2.0*ratio*( cos(phi)*realDel + sin(phi)*imagDel );

  return norm*(unpertubedVal + pertubation);

}

double AmpSqWithPerturbedRes::getUnperturbedProb( double phi, double ratio, double realDel, double imagDel, double globAys) const{
  double pertProb   = getVal(globAys, phi, ratio, realDel, imagDel);
  double unpertProb = getVal(1.0, phi, ratio, 0.0, 0.0);
  return unpertProb/(pertProb + unpertProb);
}

double AmpSqWithPerturbedRes::getPerturbedProb( double phi, double ratio, double realDel, double imagDel, double globAys) const{
  double pertProb   = getVal(globAys, phi, ratio, realDel, imagDel);
  double unpertProb = getVal(1.0, phi, ratio, 0.0, 0.0);
  return pertProb/(pertProb + unpertProb);
}


AmpSqWithPerturbedRes::~AmpSqWithPerturbedRes(){


	
}






