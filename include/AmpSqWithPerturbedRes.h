

/**
 * <B>FourBodyTools</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: May 2017
 *  
 **/

 
#ifndef AMP_SQ_W_PERTERB_RES_HH
#define AMP_SQ_W_PERTERB_RES_HH

#include <iostream>
#include <math.h>


//This is an amplude T(x) = A(x) + (1 + del) B(x) 
//where T(x) is the total amplitude, B(x) the amplude
//of a single resonance, and A(x) is the remaining
//resonances. A small complex pertuabtion delta has
//been added to B(x). 

class AmpSqWithPerturbedRes{
  
  public:

  AmpSqWithPerturbedRes();
  
  double getVal( 
  	                double norm   , // this is |A(x)|^2
  	                double phi    , // this is arg(A(x)/B(x))
  	                double ratio  , // this is |B(x)| / |A(x)|
  	                double realDel, // this is the real part of the purtubation
  	                double imagDel  // this is the imaginary part of the purtubation
  	            ) const;
  

  
  double getUnperturbedProb( double phi, double ratio, double realDel, double imagDel, double globAys = 1.0) const;

  double getPerturbedProb( double phi, double ratio, double realDel, double imagDel, double globAys = 1.0) const;


  ~AmpSqWithPerturbedRes();

};


#endif


