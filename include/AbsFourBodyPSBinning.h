/**
 * <B>FourBodyTools</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: May 2017
 *  
 **/

 
#ifndef ABS_FOUR_BODY_PS_BINNING_HH
#define ABS_FOUR_BODY_PS_BINNING_HH
#include "TLorentzVector.h"
#include "TRandom.h"
#include <iostream>


class AbsFourBodyPSBinning{
  
  int _nBins;

  public:

  AbsFourBodyPSBinning();
  
  virtual int getBinNum( DalitzEvent& evt ) const;
  
  int getNumBins() const;

  virtual ~AbsFourBodyPSBinning();

};


#endif

