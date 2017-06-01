/**
 * <B>FourBodyTools</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: May 2017
 *  
 **/

 
#ifndef MOD_INSP_TWODIM_BINNING_HH
#define MOD_INSP_TWODIM_BINNING_HH

#include "TLorentzVector.h"
#include "TRandom.h"
#include <iostream>
#include "Mint/FitAmpSum.h"
#include "Mint/DalitzEventPattern.h"
#include "Mint/counted_ptr.h"

class ModInspired2DBinning{
  
  DalitzEventPattern           _pat;
  MINT::counted_ptr<FitAmpSum> _fasA;    
  MINT::counted_ptr<FitAmpSum> _fasB;    

  public:
  
  ModInspired2DBinning(MINT::counted_ptr<FitAmpSum> fasA, MINT::counted_ptr<FitAmpSum> fasB, DalitzEventPattern& pat);
  
  double getPhaseDiff( DalitzEvent& evt ) const;
  double getAmpRatio ( DalitzEvent& evt ) const;

  int getBinNumX( DalitzEvent& evt ) const;
  int getBinNumY( DalitzEvent& evt ) const;


  virtual int getBinNum( DalitzEvent& evt ) const;

  virtual ~ModInspired2DBinning();

};


#endif

