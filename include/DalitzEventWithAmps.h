/**
 * <B>FourBodyTools</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: May 2017
 *  
 **/

 
#ifndef DALITZ_EVENT_WITH_AMPS
#define DALITZ_EVENT_WITH_AMPS
#include "TLorentzVector.h"
#include "TRandom.h"
#include <iostream>
#include "Mint/DalitzEvent.h"
#include "Mint/FitAmplitude.h"
#include "Mint/FitAmpSum.h"


class DalitzEventWithAmps{
  
  DalitzEvent _evt; 
  std::vector< std::complex<double> > _amps;
  std::complex<double> _totAmp; 

  public:

  DalitzEventWithAmps(DalitzEvent& evt);
  
  void addAmp( std::complex<double> amp );
  void addAmp( FitAmplitude* fitAmp );
  void addAmp( FitAmpSum* fas, std::vector<int> components );
  void addAmp( FitAmpSum* fas, int component1 );
  void addAmp( FitAmpSum* fas, int component1, int component2 );


  void setTotAmp( std::complex<double> amp );
  void setTotAmp( FitAmpSum* fas );


  std::complex<double> getAmp( int i ) const;
  
  std::complex<double> getAmpSum() const;

  std::complex<double> getTotAmp() const;
  
  std::complex<double> getAmpRelToTot(int i) const;


  double getProb() const;
  double getProbWithPerturb(int amp, std::complex<double> peturb) const;

  int getNumAmps() const;

  virtual ~DalitzEventWithAmps();

};


#endif

