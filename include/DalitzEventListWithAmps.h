/**
 * <B>FourBodyTools</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: May 2017
 *  
 **/

 
#ifndef DALITZ_EVENT_LIST_WITH_AMPS
#define DALITZ_EVENT_LIST_WITH_AMPS

#include "DalitzEventWithAmps.h"

class DalitzEventListWithAmps{
  
  std::vector<DalitzEventWithAmps> _evts;

  public:

  DalitzEventListWithAmps( DalitzEventList& evts );
  DalitzEventListWithAmps( );

  void addAmp( FitAmplitude* fitAmp );
  void addAmp( FitAmpSum* fas, std::vector<int> components );
  void addAmp( FitAmpSum* fas, int component1 );
  void addAmp( FitAmpSum* fas, int component1, int component2 );

  void setTotAmp( FitAmpSum* fas );

  int size() const;
  
  DalitzEventWithAmps& at(int i);

  virtual ~DalitzEventListWithAmps();

};


#endif

