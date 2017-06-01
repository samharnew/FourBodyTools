#include "DalitzEventListWithAmps.h"


DalitzEventListWithAmps::DalitzEventListWithAmps( DalitzEventList& evts ){
  
  for (unsigned i = 0; i < evts.size(); i++){
  	DalitzEventWithAmps evt( evts[i] );
  	_evts.push_back( evt );
  }

}

DalitzEventListWithAmps::DalitzEventListWithAmps( ){

}

void DalitzEventListWithAmps::addAmp( FitAmplitude* fitAmp ){
  
  for (int i = 0; i < size(); i++){
  	_evts.at(i).addAmp( fitAmp );
  }

}

void DalitzEventListWithAmps::addAmp( FitAmpSum* fas, std::vector<int> components ){

  for (int i = 0; i < size(); i++){
  	_evts.at(i).addAmp( fas, components );
  }

}

void DalitzEventListWithAmps::addAmp( FitAmpSum* fas, int component1 ){

  for (int i = 0; i < size(); i++){
  	_evts.at(i).addAmp( fas, component1 );
  }

}

void DalitzEventListWithAmps::addAmp( FitAmpSum* fas, int component1, int component2 ){

  for (int i = 0; i < size(); i++){
  	_evts.at(i).addAmp( fas, component1, component2 );
  }

}


void DalitzEventListWithAmps::setTotAmp( FitAmpSum* fas ){

  for (int i = 0; i < size(); i++){
  	_evts.at(i).setTotAmp( fas );
  }

}

DalitzEventWithAmps& DalitzEventListWithAmps::at(int i){
  return _evts.at(i);
}


int DalitzEventListWithAmps::size() const{
  return _evts.size();
}


DalitzEventListWithAmps::~DalitzEventListWithAmps(){

}


