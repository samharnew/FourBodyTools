#include "DalitzEventWithAmps.h"

DalitzEventWithAmps::DalitzEventWithAmps(DalitzEvent& evt) :
  _evt   (evt),
  _totAmp(0,0)
{


}

void DalitzEventWithAmps::addAmp( std::complex<double> amp ){
  _amps.push_back(amp);
}

void DalitzEventWithAmps::addAmp( FitAmplitude* fitAmp ){
  _amps.push_back( fitAmp->getVal(_evt) );
}

void DalitzEventWithAmps::addAmp( FitAmpSum* fas, std::vector<int> components ){
  std::complex<double> sum(0.0,0.0);
  for (auto i = components.begin(); i != components.end(); ++i ){
    sum += fas->getAmpPtr(*i)->getVal(_evt);
  }    
  addAmp(sum);
}

void DalitzEventWithAmps::addAmp( FitAmpSum* fas, int component1 ){
  std::vector<int> comp;
  comp.push_back(component1);
  addAmp(fas, comp);
}

void DalitzEventWithAmps::addAmp( FitAmpSum* fas, int component1, int component2 ){
  std::vector<int> comp;
  comp.push_back(component1);
  comp.push_back(component2);
  addAmp(fas, comp);
}

void DalitzEventWithAmps::setTotAmp( std::complex<double> amp ){
  _totAmp = amp;
}

void DalitzEventWithAmps::setTotAmp( FitAmpSum* fas ){
  setTotAmp( fas->getVal(_evt) );
}



std::complex<double> DalitzEventWithAmps::getAmp( int i ) const{
  return _amps.at(i);
}

std::complex<double> DalitzEventWithAmps::getAmpSum() const{

  std::complex<double> sum(0.0,0.0);
  for (auto i = _amps.begin(); i != _amps.end(); ++i ){
    sum += *i;
  }
  return sum;

}

std::complex<double> DalitzEventWithAmps::getTotAmp() const{
  return _totAmp;
}

int DalitzEventWithAmps::getNumAmps() const{
  return _amps.size();
}

double DalitzEventWithAmps::getProb() const{
  return norm(_totAmp);
}

std::complex<double> DalitzEventWithAmps::getAmpRelToTot(int i) const{
  return getAmp(i)/getTotAmp();
}

double DalitzEventWithAmps::getProbWithPerturb(int amp, std::complex<double> peturb) const{
  return norm(_totAmp + getAmp(amp)*peturb);
}

DalitzEventWithAmps::~DalitzEventWithAmps(){

}



