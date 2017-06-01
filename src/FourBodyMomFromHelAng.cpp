#include "FourBodyMomFromHelAng.h"




FourBodyMomFromHelAng::FourBodyMomFromHelAng(double m0, double m1, double m2, double m3, double m4):
    _m0(m0),
    _m1(m1),
    _m2(m2),
    _m3(m3),
    _m4(m4),
    _m0sq(m0*m0),
    _m1sq(m1*m1),
    _m2sq(m2*m2),
    _m3sq(m3*m3),
    _m4sq(m4*m4)
{

  _m12Min = _m1 + _m2;
  _m12Max = _m0 - _m3 - _m4;
  _m34Min = _m3 + _m4;
  _m34Max = _m0 - _m3 - _m4;
  
  //Must be smaller than this! But can't be bigger at least.
  _psMax = twoBodyBreakupMom(_m0, _m12Min, _m34Min)*twoBodyBreakupMom(_m12Max,_m1,_m2)*twoBodyBreakupMom(_m34Max,_m3,_m4);

}


double FourBodyMomFromHelAng::coshel(TLorentzVector particle, TLorentzVector parent, TLorentzVector grandparent) {
    
    TVector3 boosttoparent = -(parent.BoostVector());
     
    particle.Boost(boosttoparent);
    grandparent.Boost(boosttoparent);
     
    TVector3 particle3 = particle.Vect();
    TVector3 grandparent3 = grandparent.Vect();
    double numerator = particle3.Dot(grandparent3);
    double denominator = (particle3.Mag())*(grandparent3.Mag());

    double temp = numerator/denominator;
     
    return temp;
  
}
  

double FourBodyMomFromHelAng::angleBetweenDecayPlanes(TLorentzVector particleA1, TLorentzVector particleA2, TLorentzVector particleB1, TLorentzVector particleB2) {
    
  TLorentzVector mother = particleA1 + particleA2 + particleB1 + particleB2;
  
  TVector3 boosttomother = -(mother.BoostVector());
  particleA1.Boost(boosttomother);
  particleA2.Boost(boosttomother);
  particleB1.Boost(boosttomother);
  particleB2.Boost(boosttomother);
  
  TVector3 particleA1vect = particleA1.Vect();
  TVector3 particleA2vect = particleA2.Vect();
  TVector3 particleB1vect = particleB1.Vect();
  TVector3 particleB2vect = particleB2.Vect();
  
  TVector3 Avect = ((particleA1 + particleA2).Vect()).Unit();
  TVector3 Bvect = ((particleB1 + particleB2).Vect()).Unit();
  
  TVector3 normAvect = (particleA1vect.Cross(particleA2vect)).Unit();
  TVector3 normBvect = (particleB1vect.Cross(particleB2vect)).Unit();
  
  double cosPhi = normAvect.Dot(normBvect);
  double sinPhi = (normAvect.Cross(normBvect)).Dot(Bvect);

  double phi = atan2(sinPhi,cosPhi);
   
  return phi;
  
}
  


double FourBodyMomFromHelAng::twoBodyBreakupMomSqFromSq(double m0sq, double m1sq, double m2sq){

  double term1 = m2sq - m0sq - m1sq;
  return (term1*term1)/(4.0*m0sq) - m1sq;

}

double FourBodyMomFromHelAng::twoBodyBreakupMomSq(double m0, double m1, double m2){

  if ((m0 - m1 - m2) < 0.0 ) return -99999.9;
  return twoBodyBreakupMomSqFromSq(m0*m0, m1*m1, m2*m2);

}

double FourBodyMomFromHelAng::twoBodyBreakupMom(double m0, double m1, double m2){

  if ((m0 - m1 - m2) < 0.0 ) return -99999.9;
  return sqrt( twoBodyBreakupMomSq(m0, m1, m2) );
  
}


double FourBodyMomFromHelAng::getPhaseSpace(){
  
  double m12sq = _m12*_m12;
  double m34sq = _m34*_m34;

  double p   = twoBodyBreakupMomSqFromSq( _m0sq , m12sq , m34sq );
  double p12 = twoBodyBreakupMomSqFromSq( m12sq , _m1sq , _m2sq );
  double p34 = twoBodyBreakupMomSqFromSq( m34sq , _m3sq , _m4sq );
  
  if (p < 0 || p12 < 0 || p34 < 0) return 0.0;

  return sqrt(p*p12*p34);

}

int FourBodyMomFromHelAng::generatePhaseSpace(TRandom* rand){

  _theta12 = rand->Uniform(-1, 1);
  _theta34 = rand->Uniform(-1, 1);
  _phi     = rand->Uniform(-TMath::Pi(), TMath::Pi());
  
  int iterations = 0;
  while (1==1){
    
    iterations++;

    _m12 = rand->Uniform( _m12Min, _m12Max );
    _m34 = rand->Uniform( _m34Min, _m34Max );
  
    if ( _m12 + _m34 > _m0 ) continue;
    
    double trialPS = rand->Uniform( 0, _psMax ); 
    if ( getPhaseSpace() > trialPS ) break;
  }
  
  setHelVars(_m12, _m34, _theta12, _theta34, _phi);

  return iterations;
}

bool FourBodyMomFromHelAng::set4Vects(TLorentzVector& p0, TLorentzVector& p1, TLorentzVector& p2, TLorentzVector& p3, TLorentzVector& p4){

  _p0 = p0;
  _p1 = p1;
  _p2 = p2;
  _p3 = p3;
  _p4 = p4;
 
  _m12     = (p1 + p2).M();
  _m34     = (p3 + p4).M();
  _theta12 = coshel(p1, p1+p2, p0);
  _theta34 = coshel(p3, p3+p4, p0);
  _phi     = angleBetweenDecayPlanes(p1, p2, p3, p4);
  

  return true;
  
}


bool FourBodyMomFromHelAng::setHelVars(double m12, double m34, double coshel12, double coshel34, double phi){

  _m12      = m12;
  _m34      = m34;
  _theta12  = coshel12;
  _theta34  = coshel34;
  _phi      = phi;

  double p1sq = twoBodyBreakupMomSq(m12, _m1, _m2);
  double p3sq = twoBodyBreakupMomSq(m34, _m3, _m4);
  double p1 = sqrt(p1sq);
  double p3 = sqrt(p3sq);

  if (p1sq == -99999.9) return false;
  if (p3sq == -99999.9) return false;
  
  double p1x = -p1*coshel12;
  double p1y = -p1*sqrt(1.0 - coshel12*coshel12);
  double p1e =  sqrt(p1sq + _m1sq);
  
  double p2x = -p1x;
  double p2y = -p1y;
  double p2e = sqrt(p1sq + _m2sq);  
  
  double p3x = p3*coshel34;
  double p3y = p3*sqrt(1.0 - coshel34*coshel34);
  double p3e = sqrt(p3sq + _m3sq);  
  
  double p4x = -p3x;
  double p4y = -p3y;  
  double p4e = sqrt(p3sq + _m4sq);  
  
  double pAsq = twoBodyBreakupMomSq(_m0, m12, m34);
  double pA   = sqrt(pAsq);

  if (pAsq == -99999.9) return false;
  
  double pAe = sqrt(pAsq + m12*m12);  
  
  double pB  = -pA;
  double pBe = sqrt(pAsq + m34*m34);  
  
  TLorentzVector fourVectA(pA, 0, 0, pAe);
  TLorentzVector fourVectB(pB, 0, 0, pBe);
  
  TLorentzVector fourVect1(p1x, p1y*cos(phi), p1y*sin(phi), p1e);
  TLorentzVector fourVect2(p2x, p2y*cos(phi), p2y*sin(phi), p2e);  
  TLorentzVector fourVect3(p3x, p3y, 0, p3e);
  TLorentzVector fourVect4(p4x, p4y, 0, p4e);   
  
  fourVect1.Boost(fourVectA.BoostVector());
  fourVect2.Boost(fourVectA.BoostVector());
  
  fourVect3.Boost(fourVectB.BoostVector());
  fourVect4.Boost(fourVectB.BoostVector());
  
  _p0 = TLorentzVector(0, 0, 0, _m0);
  _p1 = fourVect1;
  _p2 = fourVect2;
  _p3 = fourVect3;
  _p4 = fourVect4;  

  return true;

  
}


bool FourBodyMomFromHelAng::setHelVars(std::vector<double> helVars){

  return setHelVars( helVars.at(0), helVars.at(1), helVars.at(2), helVars.at(3), helVars.at(4) );

}

bool FourBodyMomFromHelAng::set4Vects(std::vector<TLorentzVector> fourVects){

  return set4Vects( fourVects.at(0), fourVects.at(1), fourVects.at(2), fourVects.at(3), fourVects.at(4) );

}

TLorentzVector& FourBodyMomFromHelAng::get4Vect(int i){
  if (i == 0) return _p0;
  if (i == 1) return _p1;
  if (i == 2) return _p2;
  if (i == 3) return _p3;
  if (i == 4) return _p4;
  return _p0;
}

double FourBodyMomFromHelAng::getHelVar(int i){
  if (i == 0) return _m12    ;
  if (i == 1) return _m34    ;
  if (i == 2) return _theta12;
  if (i == 3) return _theta34;
  if (i == 4) return _phi    ;
  return _m12;
}

std::vector<TLorentzVector> FourBodyMomFromHelAng::get4Vects(){
  std::vector<TLorentzVector> vects;
  vects.push_back(_p0);
  vects.push_back(_p1);
  vects.push_back(_p2);
  vects.push_back(_p3);
  vects.push_back(_p4);
  return vects;
}

std::vector<double> FourBodyMomFromHelAng::getHelVars(){
  std::vector<double> vects;
  vects.push_back(_m12    );
  vects.push_back(_m34    );
  vects.push_back(_theta12);
  vects.push_back(_theta34);
  vects.push_back(_phi    );
  return vects;
}


void FourBodyMomFromHelAng::printHelVars(){

  std::cout << "( " << _m12 << ", " << _m34 << ", " << _theta12 << ", " << _theta34 << ", " << _phi << " )" << std::endl;


}


