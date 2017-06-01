/**
 * <B>FourBodyTools</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: May 2017
 *  
 **/

 
#ifndef FOUR_BODY_MOM_FROM_HEL_ANGLE_HH
#define FOUR_BODY_MOM_FROM_HEL_ANGLE_HH
#include "TLorentzVector.h"
#include "TRandom.h"
#include <iostream>


class FourBodyMomFromHelAng{
  
  //four-momentum of mother and daughters
  TLorentzVector _p0; 
  TLorentzVector _p1;
  TLorentzVector _p2;
  TLorentzVector _p3;
  TLorentzVector _p4;
  
  //helicity variables that describe the four-body phase space point
  double _m12;
  double _m34;
  double _theta12;
  double _theta34;
  double _phi;
  
  //mass, and mass sq, of the mother and daughters
  double _m0, _m1, _m2, _m3, _m4;
  double _m0sq, _m1sq, _m2sq, _m3sq, _m4sq;

  //Stuff that's useful to cache for generating toy MC 
  double _m12Min;
  double _m12Max;
  double _m34Min;
  double _m34Max;
  
  double _psMax;

  public:

  FourBodyMomFromHelAng(double m0, double m1, double m2, double m3, double m4);
  
  double coshel(TLorentzVector particle, TLorentzVector parent, TLorentzVector grandparent);
  double angleBetweenDecayPlanes(TLorentzVector particleA1, TLorentzVector particleA2, TLorentzVector particleB1, TLorentzVector particleB2);

  double twoBodyBreakupMom        (double m0, double m1, double m2);
  double twoBodyBreakupMomSq      (double m0, double m1, double m2);
  double twoBodyBreakupMomSqFromSq(double m0sq, double m1sq, double m2sq);

  double m12Min(){ return _m12Min; }
  double m12Max(){ return _m12Max; }
  double m34Min(){ return _m34Min; }
  double m34Max(){ return _m34Max; }

  bool set4Vects(TLorentzVector& p0, TLorentzVector& p1, TLorentzVector& p2, TLorentzVector& p3, TLorentzVector& p4);
  bool set4Vects(std::vector<TLorentzVector> fourVects);

  bool setHelVars(double m12, double m34, double coshel12, double coshel34, double phi);
  bool setHelVars(std::vector<double> helVars);
  
  void printHelVars();

  TLorentzVector& get4Vect(int i);
  std::vector<TLorentzVector> get4Vects();

  double getHelVar(int i);
  std::vector<double> getHelVars();

  double getPhaseSpace();
  int generatePhaseSpace(TRandom* rand);


};


#endif




