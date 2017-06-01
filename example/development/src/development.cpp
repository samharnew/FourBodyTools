#include <iostream>
#include <vector>

#include "PhilippeModelWrapper.h"
#include "DalitzEventListWithAmps.h"
#include "ModInspired2DBinning.h"
#include "AmpSqWithPerturbedRes.h"
#include "ModelInspCPVFitter.h"
#include "HyperPoint.h"
#include "ModelInspCPAsymmetry.h"
#include "HyperPlotStyle.h"


const int comp_a1260pToSigma  = 0; 
const int comp_a1640mToSigma  = 1; 
const int comp_pi1300p        = 2; 
const int comp_pi1300m        = 3; 
const int comp_pi1670pToSigma = 4; 
const int comp_rhoSigma       = 5; 
const int comp_f0Sigma        = 6; 
const int comp_a1260pToRho    = 7; 
const int comp_a1640pToRho    = 8; 
const int comp_rhoRhoS        = 9; 
const int comp_rhoRhoP        = 10; 
const int comp_rhoRhoD        = 11; 
const int comp_f2f2           = 12;  
const int comp_pi1670pTof2    = 13; 
const int comp_a1260mToSigma  = 14; 
const int comp_a1260mToRho    = 15; 






void generateToyMC(int seed){

  int nEvents = 100000;
  
  PhilippeModelWrapper* fasPhilWrapper = PhilippeModelWrapper::getStaticWrapper();
  fasPhilWrapper->resetModel();

  TString outdir = "toyMC/";
  outdir += fasPhilWrapper->getModelName();
  outdir += "_";
  outdir += seed;
  outdir += ".root";

  MINT::counted_ptr<FitAmpSum> fas   = fasPhilWrapper->getFitAmpSum  ();    

  DalitzEventPattern pat(421, -211, 211, 211, -211);

  TRandom3 randomGen(seed);

  SignalGenerator signalGenerator(pat, fas.get(), &randomGen);
  DalitzEventList eventList;
  signalGenerator.FillEventList(eventList, nEvents);
  eventList.save( outdir.Data() );

}

void generateCPVToyMC(int seed, int comp, std::complex<double> delta, TString uniqueTag ){

  int nEvents = 100000;
  
  PhilippeModelWrapper* fasPhilWrapper = PhilippeModelWrapper::getStaticWrapper();
  fasPhilWrapper->resetModel();
  MINT::counted_ptr<FitAmpSum> fas   = fasPhilWrapper->getFitAmpSum  ();    

  TString outdir = "toyCPVMC/";
  outdir += fasPhilWrapper->getModelName();
  outdir += "_";
  outdir += uniqueTag;
  outdir += "_comp-";
  outdir += comp;
  outdir += "_seed-";
  outdir += seed;
  outdir += ".root";
  
  std::complex<double> mult = 1.0+delta;
  
  std::cout << "Just checking" << std::endl;
  fas->getAmpPtr(comp)->print();
  fas->getAmpPtr(comp)->multiply(mult);
  std::cout << "Just checking" << std::endl;

  DalitzEventPattern pat(421, -211, 211, 211, -211);

  TRandom3 randomGen(seed);

  SignalGenerator signalGenerator(pat, fas.get(), &randomGen);
  DalitzEventList eventList;
  signalGenerator.FillEventList(eventList, nEvents);
  eventList.save( outdir.Data() );

}

  

void compareDandDbarSamples(){
  
  std::map<int, std::vector<int> > componentAmps;
  componentAmps[0 ].push_back(0);
  componentAmps[0 ].push_back(7);
  componentAmps[1 ].push_back(1);
  componentAmps[2 ].push_back(2);
  componentAmps[3 ].push_back(3);
  componentAmps[4 ].push_back(4);
  componentAmps[5 ].push_back(5);
  componentAmps[6 ].push_back(6);
  componentAmps[7 ].push_back(8);
  componentAmps[8 ].push_back(9);
  componentAmps[9 ].push_back(10);
  componentAmps[10].push_back(11);
  componentAmps[11].push_back(12);
  componentAmps[12].push_back(13);
  componentAmps[13].push_back(14);
  componentAmps[13].push_back(15);

  PhilippeModelWrapper* fasPhilWrapper = PhilippeModelWrapper::getStaticWrapper();
  fasPhilWrapper->resetModel();  
  MINT::counted_ptr<FitAmpSum> fas   = fasPhilWrapper->getFitAmpSum  ();    
  
  DalitzEventPattern pat(421, -211, 211, 211, -211);

  
  DalitzEventList eventlistdz;
  DalitzEventList eventlistdzb;
  eventlistdz .fromFile("/Users/sh7566/Documents/Work/Code/FourBodyTools/example/development/toyMC/baseLine_1.root");
  eventlistdzb.fromFile("/Users/sh7566/Documents/Work/Code/FourBodyTools/example/development/toyCPVMC/baseLine_2_9.root");
  

  DalitzEventListWithAmps evtlistwampsDz (eventlistdz );
  DalitzEventListWithAmps evtlistwampsDzb(eventlistdzb);
  
  for (auto i = componentAmps.begin(); i != componentAmps.end(); ++i){
    std::vector<int>& fasIDs = i->second;
    evtlistwampsDz .addAmp(fas.get(), fasIDs);
    evtlistwampsDzb.addAmp(fas.get(), fasIDs);
  }
  evtlistwampsDz .setTotAmp(fas.get() );
  evtlistwampsDzb.setTotAmp(fas.get() );
  
  TRandom3 random(33);

  ModelInspCPAsymmetry cpAys(evtlistwampsDz, evtlistwampsDzb);
  cpAys.createBinningSchemes(100.0, 0.3*0.25, 0.05*0.25);
  cpAys.doPeusdoExp(4000, &random);
  cpAys.makeChiSqPlot("ChiSq");
  cpAys.printChi2Breakdown();
  cpAys.makeAysPlots("Ays");
  

  //ModelInspCPVFitter fitter(evtlistwampsDz, evtlistwampsDzb);
  //
  //double neg2LLHnoCPV = fitter.getVal();
  //
  //std::vector<double> sigs;

  //for (int i = 0; i < 14; i++){
  //  fitter.setComponentToFit(i);
  //  fitter.fit(0);  
  //  double cpvVal = fitter.getVal();
  //  std::cout << "-----------------------------" << std::endl;
  //  std::cout <<   " sig = " << sqrt(neg2LLHnoCPV - cpvVal) << std::endl;
  //  sigs.push_back( sqrt(neg2LLHnoCPV - cpvVal) );
  //  
  //  TString istr = ""; istr += i;
  //  fitter.makePlot("CPVSearch_" + istr);
  //  std::cout << "-----------------------------" << std::endl;
  //}
  //
  //
  //for (int i = 0; i < 14; i++){
  //  std::cout << "-----------------------------" << std::endl;
  //  std::cout << i << " :    sig = " << sigs.at(i) << std::endl;
  //  std::cout << "-----------------------------" << std::endl;
  //}

}


int main(int argc, char** argv) {
  
  
  //const int comp_a1260pToSigma  = 0; 
  //const int comp_a1640mToSigma  = 1; 
  //const int comp_pi1300p        = 2; 
  //const int comp_pi1300m        = 3; 
  //const int comp_pi1670pToSigma = 4; 
  //const int comp_rhoSigma       = 5; 
  //const int comp_f0Sigma        = 6; 
  //const int comp_a1260pToRho    = 7; 
  //const int comp_a1640pToRho    = 8; 
  //const int comp_rhoRhoS        = 9; 
  //const int comp_rhoRhoP        = 10; 
  //const int comp_rhoRhoD        = 11; 
  //const int comp_f2f2           = 12;  
  //const int comp_pi1670pTof2    = 13; 
  //const int comp_a1260mToSigma  = 14; 
  //const int comp_a1260mToRho    = 15; 


  HyperPlotStyle::init();

  std::complex<double> mult(0.5, 0.5);
  int seed = 2;
  int component = 9;
  //generateToyMC(seed);
  
  double thrDeg = (3.0/180.0)*TMath::Pi();
  double fouDeg = (4.0/180.0)*TMath::Pi();
  
  std::complex<double> thrDegRot( cos(thrDeg)-1, sin(thrDeg) );
  std::complex<double> fouDegRot( cos(fouDeg)-1, sin(fouDeg) );
  std::complex<double> fivPerInc( 0.05, 0.0 );
  std::complex<double> fouPerInc( 0.04, 0.0 );

  generateCPVToyMC(seed, comp_a1260pToRho, fivPerInc, "a1-mag"     );
  //generateCPVToyMC(seed, comp_a1260pToRho, thrDegRot, "a1-pha"     );
  //generateCPVToyMC(seed, comp_rhoRhoD    , fivPerInc, "rhorhoD-mag");
  //generateCPVToyMC(seed, comp_rhoRhoD    , fouDegRot, "rhorhoD-pha");  
  //generateCPVToyMC(seed, comp_rhoRhoP    , fouPerInc, "rhorhoP-mag");
  //generateCPVToyMC(seed, comp_rhoRhoP    , thrDegRot, "rhorhoP-pha");  




  //compareDandDbarSamples();


  return 0;

}



