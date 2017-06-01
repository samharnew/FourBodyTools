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

TString g_outdir = "";


TString GetToyMCDir(TString modelName, TString uniqueTag, int seed ){

  TString outdir = g_outdir + "toyMC/";
  gSystem->Exec("mkdir -p " + outdir);
  outdir += modelName;
  outdir += "_";
  outdir += uniqueTag;
  if (seed != -1){
    outdir += "_seed-";
    outdir += seed;
  }
  outdir += ".root";
  
  return outdir;
}

void GenerateToyMC(int seed, int comp, std::complex<double> delta, TString uniqueTag ){
  
  int nEvents = 100000;
  
  //Get the model
  PhilippeModelWrapper* fasPhilWrapper = PhilippeModelWrapper::getStaticWrapper();
  fasPhilWrapper->resetModel();
  MINT::counted_ptr<FitAmpSum> fas   = fasPhilWrapper->getFitAmpSum  ();    
  

  TString outdir = GetToyMCDir( fasPhilWrapper->getModelName(), uniqueTag, seed );
  
  //mulitply the given component to inject CPV
  std::complex<double> mult = 1.0+delta;
  fas->getAmpPtr(comp)->print();
  fas->getAmpPtr(comp)->multiply(mult);
  
  //Generate according to this model and save

  DalitzEventPattern pat(421, -211, 211, 211, -211);

  TRandom3 randomGen(seed);

  SignalGenerator signalGenerator(pat, fas.get(), &randomGen);
  DalitzEventList eventList;
  signalGenerator.FillEventList(eventList, nEvents);
  eventList.save( outdir.Data() );

}

  

void CPVTest(TString uniqueTag1, TString uniqueTag2){
  
  PhilippeModelWrapper* fasPhilWrapper = PhilippeModelWrapper::getStaticWrapper();
  fasPhilWrapper->resetModel();  
  MINT::counted_ptr<FitAmpSum> fas   = fasPhilWrapper->getFitAmpSum  ();    
  
  DalitzEventPattern pat(421, -211, 211, 211, -211);

  DalitzEventList eventlistdz;
  DalitzEventList eventlistdzb;
  eventlistdz .fromFile( GetToyMCDir("baseLine", uniqueTag1, -1).Data() );
  eventlistdzb.fromFile( GetToyMCDir("baseLine", uniqueTag2, -1).Data() );
  
  DalitzEventListWithAmps evtlistwampsDz (eventlistdz );
  DalitzEventListWithAmps evtlistwampsDzb(eventlistdzb);
  
  for (unsigned i = 0; i < fas->size(); i++){
    evtlistwampsDz .addAmp(fas.get(), i );
    evtlistwampsDzb.addAmp(fas.get(), i );
  }
  evtlistwampsDz .setTotAmp(fas.get() );
  evtlistwampsDzb.setTotAmp(fas.get() );
  
  TRandom3 random(33);
  
  TString outdir = g_outdir;
  outdir += "cpvTest/";
  gSystem->Exec("mkdir -p " + outdir);
  outdir += uniqueTag1 + "_" + uniqueTag2 + "/";
  gSystem->Exec("mkdir -p " + outdir);

  ModelInspCPAsymmetry cpAys(evtlistwampsDz, evtlistwampsDzb);
  cpAys.createBinningSchemes(100.0, 0.3, 0.05);
  cpAys.doPeusdoExp(4000, &random);
  cpAys.makeChiSqPlot(outdir + "ChiSq");
  cpAys.printChi2Breakdown();
  cpAys.makeAysPlots(outdir + "Ays");
  

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
  
  HyperPlotStyle::init();  

  bool generate    =  0;
  bool cpvTest     =  0;

  int cpvOption    =  0; 
  int seed         =  1; 
  
   

  for(int i = 1; i<argc; i=i+2){
  
    //Options to do with offline selection
    if       (std::string(argv[i])=="--generate"       ) { generate    =  1         ; i--; }
    else if  (std::string(argv[i])=="--cpv-test"       ) { cpvTest     =  1         ; i--; }
    else if  (std::string(argv[i])=="--cpv-opt"        ) { cpvOption   =  atoi(argv[i+1]); }
    else if  (std::string(argv[i])=="--seed"           ) { seed        =  atoi(argv[i+1]); }
    else if  (std::string(argv[i])=="--outdir"         ) { g_outdir    =  argv[i+1];       }
    else { 
      std::cout << "Entered invalid argument " << argv[i] << std::endl;
      return 0;
    }
  }
  
  //Some constants that are useful later
  const double thrDeg = (3.0/180.0)*TMath::Pi();
  const double fouDeg = (4.0/180.0)*TMath::Pi();  
  const std::complex<double> thrDegRot( cos(thrDeg)-1, sin(thrDeg) );
  const std::complex<double> fouDegRot( cos(fouDeg)-1, sin(fouDeg) );
  const std::complex<double> fivPerInc( 0.05, 0.0 );
  const std::complex<double> fouPerInc( 0.04, 0.0 );

  //What each cpvOption means 
  std::complex<double> delta(0.0, 0.0);
  int component    = 0;
  TString sampleID = "nominal";

  if (cpvOption == 1){
    delta     = fivPerInc;
    component = comp_a1260pToRho;
    sampleID  = "a1-mag";
  }
  if (cpvOption == 2){
    delta     = thrDegRot;
    component = comp_a1260pToRho;
    sampleID  = "a1-pha";
  }
  if (cpvOption == 3){
    delta     = fivPerInc;
    component = comp_rhoRhoD;
    sampleID  = "rhorhoD-mag";
  }
  if (cpvOption == 4){
    delta     = fouDegRot;
    component = comp_rhoRhoD;
    sampleID  = "rhorhoD-pha";
  }
  if (cpvOption == 5){
    delta     = fouPerInc;
    component = comp_rhoRhoP;
    sampleID  = "rhorhoP-mag";
  }
  if (cpvOption == 6){
    delta     = thrDegRot;
    component = comp_rhoRhoP;
    sampleID  = "rhorhoP-pha";
  }

  if (generate) GenerateToyMC(seed, component, delta, sampleID );
  if (cpvTest ) CPVTest("nominal", sampleID);



  return 0;

}



