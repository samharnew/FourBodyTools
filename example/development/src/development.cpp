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
  
  int modSeed = comp*1000 + seed;

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

  TRandom3 randomGen(modSeed);

  SignalGenerator signalGenerator(pat, fas.get(), &randomGen);
  DalitzEventList eventList;
  signalGenerator.FillEventList(eventList, nEvents);
  eventList.save( outdir.Data() );

}

void HAddSamples(int lowseed, int highseed, TString sampleID){
  
  TString outputFile = GetToyMCDir("baseLine", sampleID, -1);

  TString fileListToMerge = "";
  for (int i = lowseed; i <= highseed; i++){
    fileListToMerge += GetToyMCDir("baseLine", sampleID, i);
    fileListToMerge += " ";
  }
  
  gSystem->Exec("hadd -f " + outputFile + " " + fileListToMerge);

}  

void CPVTest(TString uniqueTag1, TString uniqueTag2, int ntoys, int seed ){
  
  PhilippeModelWrapper* fasPhilWrapper = PhilippeModelWrapper::getStaticWrapper();
  fasPhilWrapper->resetModel();  
  MINT::counted_ptr<FitAmpSum> fas   = fasPhilWrapper->getFitAmpSum  ();    
  
  TString modelName = fasPhilWrapper->getModelName();

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
  
  TRandom3 random(seed);
  
  TString outdir = g_outdir;
  outdir += "cpvTest/";
  gSystem->Exec("mkdir -p " + outdir);
  outdir += uniqueTag1 + "_" + uniqueTag2 + "/";
  gSystem->Exec("mkdir -p " + outdir);
  outdir += modelName + "/";
  gSystem->Exec("mkdir -p " + outdir);

  ModelInspCPAsymmetry cpAys(evtlistwampsDz, evtlistwampsDzb);
  cpAys.createBinningSchemes(100.0, 0.35, 0.1);
  cpAys.doPeusdoExp(ntoys, &random);
  cpAys.makeChiSqPlot(outdir + "ChiSq");
  cpAys.printChi2Breakdown();
  cpAys.makeAysPlots(outdir + "Ays");
  cpAys.saveResults(outdir + "results");
  

  CPAsyChi2ResSet results(outdir + "results");
  results.print();

}

inline bool fileExists(const std::string& name) {
  std::ifstream f(name.c_str());
  return f.good();
}

void ResultsTable( std::ostream& os, TString uniqueTag1, TString uniqueTag2, TString genModel, std::vector<TString> altModels, bool header, bool footer ){
  
  TString outdir = g_outdir;
  outdir += "cpvTest/";
  outdir += uniqueTag1 + "_" + uniqueTag2 + "/";

  CPAsyChi2ResSet results(outdir + genModel + "/results");
  
  double pvalue = results.getProb(-1);
  
  WidthFinder stats;
  for (unsigned i = 0; i < altModels.size(); i++){
    TString filename = outdir + altModels.at(i) + "/results";
    if ( fileExists(filename.Data()) == false ){
      std::cout << "Could not find a results file at " << filename << std::endl;
      continue;
    }
    CPAsyChi2ResSet resultsAlt(filename);
    stats.add( resultsAlt.getProb(-1) );
  }

  if (header){
    os << "\\begin{tabular}{ l | l | l }" << std::endl;
    os << "\\hline" << std::endl;
    os << "CPV Sample & p-value (nominal model) & p-value range (alt. models) \\\\" << std::endl;
    os << "\\hline" << std::endl;
  }

  os << std::setw(20) << std::left << uniqueTag2 << " & " 
     << std::setw(20) << std::left << pvalue << " & " 
     << std::setw(20) << std::left << stats.getMin() << " , " << std::left << stats.getMax() << "\\\\" << std::endl;

  if (footer){
    os << "\\end{tabular}" << std::endl;    
  }

}

void ResultsTable( ){
  
  std::ofstream myfiletex;
  myfiletex.open (g_outdir + "ResultsTable.tex");   

  TString nomModel = "baseLine";

  std::vector<TString> altModels;
  altModels.push_back("altModel1");
  altModels.push_back("altModel2");
  altModels.push_back("altModel3");
  altModels.push_back("altModel4");
  altModels.push_back("altModel5");
  altModels.push_back("extended" );
  altModels.push_back("noA11640" );
  altModels.push_back("noPi1300" );
  
  std::vector<TString> cpvSamples;
  cpvSamples.push_back("a1-mag" );
  cpvSamples.push_back("a1-pha" );
  cpvSamples.push_back("rhorhoD-mag" );
  cpvSamples.push_back("rhorhoD-pha" );
  cpvSamples.push_back("rhorhoP-mag" );
  cpvSamples.push_back("rhorhoP-pha" );

  int nsamples = cpvSamples.size();
  
  for (int i = 0; i < nsamples; i++){
    ResultsTable(myfiletex, "nominal", cpvSamples.at(i), nomModel, altModels, i==0, i==nsamples-1);
  }
  
  myfiletex.close();

}

int main(int argc, char** argv) {
  
  HyperPlotStyle::init();  

  bool generate    =  0;
  bool cpvTest     =  0;
  bool hadd        =  0;
  bool resultsTab  =  0;

  int cpvOption    =  0; 
  int seed         =  1; 
  int ntoys        =  5000;
   

  for(int i = 1; i<argc; i=i+2){
  
    //Options to do with offline selection
    if       (std::string(argv[i])=="--generate"       ) { generate    =  1         ; i--; }
    else if  (std::string(argv[i])=="--cpv-test"       ) { cpvTest     =  1         ; i--; }
    else if  (std::string(argv[i])=="--hadd"           ) { hadd        =  1         ; i--; }
    else if  (std::string(argv[i])=="--cpv-opt"        ) { cpvOption   =  atoi(argv[i+1]); }
    else if  (std::string(argv[i])=="--seed"           ) { seed        =  atoi(argv[i+1]); }
    else if  (std::string(argv[i])=="--outdir"         ) { g_outdir    =  argv[i+1];       }
    else if  (std::string(argv[i])=="--ntoys"          ) { ntoys       =  atoi(argv[i+1]);       }
    else if  (std::string(argv[i])=="--results-table"  ) { resultsTab  =  atoi(argv[i+1]);       }

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

  if (generate  ) GenerateToyMC(seed, component, delta, sampleID );
  if (cpvTest   ) CPVTest("nominal", sampleID, ntoys, seed);
  if (hadd      ) HAddSamples(1, seed, sampleID);
  if (resultsTab) ResultsTable();


  return 0;

}



