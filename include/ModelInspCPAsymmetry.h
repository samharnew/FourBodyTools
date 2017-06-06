/**
 * <B>FourBodyTools</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: May 2017
 *  
 **/

 
#ifndef MODEL_INSPIRED_AYS_HH
#define MODEL_INSPIRED_AYS_HH

#include "DalitzEventListWithAmps.h"

#include "Mint/Minimisable.h"
#include "Mint/Minimiser.h"
#include "Mint/FitParameter.h"
#include "TH2D.h"
#include "HyperHistogram.h"
#include "HyperPlotStyle.h"
#include "StatisticsFinder.h"
#include "TF1.h"

#include "Riostream.h"
#include "TROOT.h"
#include "TColor.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include "TVirtualPad.h"
#include "TVirtualX.h"
#include "TError.h"
#include "TMathBase.h"
#include "TApplication.h"


class CPAsyChi2Res {

  double              _totChisq;
  std::vector<double> _indChisq;   
  
  public:
  
  CPAsyChi2Res(int nAmps);

  void setChiSq(int amp, double chi2);
  
  void setTotChiSq(double chi2 = -1);
  double getChiSq(int amp);
  double getTotChiSq();

  int getNumAmps();
  
  ~CPAsyChi2Res();

};

class CPAsyChi2ResSet {

  std::vector<CPAsyChi2Res> _toyResults;
  CPAsyChi2Res _dataResults;  
  
  std::vector<WidthFinder> _toyStatsInd;
  WidthFinder _toyStatsTot;

  public:

  CPAsyChi2ResSet(int nAmps);
  CPAsyChi2ResSet(TString filename);

  int getNumAmps();
  
  void setDataChiSq(CPAsyChi2Res results);

  void addToyChiSq(CPAsyChi2Res results);
  
  double getToyChiSq(int component, int toy);

  WidthFinder& getToyStats(int component);


  double getDataChiSq(int component);
  
  int getNumToys();
    
  double getNDOF(int component);
  
  double getScale(int component);
  
  double getVariance(int component);
  
  double getMean(int component);
  double getProb(int component);
  double getSig(int component);


  void makeChiSqPlot(int component, TString outdir, bool incMeas);
  void makeChiSqPlot( TString outdir );

  void saveResults(TString outdir);
  void loadResults(TString outdir);

  void print();


  ~CPAsyChi2ResSet();


};


class ModelInspCPAsymmetry {
  
  DalitzEventListWithAmps _dzEvts ; 
  DalitzEventListWithAmps _dzbEvts; 
  
  //This stores the bin numbers for the dz and dzb samples for each binning scheme
  std::vector< std::vector< int > > _dzBinNums;
  std::vector< std::vector< int > > _dzbBinNums;
  
  //Adaptive binning hsitograms, one for each resonance
  std::vector<HyperHistogram*> _histsAll;
  std::vector<HyperHistogram*> _histsDz ;
  std::vector<HyperHistogram*> _histsDzb;
  
  //Adaptive binning properties
  double _minBinContent;
  double _minPhaseWidth;
  double _minAmpRatioWidth;

  //save the chi2 values that are found from peusdo experiements
  CPAsyChi2ResSet _chiSqRes;
  
  //save some limits for the 2D plotting so they can be common among 
  //all plots
  double _maxAys ;
  double _maxPull;

  public: 

  ModelInspCPAsymmetry(DalitzEventListWithAmps& dzEvts, DalitzEventListWithAmps& dzbEvts);
  
  int getNumAmps();

  void createBinningSchemes(double minBinContent, double minPhaseWidth, double minAmpRatioWidth);
  void createBinningSchemes(int ampComponent);

  void fillBinningSchemes(TRandom* random = 0);
  void addDzEvtToBinningSchemes(int evtNum, bool dz = true);
  void addDzbEvtToBinningSchemes(int evtNum, bool dz = true);
  
  void normaliseDzbHistsToDz();

  void clearBinningSchemes();
  
  double getChiSq(int component);
  CPAsyChi2Res getChiSqContainer();
  

  void doPeusdoExp(TRandom* random);
  void doPeusdoExp(int nExp, TRandom* random);

  void makeChiSqPlot( TString outdir );

  int getNumBins(int component);


  HyperHistogram getAysHist (int component);
  HyperHistogram getPullHist(int component);
  
  double getMaxAbsAys(int component);
  double getMaxAbsPull(int component);
  
  double getMaxAbsAys();
  double getMaxAbsPull();

  void updateHistLimits();

  void saveResults(TString outdir);
  void printChi2Breakdown();

  void makeAysPlots(int component, TString outdir);
  void makeAysPlots(TString outdir);


  ~ModelInspCPAsymmetry();

};







#endif

