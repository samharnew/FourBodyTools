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
  std::vector<double>                _totChisq;
  std::vector< std::vector<double> > _indChisq;
  WidthFinder                        _totChisqStats;
  std::vector<WidthFinder>           _indChisqStats;


  public: 

  ModelInspCPAsymmetry(DalitzEventListWithAmps& dzEvts, DalitzEventListWithAmps& dzbEvts);
  
  int getNumAmps();

  void createBinningSchemes(double minBinContent, double minPhaseWidth, double minAmpRatioWidth);
  void createBinningSchemes(int ampComponent);

  void fillBinningSchemes(TRandom* random = 0);
  void addDzEvtToBinningSchemes(int evtNum, bool dz = true);
  void addDzbEvtToBinningSchemes(int evtNum, bool dz = true);

  void clearBinningSchemes();
  
  double getChiSq(int component);
  
  void doPeusdoExp(TRandom* random);
  void doPeusdoExp(int nExp, TRandom* random);

  void makeChiSqPlot( std::vector<double> chi2vals, double measChisq, TString outdir, bool incMeas );
  void makeChiSqPlot( TString outdir );

  double getNDOF(int component);
  double getProb(int component);
  double getSig(int component);
  int getNumBins(int component);
  double getVariance(int component);
  double getMean(int component);
  double getScale(int component);

  void printChi2Breakdown();

  void makeAysPlots(int component, TString outdir);
  void makeAysPlots(TString outdir);


  ~ModelInspCPAsymmetry();

};


#endif

