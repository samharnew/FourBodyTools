/**
 * <B>FourBodyTools</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: May 2017
 *  
 **/

 
#ifndef MODEL_INSPIRED_CPV_FITTER_HH
#define MODEL_INSPIRED_CPV_FITTER_HH

#include "DalitzEventListWithAmps.h"

#include "Mint/Minimisable.h"
#include "Mint/Minimiser.h"
#include "Mint/FitParameter.h"
#include "TH2D.h"
#include "HyperHistogram.h"
#include "HyperPlotStyle.h"


class ModelInspCPVFitter : public MINT::Minimisable {
  
  MINT::MinuitParameterSet*          _parSet;

  MINT::Minimiser* _minimiser;

  DalitzEventListWithAmps _dzEvts; 
  DalitzEventListWithAmps _dzbEvts; 
  
  MINT::FitParameter _globAys;
  std::vector< MINT::counted_ptr<MINT::FitParameter> > _realPerts; 
  std::vector< MINT::counted_ptr<MINT::FitParameter> > _imagPerts; 
  
  int _currentComponent;

  public:

  ModelInspCPVFitter(DalitzEventListWithAmps& dzEvts, DalitzEventListWithAmps& dzbEvts);
  
  virtual double getVal();

  double getProbDzNoNorm( int evt, int ampComp );

  double getProbDzbNoNorm( int evt, int ampComp );
  
  void setComponentToFit(int i);

  MINT::Minimiser* getMinimiser();

  void fixComponents(int compNum = -1, bool unfix = false);
  void resetComponents();

  double getNeg2LLHDz (int comp);
  double getNeg2LLHDzb(int comp);

  void makePlot(TString outdir);


  int fit(bool minos = 0);
  
  ~ModelInspCPVFitter();

};


#endif

