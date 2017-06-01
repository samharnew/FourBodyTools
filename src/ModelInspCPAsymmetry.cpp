#include "ModelInspCPAsymmetry.h"




ModelInspCPAsymmetry::ModelInspCPAsymmetry(DalitzEventListWithAmps& dzEvts, DalitzEventListWithAmps& dzbEvts) :
  _dzEvts          ( dzEvts ),
  _dzbEvts         ( dzbEvts),
  _dzBinNums       ( dzEvts .size(), std::vector<int>( _dzEvts.at(0).getNumAmps(), -1 )   ),
  _dzbBinNums      ( dzbEvts.size(), std::vector<int>( _dzEvts.at(0).getNumAmps(), -1 )   ),
  _histsAll        ( _dzEvts.at(0).getNumAmps(), 0 ),
  _histsDz         ( _dzEvts.at(0).getNumAmps(), 0 ),
  _histsDzb        ( _dzEvts.at(0).getNumAmps(), 0 ),
  _minBinContent   ( 100  ),
  _minPhaseWidth   ( 0.1  ),
  _minAmpRatioWidth( 0.05 ),
  _indChisq        ( _dzEvts.at(0).getNumAmps(), std::vector<double>() ),
  _indChisqStats   ( _dzEvts.at(0).getNumAmps(), WidthFinder()         )
{
  
} 

void ModelInspCPAsymmetry::createBinningSchemes(double minBinContent, double minPhaseWidth, double minAmpRatioWidth){
  
  _minBinContent    = minBinContent;
  _minPhaseWidth    = minPhaseWidth;
  _minAmpRatioWidth = minAmpRatioWidth;
  
  for (unsigned i = 0; i < _histsAll.size(); i++){
  	createBinningSchemes(i);
  }

  _totChisq.clear();
  _totChisqStats = WidthFinder();

  _indChisq      = std::vector< std::vector<double> >( getNumAmps(), std::vector<double>() );
  _indChisqStats = std::vector<WidthFinder>( getNumAmps(), WidthFinder() );

}

void ModelInspCPAsymmetry::createBinningSchemes(int component){

  if (_histsAll.at(component) != 0) delete _histsAll.at(component);
  if (_histsDz .at(component) != 0) delete _histsDz .at(component);
  if (_histsDzb.at(component) != 0) delete _histsDzb.at(component);

  HyperPointSet pointsBoth(2);
  HyperPointSet pointsDz  (2);
  HyperPointSet pointsDzb (2);

  for (int i = 0; i < _dzEvts.size(); i++){
  	double phase = arg(_dzEvts.at(i).getAmpRelToTot(component));
  	double ratio = abs(_dzEvts.at(i).getAmpRelToTot(component));
  	HyperPoint point(phase, ratio);
  	pointsBoth.push_back(point);
  	pointsDz  .push_back(point);
  }
  for (int i = 0; i < _dzbEvts.size(); i++){
  	double phase = arg(_dzbEvts.at(i).getAmpRelToTot(component));
  	double ratio = abs(_dzbEvts.at(i).getAmpRelToTot(component));
  	HyperPoint point(phase, ratio);
  	pointsBoth.push_back(point);
  	pointsDzb .push_back(point);
  }
  

  HyperPoint  minWidth( _minPhaseWidth, _minAmpRatioWidth ); 
  HyperCuboid limits( HyperPoint(-TMath::Pi(), 0.0), HyperPoint(+TMath::Pi(), 1.0) );

  HyperHistogram hyperHist(limits, pointsBoth, HyperBinningAlgorithms::SMART , 
  	                                           AlgOption::MinBinWidth(minWidth), 
  	                                           AlgOption::MinBinContent(_minBinContent) );

  _histsAll.at(component) = new HyperHistogram(hyperHist.getBinning());
  _histsDz .at(component) = new HyperHistogram(hyperHist.getBinning());
  _histsDzb.at(component) = new HyperHistogram(hyperHist.getBinning());

  for (int i = 0; i < _dzEvts.size(); i++){
  	int binNum = hyperHist.getBinning().getBinNum( pointsDz.at(i) );
  	_dzBinNums.at(i).at(component) = binNum;
  	_histsAll.at(component)->fill ( pointsDz.at(i) );
  }  

  for (int i = 0; i < _dzEvts.size(); i++){
  	int binNum = hyperHist.getBinning().getBinNum( pointsDzb.at(i) );
  	_dzbBinNums.at(i).at(component) = binNum;
  	_histsAll.at(component)->fill ( pointsDzb.at(i) );
  }  

}

void ModelInspCPAsymmetry::clearBinningSchemes(){
 
  for (unsigned i = 0; i < _histsAll.size(); i++){
    //_histsAll.at(i)->clear();
    _histsDz .at(i)->clear();
    _histsDzb.at(i)->clear();
  }

}


void ModelInspCPAsymmetry::fillBinningSchemes(TRandom* random){

  clearBinningSchemes();

  for (int i = 0; i < _dzEvts.size(); i++){
    bool dz = true;
    if (random != 0){
    	if (random->Uniform(0.0,1.0) > 0.5) dz = false;
    }
  	addDzEvtToBinningSchemes(i, dz);
  }
  for (int i = 0; i < _dzbEvts.size(); i++){
    bool dz = false;
    if (random != 0){
    	if (random->Uniform(0.0,1.0) > 0.5) dz = true;
    }
  	addDzbEvtToBinningSchemes(i, dz);
  }


}

void ModelInspCPAsymmetry::addDzEvtToBinningSchemes(int evtNum, bool dz){

  for (int i = 0; i < getNumAmps(); i++){
    int binNum = _dzBinNums.at(evtNum).at(i);
    if (dz) _histsDz .at(i)->fillBase(binNum, 1.0);
    else    _histsDzb.at(i)->fillBase(binNum, 1.0);
  }

}

void ModelInspCPAsymmetry::addDzbEvtToBinningSchemes(int evtNum, bool dz){

  for (int i = 0; i < getNumAmps(); i++){
    int binNum = _dzbBinNums.at(evtNum).at(i);
    if (dz) _histsDz .at(i)->fillBase(binNum, 1.0);
    else    _histsDzb.at(i)->fillBase(binNum, 1.0);
  }

}

double ModelInspCPAsymmetry::getChiSq(int component){
  
  if (component != -1){
    HyperHistogram* dzHist  = _histsDz .at(component);
    HyperHistogram* dzbHist = _histsDzb.at(component);
    return dzHist->chi2(*dzbHist);
  }

  double chisq = 0.0;

  for (int i = 0; i < getNumAmps(); i++){
    HyperHistogram* dzHist  = _histsDz .at(i);
    HyperHistogram* dzbHist = _histsDzb.at(i);
    chisq += dzHist->chi2(*dzbHist);
  }

  return chisq;
}


int ModelInspCPAsymmetry::getNumAmps(){
  return _histsDz.size();
}


void ModelInspCPAsymmetry::doPeusdoExp(TRandom* random){

  fillBinningSchemes(random);

  std::vector<double> chisqvals;

  double totChi2 = 0.0;
  for (int i = 0; i < getNumAmps(); i++){
  	double chi2 = getChiSq(i);

  	_indChisq     .at(i).push_back(chi2);
  	_indChisqStats.at(i).add      (chi2);

  	totChi2 += chi2;
  }

  _totChisq     .push_back( totChi2 );
  _totChisqStats.add      ( totChi2 );

}

void ModelInspCPAsymmetry::doPeusdoExp(int nExp, TRandom* random){
  
  for (int i = 0; i < nExp; i++){
  	if (i % 100 == 0) std::cout << "Finished toy " << i << " of " << nExp << std::endl;
  	doPeusdoExp(random);
  }

}





void ModelInspCPAsymmetry::makeAysPlots(int component, TString outdir){

  HyperHistogram* dzHist  = _histsDz .at(component);
  HyperHistogram* dzbHist = _histsDzb.at(component);  
  
  dzHist ->setNames( HyperName("#phi [radians]", "r") );
  dzbHist->setNames( HyperName("#phi [radians]", "r") );

  HyperPlotStyle::setPalette("birdy");
  
  double maxDen1 = dzHist ->getMaxDensity();
  double maxDen2 = dzbHist->getMaxDensity();
  double maxDen  = maxDen1>maxDen2?maxDen1:maxDen2;
  
  dzHist ->setMaxDensity(maxDen);
  dzbHist->setMaxDensity(maxDen);
  dzHist ->setMinDensity(0.0   );
  dzbHist->setMinDensity(0.0   );

  dzHist ->drawDensity(outdir + "_Dz" , "COLZ Edges1");
  dzbHist->drawDensity(outdir + "_Dzb", "COLZ Edges1");

  HyperPlotStyle::setPalette("pulls");
  
  HyperHistogram hist(dzHist->getBinning());
  hist.setNames( HyperName("#phi [radians]", "r") );
  hist.asymmetry   (*dzHist, *dzbHist);
  hist.setMin(-1.0);
  hist.setMax( 1.0);
  hist.draw(outdir + "_Ays", "COLZ Edges1");

  //gStyle->SetPalette(kCoffee);
  //TColor::InvertPalette();

  HyperHistogram histSig(dzHist->getBinning());
  histSig.setNames( HyperName("#phi [radians]", "r") );
  histSig.pulls   (*dzHist, *dzbHist);
  double max = fabs(histSig.getMax());
  double min = fabs(histSig.getMin());
  double minmax = max>min?max:min;
  histSig.setMin(-minmax);
  histSig.setMax( minmax);
  histSig.draw(outdir + "_Pull", "COLZ Edges1");
 
  HyperPlotStyle::setPalette("birdy");

}


void ModelInspCPAsymmetry::makeAysPlots(TString outdir){

  for (int i = 0; i < getNumAmps(); i++){
  	TString istr = ""; istr += i;
    makeAysPlots( i , outdir + "_amp" + istr );
  }  

}


void ModelInspCPAsymmetry::makeChiSqPlot(std::vector<double> chi2vals, double measChisq, TString outdir, bool incMeas){
  

  WidthFinder stats;
  
  int ntoys = chi2vals.size();

  for (int i = 0; i < ntoys; i++){
    stats.add( chi2vals.at(i) );
  }

  WidthFinder stats2(stats);
  if (incMeas) stats2.add(measChisq);
  
  int nbins = 100;
  double min = stats2.getMin() - stats2.range()*0.05;
  double max = stats2.getMax() + stats2.range()*0.05;
  double binWid = (max - min)/double(nbins);

  TH1D hist( "chisq", "chisq", nbins, min, max );
  for (int i = 0; i < ntoys; i++){
    hist.Fill( chi2vals.at(i) );
  }
  
  double var  = stats.varience();
  double mean = stats.mean    ();

  double scale = 0.5*(var/mean);
  double ndof  = mean/scale;

  TF1* chi2_naive = new TF1("chi2Dist","[0]*ROOT::Math::chisquared_pdf(x,[1],0)",min, max);
  chi2_naive->SetParameters( ntoys*binWid, mean );
  chi2_naive->SetLineColor(kBlue);

  TF1* chi2_scaled = new TF1("chi2Dist","[0]*ROOT::Math::chisquared_pdf(x/[2],[1],0)",min, max);
  chi2_scaled->SetParameters( (ntoys*binWid)/scale, ndof, scale);
  chi2_scaled->SetLineColor(kRed);
  
  TString strNdofNaive  = ""; strNdofNaive  += mean;
  TString strNdofScaled = ""; strNdofScaled += ndof;
  TString strScale      = ""; strScale      += scale;
  
  strNdofNaive  = strNdofNaive (0, ceil(log10(mean ))+2);
  strNdofScaled = strNdofScaled(0, ceil(log10(ndof ))+2);
  strScale      = strScale     (0, ceil(log10(scale))+3);

  RootPlotter1D plotter(&hist);
  plotter.addObject(chi2_naive );
  plotter.addObject(chi2_scaled);
  plotter.addText( "ndof  = " + strNdofNaive , 0.63, 0.79    , 1, 2, 0.065, true, kBlue );
  plotter.addText( "ndof  = " + strNdofScaled, 0.63, 0.79-0.1, 1, 2, 0.065, true, kRed  );
  plotter.addText( "scale = " + strScale     , 0.63, 0.79-0.2, 1, 2, 0.065, true, kRed  );

  if (incMeas) plotter.addVerticalLine(measChisq, 1, kRed );
  plotter.plot(outdir);
  plotter.setMin(0.5);
  plotter.logY();
  plotter.plot(outdir + "_log");


}

void ModelInspCPAsymmetry::makeChiSqPlot( TString outdir ){
  
  std::vector<double> indChi2;

  fillBinningSchemes();

  double totChi2 = 0.0;

  for (int i = 0; i < getNumAmps(); i++){
  	TString istr = ""; istr += i;
    double chisq = getChiSq(i);
    makeChiSqPlot( _indChisq.at(i), chisq , outdir + "_amp" + istr + "_incMeas", true  );
    makeChiSqPlot( _indChisq.at(i), chisq , outdir + "_amp" + istr , false );

    totChi2 += chisq;
  }

  makeChiSqPlot( _totChisq, totChi2 , outdir + "_incMeas", true   );
  makeChiSqPlot( _totChisq, totChi2 , outdir, false  );


}

double ModelInspCPAsymmetry::getNDOF(int component){
  double mean = getMean(component);
  double var  = getVariance(component);
  double scale = 0.5*(var/mean);
  double ndof  = mean/scale;  
  return ndof;
}

double ModelInspCPAsymmetry::getScale(int component){
  double mean = getMean(component);
  double var  = getVariance(component);
  double scale = 0.5*(var/mean);
  return scale;
}

double ModelInspCPAsymmetry::getVariance(int component){
  if (component == -1) return _totChisqStats      .varience();
  return _indChisqStats.at(component).varience();
}

double ModelInspCPAsymmetry::getMean(int component){
  if (component == -1) return _totChisqStats      .mean();
  return _indChisqStats.at(component).mean();
}

double ModelInspCPAsymmetry::getProb(int component){
  double ndof  = getNDOF  (component);
  double scale = getScale (component);
  double chi2  = getChiSq (component);
  double prob  = TMath::Prob(chi2/scale, floor(ndof) );
  return prob;
}

double ModelInspCPAsymmetry::getSig(int component){
  double prob = getProb(component);
  return fabs(TMath::NormQuantile( prob/2.0 ));
}

int ModelInspCPAsymmetry::getNumBins(int component){
  if (component != -1) return _histsAll.at(component)->getNBins();
  
  int nbins = 0;
  for (int i = 0; i < getNumAmps(); i++){
  	nbins += _histsAll.at(i)->getNBins();
  }
  return nbins;

}


void ModelInspCPAsymmetry::printChi2Breakdown(){

  std::vector<double> indChi2;

  fillBinningSchemes();
  
  std::cout << "--------------- ALL HISTS ----------------" << std::endl;
  std::cout << "chi2  = " << getChiSq(-1) << std::endl; 
  std::cout << "ndof  = " << getNDOF(-1) << std::endl;
  std::cout << "scale = " << getScale(-1) << std::endl;
  std::cout << "pval  = " << getProb(-1) << std::endl;
  std::cout << "sigma = " << getSig (-1) << std::endl;
  std::cout << "nbins = " << getNumBins (-1) << std::endl;

  for (int i = 0; i < getNumAmps(); i++){
    std::cout << "------------- HIST " << i << " ----------------" << std::endl;
    std::cout << "chi2  = " << getChiSq(i) << std::endl; 
    std::cout << "ndof  = " << getNDOF (i) << std::endl;
    std::cout << "scale = " << getScale(i) << std::endl;    
    std::cout << "pval  = " << getProb(i) << std::endl;
    std::cout << "sigma = " << getSig (i) << std::endl;
    std::cout << "nbins = " << getNumBins (i) << std::endl;
  }



}

ModelInspCPAsymmetry::~ModelInspCPAsymmetry(){


}


