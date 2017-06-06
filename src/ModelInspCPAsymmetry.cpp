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
  _chiSqRes        ( _dzEvts.at(0).getNumAmps() )
{
  
} 

void ModelInspCPAsymmetry::createBinningSchemes(double minBinContent, double minPhaseWidth, double minAmpRatioWidth){
  
  _minBinContent    = minBinContent;
  _minPhaseWidth    = minPhaseWidth;
  _minAmpRatioWidth = minAmpRatioWidth;
  
  for (unsigned i = 0; i < _histsAll.size(); i++){
    createBinningSchemes(i);
  }
  
  fillBinningSchemes(0);

  _chiSqRes      = CPAsyChi2ResSet( getNumAmps() ); 
  _chiSqRes.setDataChiSq( getChiSqContainer() );
  
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
    if (ratio == 0.0) {
      std::cout << "Error: The ratio is zero, and I want log(ratio)" << std::endl;
      continue;
    }
    HyperPoint point(phase, log10(ratio) );
    pointsBoth.push_back(point);
    pointsDz  .push_back(point);
  }
  for (int i = 0; i < _dzbEvts.size(); i++){
    double phase = arg(_dzbEvts.at(i).getAmpRelToTot(component));
    double ratio = abs(_dzbEvts.at(i).getAmpRelToTot(component));
    if (ratio == 0.0) {
      std::cout << "Error: The ratio is zero, and I want log(ratio)" << std::endl;
      continue;
    }
    HyperPoint point(phase, log10(ratio) );
    pointsBoth.push_back(point);
    pointsDzb .push_back(point);
  }
  
  HyperCuboid limits( HyperPoint(-TMath::Pi(), pointsBoth.getMin().at(1) - 0.1 ), HyperPoint(+TMath::Pi(), pointsBoth.getMax().at(1) + 0.1 ) );

  HyperPoint  minWidth( _minPhaseWidth, _minAmpRatioWidth ); 
  std::vector<int> binningDim; binningDim.push_back(0);

  HyperHistogram hyperHistA(limits, pointsBoth, HyperBinningAlgorithms::MINT , 
                                               AlgOption::MinBinWidth(TMath::Pi()*0.95), 
                                               AlgOption::MinBinContent(_minBinContent),
                                               AlgOption::BinningDimensions(binningDim)
                                                );
  
  HyperHistogram hyperHist(limits, pointsBoth, HyperBinningAlgorithms::MINT , 
                                               AlgOption::MinBinWidth(minWidth), 
                                               AlgOption::MinBinContent(_minBinContent),
                                               AlgOption::StartBinning((HyperBinning&)hyperHistA.getBinning())
                                               //AlgOption::SnapToGrid(true),
                                               //AlgOption::GridMultiplier(2)
                                                );

  _histsAll.at(component) = new HyperHistogram(hyperHist.getBinning());
  _histsDz .at(component) = new HyperHistogram(hyperHist.getBinning());
  _histsDzb.at(component) = new HyperHistogram(hyperHist.getBinning());

  for (int i = 0; i < _dzEvts.size(); i++){
    int binNum = hyperHist.getBinning().getBinNum( pointsDz.at(i) );
    _dzBinNums.at(i).at(component) = binNum;
    _histsAll.at(component)->fill ( pointsDz.at(i) );
  }  

  for (int i = 0; i < _dzbEvts.size(); i++){
    int binNum = hyperHist.getBinning().getBinNum( pointsDzb.at(i) );
    _dzbBinNums.at(i).at(component) = binNum;
    _histsAll.at(component)->fill ( pointsDzb.at(i) );
  }  


}

void ModelInspCPAsymmetry::clearBinningSchemes(){
 
  for (unsigned i = 0; i < _histsAll.size(); i++){
    _histsDz .at(i)->clear();
    _histsDzb.at(i)->clear();
  }

}


void ModelInspCPAsymmetry::fillBinningSchemes(TRandom* random){

  clearBinningSchemes();

  for (int i = 0; i < _dzEvts.size(); i++){
    bool dz = true;
    if (random != 0){
      if (random->Rndm() > 0.5) dz = false;
    }
    addDzEvtToBinningSchemes(i, dz);
  }
  for (int i = 0; i < _dzbEvts.size(); i++){
    bool dz = false;
    if (random != 0){
      if (random->Rndm() > 0.5) dz = true;
    }
    addDzbEvtToBinningSchemes(i, dz);
  }
  
  normaliseDzbHistsToDz();

}

void ModelInspCPAsymmetry::normaliseDzbHistsToDz(){
  
  for (int i = 0; i < getNumAmps(); i++){
    _histsDzb .at(i)->normalise( _histsDz.at(i)->integral() );
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

CPAsyChi2Res ModelInspCPAsymmetry::getChiSqContainer(){

  CPAsyChi2Res results(getNumAmps());
  
  for (int i = 0; i < getNumAmps(); i++){
    results.setChiSq(i, getChiSq(i));
  }
  results.setTotChiSq();

  return results;
  
}

int ModelInspCPAsymmetry::getNumAmps(){
  return _histsDz.size();
}




void ModelInspCPAsymmetry::doPeusdoExp(TRandom* random){

  fillBinningSchemes(random);

  _chiSqRes.addToyChiSq( getChiSqContainer() );

}

void ModelInspCPAsymmetry::doPeusdoExp(int nExp, TRandom* random){
  
  for (int i = 0; i < nExp; i++){
    if (i % 100 == 0) std::cout << "Finished toy " << i << " of " << nExp << std::endl;
    doPeusdoExp(random);
  }

}


HyperHistogram ModelInspCPAsymmetry::getAysHist(int component){

  HyperHistogram* dzHist  = _histsDz .at(component);
  HyperHistogram* dzbHist = _histsDzb.at(component);  
  
  HyperHistogram hist(dzHist->getBinning());
  hist.setNames( HyperName("#phi [radians]", "r") );
  hist.asymmetry   (*dzHist, *dzbHist);
  
  return hist;

}

HyperHistogram ModelInspCPAsymmetry::getPullHist(int component){

  HyperHistogram* dzHist  = _histsDz .at(component);
  HyperHistogram* dzbHist = _histsDzb.at(component);  
  
  HyperHistogram hist(dzHist->getBinning());
  hist.setNames( HyperName("#phi [radians]", "r") );
  hist.pulls   (*dzHist, *dzbHist);
  
  return hist;

}

double ModelInspCPAsymmetry::getMaxAbsAys(int component){
  
  HyperHistogram hist = getAysHist(component);
  double min = fabs(hist.getMin());
  double max = fabs(hist.getMax());
  return max>min?max:min;

}

double ModelInspCPAsymmetry::getMaxAbsPull(int component){
  
  HyperHistogram hist = getPullHist(component);
  double min = fabs(hist.getMin());
  double max = fabs(hist.getMax());
  return max>min?max:min;

}

double ModelInspCPAsymmetry::getMaxAbsAys(){
  
  MinMaxFinder stats;
  for (int i = 0; i < getNumAmps(); i++){
    stats.add( getMaxAbsAys(i) );
  }
  return stats.getMax();

}

double ModelInspCPAsymmetry::getMaxAbsPull(){
  
  MinMaxFinder stats;
  for (int i = 0; i < getNumAmps(); i++){
    stats.add( getMaxAbsPull(i) );
  }
  return stats.getMax();

}

void ModelInspCPAsymmetry::updateHistLimits(){

  _maxAys  = getMaxAbsAys ();
  _maxPull = getMaxAbsPull();

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
  hist.setMin(-_maxAys);
  hist.setMax( _maxAys);
  hist.draw(outdir + "_Ays", "COLZ Edges1");

  //gStyle->SetPalette(kCoffee);
  //TColor::InvertPalette();

  HyperHistogram histSig(dzHist->getBinning());
  histSig.setNames( HyperName("#phi [radians]", "r") );
  histSig.pulls   (*dzHist, *dzbHist);
  double max = fabs(histSig.getMax());
  double min = fabs(histSig.getMin());
  double minmax = max>min?max:min;
  histSig.setMin(-_maxPull);
  histSig.setMax( _maxPull);
  histSig.draw(outdir + "_Pull", "COLZ Edges1");
 
  HyperPlotStyle::setPalette("birdy");

}


void ModelInspCPAsymmetry::makeAysPlots(TString outdir){
  
  updateHistLimits();

  for (int i = 0; i < getNumAmps(); i++){
    TString istr = ""; istr += i;
    makeAysPlots( i , outdir + "_amp" + istr );
  }  

}




void ModelInspCPAsymmetry::makeChiSqPlot( TString outdir ){
  
  _chiSqRes.makeChiSqPlot(outdir);

}



int ModelInspCPAsymmetry::getNumBins(int component){
  if (component != -1) return _histsAll.at(component)->getNBins();
  
  int nbins = 0;
  for (int i = 0; i < getNumAmps(); i++){
    nbins += _histsAll.at(i)->getNBins();
  }
  return nbins;

}


void ModelInspCPAsymmetry::saveResults(TString outdir){
  _chiSqRes.saveResults(outdir);
}

void ModelInspCPAsymmetry::printChi2Breakdown(){
  _chiSqRes.print();
}

ModelInspCPAsymmetry::~ModelInspCPAsymmetry(){


}





























CPAsyChi2Res::CPAsyChi2Res(int nAmps) :
  _totChisq(0.0),
  _indChisq(nAmps, 0.0)
{
  
}

void CPAsyChi2Res::setChiSq(int amp, double chi2){
  _indChisq.at(amp) = chi2;
}

void CPAsyChi2Res::setTotChiSq(double chi2){

  if (chi2 != -1){
    _totChisq = chi2;
    return;
  }
  _totChisq = 0.0;
  for (auto iter = _indChisq.begin(); iter != _indChisq.end(); ++iter){
    _totChisq += *iter;
  }

}

double CPAsyChi2Res::getChiSq(int amp){
  return _indChisq.at(amp);
}

double CPAsyChi2Res::getTotChiSq(){
  return _totChisq;
}





int CPAsyChi2Res::getNumAmps(){
  return _indChisq.size();
}
  
CPAsyChi2Res::~CPAsyChi2Res(){
  
}












CPAsyChi2ResSet::CPAsyChi2ResSet(int nAmps) :
  _dataResults(nAmps),
  _toyStatsInd(nAmps, WidthFinder() )
{

}

CPAsyChi2ResSet::CPAsyChi2ResSet(TString filename) :
  _dataResults(0),
  _toyStatsInd(0, WidthFinder() )
{
  
  loadResults(filename);

}

int CPAsyChi2ResSet::getNumAmps(){
  return _dataResults.getNumAmps();
}

void CPAsyChi2ResSet::setDataChiSq(CPAsyChi2Res results){
  _dataResults = results;
}

void CPAsyChi2ResSet::addToyChiSq(CPAsyChi2Res results){
  _toyResults.push_back(results);

  for (int i = 0; i < getNumAmps(); i++){
    _toyStatsInd.at(i).add( results.getChiSq(i) );
  }

  _toyStatsTot.add( results.getTotChiSq() );

}

double CPAsyChi2ResSet::getDataChiSq(int component){
  if (component == -1) return _dataResults.getTotChiSq();
  return _dataResults.getChiSq(component);
}

double CPAsyChi2ResSet::getToyChiSq(int component, int toy){
  if (component == -1) return _toyResults.at(toy).getTotChiSq();
  return _toyResults.at(toy).getChiSq(component);
}

WidthFinder& CPAsyChi2ResSet::getToyStats(int component){
  if (component == -1) return _toyStatsTot;
  return _toyStatsInd.at(component);
}

int CPAsyChi2ResSet::getNumToys(){
  return _toyResults.size();
}
  
double CPAsyChi2ResSet::getNDOF(int component){
  double mean = getMean(component);
  double var  = getVariance(component);
  double scale = 0.5*(var/mean);
  double ndof  = mean/scale;  
  return ndof;
}

double CPAsyChi2ResSet::getScale(int component){
  double mean = getMean(component);
  double var  = getVariance(component);
  double scale = 0.5*(var/mean);
  return scale;
}

double CPAsyChi2ResSet::getVariance(int component){
  return getToyStats(component).varience();
}

double CPAsyChi2ResSet::getMean(int component){
  return getToyStats(component).mean();
}

double CPAsyChi2ResSet::getProb(int component){
  double ndof  = getNDOF      (component);
  double scale = getScale     (component);
  double chi2  = getDataChiSq (component);
  double prob  = TMath::Prob(chi2/scale, floor(ndof) );
  return prob;
}

double CPAsyChi2ResSet::getSig(int component){
  double prob = getProb(component);
  return fabs(TMath::NormQuantile( prob/2.0 ));
}


void CPAsyChi2ResSet::makeChiSqPlot(int component, TString outdir, bool incMeas){
  
  WidthFinder stats( getToyStats(component) );
  if (incMeas) stats.add( getDataChiSq( component ) );

  int nbins = 100;
  double min = stats.getMin() - stats.range()*0.05;
  double max = stats.getMax() + stats.range()*0.05;
  double binWid = (max - min)/double(nbins);
  
  int ntoys = getNumToys();

  TH1D hist( "chisq", "chisq", nbins, min, max );
  for (int i = 0; i < ntoys; i++){
    hist.Fill( getToyChiSq(component, i) );
  }
  
  double var  = getVariance(component);
  double mean = getMean    (component);

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

  if (incMeas) plotter.addVerticalLine( getDataChiSq(component), 1, kRed );
  plotter.plot(outdir);
  plotter.setMin(0.5);
  plotter.logY();
  plotter.plot(outdir + "_log");


}

void CPAsyChi2ResSet::makeChiSqPlot( TString outdir ){
  
  for (int i = 0; i < getNumAmps(); i++){
    TString istr = ""; istr += i;
    makeChiSqPlot( i , outdir + "_amp" + istr + "_incMeas", true  );
    makeChiSqPlot( i , outdir + "_amp" + istr             , false );
  }

  makeChiSqPlot( -1 , outdir + "_incMeas", true   );
  makeChiSqPlot( -1 , outdir             , false  );

}

void CPAsyChi2ResSet::saveResults(TString outdir){

  TFile file(outdir + ".root", "RECREATE");

  TTree* treeToys = new TTree("ChiSqFromToys", "ChiSqFromToys");
  TTree* treeData = new TTree("ChiSqFromData", "ChiSqFromData");
  
  std::vector<double*> chi2val;

  for (int i = 0; i < getNumAmps(); i++){
    TString branchname = "amp";
    branchname += i;
    chi2val.push_back(new double(0.0));
    treeToys->Branch(branchname, chi2val.back());
    treeData->Branch(branchname, chi2val.back());
  }
  
    
  for (unsigned j = 0; j < getNumToys(); j++){
    for (int i = 0; i < getNumAmps(); i++){
      *chi2val.at(i) = getToyChiSq(i, j);
    }
    treeToys->Fill();
  }

  for (int i = 0; i < getNumAmps(); i++){
    *chi2val.at(i) = getDataChiSq(i);
  }  
  treeData->Fill();

  treeToys->Write();
  treeData->Write();  

  file.Close();

}


void CPAsyChi2ResSet::loadResults(TString outdir){

  TFile file(outdir + ".root", "READ");

  TTree* treeToys = (TTree*)file.Get("ChiSqFromToys");
  TTree* treeData = (TTree*)file.Get("ChiSqFromData");
  
  int nAmps = 0;
  bool branchExists = true;
  while (branchExists){
    TString branchname = "amp";
    branchname += nAmps;
    if (treeToys->GetListOfBranches()->FindObject(branchname) == 0) break;
    nAmps++;
  }
  
  std::cout << "Loading file with " << nAmps << " components" << std::endl;

  _dataResults = CPAsyChi2Res(nAmps);
  _toyStatsInd = std::vector<WidthFinder>(nAmps, WidthFinder());

  std::vector<double*> chi2val;
  
  for (int i = 0; i < getNumAmps(); i++){
    TString branchname = "amp";
    branchname += i;
    chi2val.push_back(new double(0.0));
    treeToys->SetBranchAddress(branchname, chi2val.back());
    treeData->SetBranchAddress(branchname, chi2val.back());
  }
  

  for (int j = 0; j < treeToys->GetEntries(); j++){
    treeToys->GetEntry(j);
    
    CPAsyChi2Res toyRes( getNumAmps() );
    for (int i = 0; i < getNumAmps(); i++){
      toyRes.setChiSq(i, *chi2val.at(i));
    }
    toyRes.setTotChiSq();
    
    addToyChiSq(toyRes);

  }

  treeData->GetEntry(0);
  
  CPAsyChi2Res dataRes( getNumAmps() );
  for (int i = 0; i < getNumAmps(); i++){
    dataRes.setChiSq(i, *chi2val.at(i));
  }
  dataRes.setTotChiSq();
  
  setDataChiSq(dataRes); 


  file.Close();

}



void CPAsyChi2ResSet::print(){

  
  std::cout << "--------------- ALL HISTS ----------------" << std::endl;
  std::cout << "chi2  = " << getDataChiSq(-1) << std::endl; 
  std::cout << "ndof  = " << getNDOF(-1) << std::endl;
  std::cout << "scale = " << getScale(-1) << std::endl;
  std::cout << "pval  = " << getProb(-1) << std::endl;
  std::cout << "sigma = " << getSig (-1) << std::endl;

  for (int i = 0; i < getNumAmps(); i++){
    std::cout << "------------- HIST " << i << " ----------------" << std::endl;
    std::cout << "chi2  = " << getDataChiSq(i) << std::endl; 
    std::cout << "ndof  = " << getNDOF (i) << std::endl;
    std::cout << "scale = " << getScale(i) << std::endl;    
    std::cout << "pval  = " << getProb(i) << std::endl;
    std::cout << "sigma = " << getSig (i) << std::endl;
  }



}





CPAsyChi2ResSet::~CPAsyChi2ResSet(){

}





