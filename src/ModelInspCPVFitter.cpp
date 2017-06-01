#include "ModelInspCPVFitter.h"







ModelInspCPVFitter::ModelInspCPVFitter(DalitzEventListWithAmps& dzEvts, DalitzEventListWithAmps& dzbEvts) :
  _parSet( new MINT::MinuitParameterSet() ),
  _minimiser(0),
  _dzEvts (dzEvts ),
  _dzbEvts(dzbEvts),
  _globAys( "GlobalAys" , 0, 0.0 , 0.01  ,  -1.0 , 1.0, _parSet),
  _currentComponent(0)
{

  setPset(_parSet);


  int nComps = _dzEvts.at(0).getNumAmps();
  
  for (int i = 0; i < nComps; i++){
  	TString istr = ""; istr += i;

    _realPerts.push_back( new MINT::FitParameter(("AmpPert_" + istr + "_Re").Data(), 2, 0.0, 0.01, 0.0, 0.0, _parSet) );
    _imagPerts.push_back( new MINT::FitParameter(("AmpPert_" + istr + "_Im").Data(), 2, 0.0, 0.01, 0.0, 0.0, _parSet) );
  }



}
  

MINT::Minimiser* ModelInspCPVFitter::getMinimiser(){
  
  if (_minimiser == 0){
    _minimiser = new MINT::Minimiser(this);
  }
  return _minimiser;

}

double ModelInspCPVFitter::getProbDzNoNorm( int evt, int ampComp ){
  
  std::complex<double> dzbpeturb(  
    *_realPerts.at(ampComp),
    *_imagPerts.at(ampComp)
  ); 

  double globalAys = _globAys;

  double probDzNoNorm  = (1.0 - globalAys)*_dzEvts.at(evt).getProbWithPerturb(ampComp, -dzbpeturb);
  double probDzbNoNorm = (1.0 + globalAys)*_dzEvts.at(evt).getProbWithPerturb(ampComp,  dzbpeturb);

  return probDzNoNorm/(probDzNoNorm + probDzbNoNorm);

}

double ModelInspCPVFitter::getProbDzbNoNorm( int evt, int ampComp ){
  
  std::complex<double> dzbpeturb(  
    *_realPerts.at(ampComp),
    *_imagPerts.at(ampComp)
  ); 

  double globalAys = _globAys;

  double probDzNoNorm  = (1.0 - globalAys)*_dzbEvts.at(evt).getProbWithPerturb(ampComp, -dzbpeturb);
  double probDzbNoNorm = (1.0 + globalAys)*_dzbEvts.at(evt).getProbWithPerturb(ampComp,  dzbpeturb);

  return probDzbNoNorm/(probDzNoNorm + probDzbNoNorm);

}

double ModelInspCPVFitter::getNeg2LLHDz(int comp){
  
  double neg2LLH = 0.0;
  for (int i = 0; i < _dzEvts.size(); i++){
    neg2LLH += 2.0*log ( getProbDzNoNorm( i, comp ) );
  }
  return -neg2LLH;

}

double ModelInspCPVFitter::getNeg2LLHDzb(int comp){
  
  double neg2LLH = 0.0;
  for (int i = 0; i < _dzbEvts.size(); i++){
    neg2LLH += 2.0*log ( getProbDzbNoNorm( i, comp ) );
  }
  return -neg2LLH;
  
}


void ModelInspCPVFitter::makePlot(TString outdir){
    
  HyperPointSet pointsDz  (2);
  HyperPointSet pointsDzb (2);
  HyperPointSet pointsBoth(2);

  for (int i = 0; i < _dzEvts.size(); i++){
  	double phase = arg(_dzEvts.at(i).getAmpRelToTot(_currentComponent));
  	double ratio = abs(_dzEvts.at(i).getAmpRelToTot(_currentComponent));
  	HyperPoint point(phase, ratio);
  	pointsDz.push_back(point);
  	pointsBoth.push_back(point);
  }
  for (int i = 0; i < _dzbEvts.size(); i++){
  	double phase = arg(_dzbEvts.at(i).getAmpRelToTot(_currentComponent));
  	double ratio = abs(_dzbEvts.at(i).getAmpRelToTot(_currentComponent));
  	HyperPoint point(phase, ratio);
  	pointsDzb.push_back(point);
  	pointsBoth.push_back(point);
  }
  
  HyperPoint  minWidth( TMath::Pi()*0.05, 0.01 ); 
  HyperCuboid limits( HyperPoint(-TMath::Pi(), 0.0), HyperPoint(+TMath::Pi(), 1.0) );

  HyperHistogram hyperHist(limits, pointsBoth, HyperBinningAlgorithms::SMART , 
  	                                           AlgOption::MinBinWidth(minWidth), 
  	                                           AlgOption::MinBinContent(100) );

  HyperHistogram hyperHistDz (hyperHist.getBinning());
  HyperHistogram hyperHistDzb(hyperHist.getBinning());
  HyperHistogram hyperHistAys(hyperHist.getBinning());

  hyperHistDz .fill(pointsDz );
  hyperHistDzb.fill(pointsDzb);

  HyperPlotStyle::setPalette("birdy");

  hyperHistDz .drawDensity(outdir + "_DzHyp.pdf" , "COLZ Edges1");
  hyperHistDzb.drawDensity(outdir + "_DzbHyp.pdf", "COLZ Edges1");
  
  double significance = hyperHistDz.chi2sig(hyperHistDzb);
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "significance = " << significance << std::endl;
  std::cout << "-----------------------------------------" << std::endl;

  HyperPlotStyle::setPalette("pulls");

  hyperHistAys.asymmetry   (hyperHistDzb, hyperHistDz);
  hyperHistAys.setMin(-1.0);
  hyperHistAys.setMax( 1.0);
  hyperHistAys.draw(outdir + "_divHyp.pdf", "COLZ Edges1");


  //std::complex<double> dzbpeturb(  
  //  *_realPerts.at(_currentComponent),
  //  *_imagPerts.at(_currentComponent)
  //); 
  //
  //for (int i = 1; i <= nbins; i++){
  //  for (int j = 1; j <= nbins; j++){
  //  	double nDz  = histDz .GetBinContent(i,j);
  //  	double nDzb = histDzb.GetBinContent(i,j);
  //    if (nDz == 0.0 && nDzb == 0.0) histDiv.SetBinContent(i,j, 0.0);
  //    else histDiv.SetBinContent( i,j, (nDzb - nDz)/(nDz) );
  //    
  //    double phi = histExpDiv.GetXaxis()->GetBinCenter(i);
  //    double r   = histExpDiv.GetYaxis()->GetBinCenter(j);
  //    double expect = 2.0*r*( cos(phi)*dzbpeturb.real() - sin(phi)*dzbpeturb.imag() ) / (1.0 + r*r);
  //    histExpDiv.SetBinContent(i,j, expect );
  //  }  	
  //}
  


}

double ModelInspCPVFitter::getVal(){
  
  return getNeg2LLHDz(_currentComponent) + getNeg2LLHDzb(_currentComponent);

}
 
void ModelInspCPVFitter::setComponentToFit(int i){
  _currentComponent = i;
}

void ModelInspCPVFitter::fixComponents(int compNum, bool unfix){

  for (unsigned i = 0; i < _realPerts.size(); i++){
  	if (compNum == -1 || i == compNum){
      if (unfix == true){
      	_realPerts.at(i)->unFix();
      	_imagPerts.at(i)->unFix();
      }
      else{
      	_realPerts.at(i)->fix();
      	_imagPerts.at(i)->fix();    	
      }
    }
  }

}

void ModelInspCPVFitter::resetComponents(){

  for (unsigned i = 0; i < _realPerts.size(); i++){

    _realPerts.at(i)->setCurrentFitVal(0.0);
    _imagPerts.at(i)->setCurrentFitVal(0.0);

  }  
  
  _globAys.setCurrentFitVal(0.0);
}



int ModelInspCPVFitter::fit(bool minos ){
  
  resetComponents();  
  fixComponents(-1               , false);
  fixComponents(_currentComponent, true);

  MINT::Minimiser* mini = 0;
  if (mini == 0) mini = getMinimiser();

  bool migradHesseSuccess = 1;

  if (minos == 1) migradHesseSuccess = mini->continueMinosFit();
  else            migradHesseSuccess = mini->continueFit();  

  int errCode = 0;
  TString errMsg = "FIT SUCCESSFUL";

  double neg2llh = 0.0;
  double edm     = 0.0;
  double errdef  = 0.0;
  int nparfloat  = 0;
  int npar       = 0;
  int istat      = 0;  
  //0 - not calculated
  //1 - only approximate
  //2 - not pos def
  //3 - full accurate

  mini->mnstat(neg2llh, edm, errdef, nparfloat, npar, istat);

  if      ( migradHesseSuccess == false   ) { errCode = 1; errMsg = "MIGRAD / HESSE FAILED";         }
  else if ( edm                > 0.01     ) { errCode = 4; errMsg = "BAD EDM";                       }
  else if ( neg2llh            != neg2llh ) { errCode = 5; errMsg = "NEG2LLH IS NaN";                }
  else if ( istat              == 0       ) { errCode = 6; errMsg = "COVARIENCE MATRIX NOT CALC";    }
  else if ( istat              == 1       ) { errCode = 7; errMsg = "COVARIENCE MATRIX ONLY APPROX"; }
  else if ( istat              == 2       ) { errCode = 8; errMsg = "COVARIENCE MATRIX NOT POS DEF"; }
  

  std::cout << "********************** ObservableFitter::fit() **********************" << std::endl;
  std::cout << "*********************************************************************" << std::endl;
  std::cout << "                         " << errMsg << std::endl;
  std::cout << "*********************************************************************" << std::endl;

  return errCode;

}


ModelInspCPVFitter::~ModelInspCPVFitter(){



}

