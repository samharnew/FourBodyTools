#include "PhilippeModelWrapper.h"


PhilippeModelWrapper* PhilippeModelWrapper::s_wrapper = 0;

PhilippeModelWrapper::PhilippeModelWrapper() :
  pat    (421, -211, 211, 211, -211),
  fas    ( new FitAmpSum(pat) ),
  fascc  (0),
  a1_re  ("a1_Re"  ),
  a1_im  ("a1_Im"  ),
  pi_re  ("pi_Re"  ),
  pi_im  ("pi_Im"  ),
  pi2_re ("pi2_Re" ),
  pi2_im ("pi2_Im" ),
  a1p_re ("a1p_Re" ),
  a1p_im ("a1p_Im" ),
  fas_a1 (0        ),
  fas_pi (0        ),
  fas_pi2(0        ),
  fas_a1p(0        ),
  ModelName("ModelName", (std::string) "", (char*) 0)

{
  
  //Just need to pick some seed to make the normalised amplitude 
  //parameters reproducable 

  TRandom3 random(4);

  int nPhasespace = 200000; 

  eventListMC.generatePhaseSpaceEvents(nPhasespace,pat, &random);  

  fas->getVal(eventListMC[0]);   
  
  if (a1_re.getCurrentFitVal() != 0.0 || a1_im.getCurrentFitVal() != 0.0){
    a1List = fas->GetCloneOfSubsetSameFitParameters("a(1)(1260)+");
    fas_a1 = new FitAmpSum(*a1List);
    fas_a1->getVal(eventListMC[0]);
    fas_a1->CConjugateFinalStateSameFitParameters();
    ar_a1 = new AmpRatio(a1_re,a1_im);
    fas_a1->multiply(ar_a1);
    fas->addAsList(*fas_a1,1.);
  }
  
  if (pi_re.getCurrentFitVal() != 0.0 || pi_im.getCurrentFitVal() != 0.0){
    piList = fas->GetCloneOfSubsetSameFitParameters("pi(1300)");
    fas_pi = new FitAmpSum(*piList);
    fas_pi->getVal(eventListMC[0]);
    fas_pi->CConjugateFinalStateSameFitParameters();
    ar_pi = new AmpRatio(pi_re,pi_im);
    fas_pi->multiply(ar_pi);
    fas->addAsList(*fas_pi,1.);
  }

  if (pi2_re.getCurrentFitVal() != 0.0 || pi2_im.getCurrentFitVal() != 0.0){
    pi2List = fas->GetCloneOfSubsetSameFitParameters("pi(2)(1670)");
    fas_pi2 = new FitAmpSum(*pi2List);
    fas_pi2->getVal(eventListMC[0]);
    fas_pi2->CConjugateFinalStateSameFitParameters();
    ar_pi2 = new AmpRatio(pi2_re,pi2_im);
    fas_pi2->multiply(ar_pi2);
    fas->addAsList(*fas_pi2,1.);
  }

  if (a1p_re.getCurrentFitVal() != 0.0 || a1p_im.getCurrentFitVal() != 0.0){
    a1pList = fas->GetCloneOfSubsetSameFitParameters("a(1)(1640)");
    fas_a1p = new FitAmpSum(*a1pList);
    fas_a1p->getVal(eventListMC[0]);
    fas_a1p->CConjugateFinalStateSameFitParameters();
    ar_a1p = new AmpRatio(a1p_re,a1p_im);
    fas_a1p->multiply(ar_a1p);
    fas->addAsList(*fas_a1p,1.);
  }

  fas->getVal(eventListMC[0]);
  fas->print();
  fas->normalizeAmps(eventListMC);
  

  copyList = fas->GetCloneSameFitParameters();
  fascc = new FitAmpSum(*copyList);
  fascc->getVal(eventListMC[0]);
  fascc->CPConjugateSameFitParameters();
  fascc->getVal(eventListMC[0]);

}

TString PhilippeModelWrapper::getModelName(){

  return (TString)ModelName.getVal();

}



MINT::counted_ptr<FitAmpSum> PhilippeModelWrapper::getFitAmpSum(){

  return fas;

}

MINT::counted_ptr<FitAmpSum> PhilippeModelWrapper::getFitAmpSumCC(){

  return fascc;


}


PhilippeModelWrapper* PhilippeModelWrapper::getStaticWrapper(){
  if ( s_wrapper == 0 ) s_wrapper = new PhilippeModelWrapper();
  return s_wrapper;
}


void PhilippeModelWrapper::radomiseModel(int seed){

  TRandom3 random(seed);

  for(unsigned int i=0; i < fas->size(); i++){

    MINT::FitParameter& realPar = fas->getAmpPtr(i)->FitAmpPhase().p1();
    MINT::FitParameter& imagPar = fas->getAmpPtr(i)->FitAmpPhase().p2();
    
    realPar.setCurrentValToInit();
    imagPar.setCurrentValToInit();

    double realVal = realPar.mean();
    double realErr = realPar.err ();

    double imagVal = imagPar.mean();
    double imagErr = imagPar.err ();

    double randRealVal = random.Gaus(realVal, realErr);
    double randImagVal = random.Gaus(imagVal, imagErr);

    realPar.setCurrentFitVal(randRealVal);
    imagPar.setCurrentFitVal(randImagVal);

  }  
  
  {
    a1_re .setCurrentValToInit();
    double val    = a1_re.mean();
    double err    = a1_re.err ();
    double ranVal = random.Gaus(val, err);
    a1_re.setCurrentFitVal(ranVal);
  }

  {
    a1_im .setCurrentValToInit();
    double val    = a1_im.mean();
    double err    = a1_im.err ();
    double ranVal = random.Gaus(val, err);
    a1_im.setCurrentFitVal(ranVal);
  }

  {
    pi_re .setCurrentValToInit();
    double val    = pi_re.mean();
    double err    = pi_re.err ();
    double ranVal = random.Gaus(val, err);
    pi_re.setCurrentFitVal(ranVal);
  }

  {
    pi_im .setCurrentValToInit();
    double val    = pi_im.mean();
    double err    = pi_im.err ();
    double ranVal = random.Gaus(val, err);
    pi_im.setCurrentFitVal(ranVal);
  }

  {
    pi2_re .setCurrentValToInit();
    double val    = pi2_re.mean();
    double err    = pi2_re.err ();
    double ranVal = random.Gaus(val, err);
    pi2_re.setCurrentFitVal(ranVal);
  }

  {
    pi2_im .setCurrentValToInit();
    double val    = pi2_im.mean();
    double err    = pi2_im.err ();
    double ranVal = random.Gaus(val, err);
    pi2_im.setCurrentFitVal(ranVal);
  }

  {
    a1p_re .setCurrentValToInit();
    double val    = a1p_re.mean();
    double err    = a1p_re.err ();
    double ranVal = random.Gaus(val, err);
    a1p_re.setCurrentFitVal(ranVal);
  }

  {
    a1p_im .setCurrentValToInit();
    double val    = a1p_im.mean();
    double err    = a1p_im.err ();
    double ranVal = random.Gaus(val, err);
    a1p_im.setCurrentFitVal(ranVal);
  }            

  //MINT::MinuitParameterSet::getDefaultSet()->print();

}

void PhilippeModelWrapper::resetModel(){

  for(unsigned int i=0; i < fas->size(); i++){

    MINT::FitParameter& realPar = fas->getAmpPtr(i)->FitAmpPhase().p1();
    MINT::FitParameter& imagPar = fas->getAmpPtr(i)->FitAmpPhase().p2();
    
    realPar.setCurrentValToInit();
    imagPar.setCurrentValToInit();

  }  
  
  a1_re .setCurrentValToInit();
  a1_im .setCurrentValToInit();
  pi_re .setCurrentValToInit();
  pi_im .setCurrentValToInit();
  pi2_re.setCurrentValToInit();
  pi2_im.setCurrentValToInit();
  a1p_re.setCurrentValToInit();
  a1p_im.setCurrentValToInit();

  //MINT::MinuitParameterSet::getDefaultSet()->print();

}


PhilippeModelWrapper::~PhilippeModelWrapper(){
  
  fascc.release();
  copyList.release();
  

  fas    .release();

  ar_a1 .release();
  ar_pi .release();
  ar_pi2.release();
  ar_a1p.release();

  delete fas_a1 ;
  delete fas_pi ;
  delete fas_pi2;
  delete fas_a1p;

  a1List .release();
  piList .release();
  pi2List.release();
  a1pList.release();


}




