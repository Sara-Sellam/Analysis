#include "Riostream.h" 

#include "Bkg.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

// Define some constants
Double_t mK1 = 493.667;
//Double_t mPi = 139.57018; // not needed so far.


// Momentum of a daughter in the resonance rest frame
Double_t get_q2( Double_t m, Double_t mA, Double_t mB)
{
  Double_t num = (m*m - (mA + mB)*(mA + mB)) * (m*m - (mA - mB)*(mA - mB));
  if (num<0) return 0.;

  return sqrt(num)/(2.*m);
}


ClassImp(Bkg) 

 Bkg::Bkg(const char *name, const char *title, 
			RooAbsReal& _m,
			RooAbsReal& _a1,
                        RooAbsReal& _a2) :
  RooAbsPdf(name,title), 
  m("m","m",this,_m),
  a1("a1","a1",this,_a1),
  a2("a2","a2",this,_a2)
 { 
 }

 Bkg::Bkg(const Bkg& other, const char* name) :  
   RooAbsPdf(other,name), 
   m("m",this,other.m),
   a1("a1",this,other.a1),
   a2("a2",this,other.a2)
 { 
 }

 Double_t Bkg::evaluate() const 
 { 

  Double_t q = get_q2(m,mK1,mK1);
  Double_t mb = m-2.*mK1;

  Double_t pdf = q*pow(mb,a1)*exp(-a2*mb);
  //Double_t pdf = q*(1+a1*mb+a2*mb*mb);

  if (pdf<0.) return 0.;
  return pdf;
}
