#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

//' C version of the VSEM model
//' @param par parameter vector
//' @param PAR 	Photosynthetically active radiation (PAR) MJ /m2 /day
//' @export
// [[Rcpp::export]]
NumericMatrix vsemC(NumericVector par, NumericVector PAR){
    
  int numObs = PAR.size();

  // Parameter definitions
  
  double KEXT  = par[0];
  double LAR   = par[1];
  double LUE   = par[2];
  double GAMMA = par[3];
  double tauV  = par[4];
  double tauS  = par[5];
  double tauR  = par[6]; // 1440;
  double Av    = par[7]; //0.5;
  double Cv    = par[8]; //3.0;
  double Cs    = par[9]; //15;
  double Cr    = par[10]; //3.0r;
  double G;
  double NEE;
  double NPP;
  
  NumericMatrix out(numObs, 4);
  
  for (int i = 0; i < numObs; ++i){
    G   = PAR[i] * LUE * (1 - exp(-KEXT*LAR*Cv)) ;
    NPP = (1-GAMMA)*G;
    Cv  = Cv + Av*NPP - Cv/tauV;
    Cr  = Cr + (1.0-Av)*NPP - Cr/tauR;
    Cs  = Cs + Cr/tauR + Cv/tauV - Cs/tauS;
    NEE = (Cs/tauS + GAMMA*G) - G;    
    out(i,0) = NEE;
    out(i,1) = Cv;
    out(i,2) = Cs;
    out(i,3) = Cr;
  }
  
  return out;
}