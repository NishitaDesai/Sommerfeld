
#include<iostream>
#include<cmath>
#include<complex>
#include<fstream>
#include<stdlib.h>

using namespace std;

typedef complex<double> cplx;

class Sommerfeld{

public:
  Sommerfeld(double Min, double dM, double vin, bool adapt){
    init(Min, dM, vin, adapt);
  };
  ~Sommerfeld();

  void eulerN(double xmax);
  void rk4N(double xmax);
  void midpointN(double xmax);
  void printOutput(bool shortOP = false);
  void setV(double x);

  cplx A[4][4], N[4][4], GammaP[4][4];
  double V[4][4];
  double maxErr;
  bool  doAdaptive;
  double tol;
  
private:
  void getGammaP();
  void calculate(int ik);
  void init(double Min, double dM, double vin, bool adapt);
  void zero();
  
  double M, dM, v, KE, alp2, alpEM, mZ, mW, mH;
  double rW, rZ, xbar, cW2, eps;

  int Ntot;
  int Nsize;

  double II[4][4];
  cplx GG[4][4];

  cplx  Nin[4][4], Nold[4][4];
  cplx k[5][4][4];

  cplx  Ain[4][4];
  cplx ka[5][4][4];
  
  double kk1, kk2;

};
