// Toy model for testing Sommerfeld Corrections

//#define DEBUG
#include "SommCorr.hh"

Sommerfeld::~Sommerfeld(){
}


void Sommerfeld::init(double Min, double dMin, double vin, bool adapt){

  // M = 2749.4;
  // dM = 0.21;
  // v = 0.012;
  
  M = Min;
  dM = dMin;
  v = vin/2.0; // Beneke notation vrel = 2*v
  doAdaptive = adapt;
  
  alp2 = 0.0315167;
  alpEM = 7.30e-3;
  mZ = 91.1987;
  mW = 80.423;
  mH = 125.5;

  rW = mW / (M * v);
  rZ = mZ / (M * v);
  xbar = alp2 / v * sqrt(2);
  cW2 = 0.769;

  GG[0][1] = 0.0;
  GG[1][0] = 0.0;
  kk1 = 1.0;
  double kk2temp =  2 * dM / (M * v * v);
  if(kk2temp > 1){
    kk2 = sqrt(kk2temp-1);
    GG[0][0] = cplx(0, 1.0) * kk1;
    GG[1][1] = cplx(-1.0, 0.0) * kk2;

  }else{
    kk2 = sqrt(1-kk2temp);
    GG[0][0] = cplx(0,1.0) * kk1;
    GG[1][1] = cplx(0.0,1.0) * kk2;

  }

  eps = 1.0e-7;
  maxErr = 0.0;
  Ntot = 0;
  Nsize = 2;

#ifdef DEBUG
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++){
      cout << "GG[" << i << "][" << j << "]="<< GG[i][j] <<endl;
    }
  cout << "dM = "  << dM;
  cout << " kk2^2 = " << 1 - 2 * dM / (M*v*v) << endl;
  cout << "rW = " << rW << " rZ = " << rZ << " xbar = " << xbar << endl;
  cout << "coeff of Z-term =" << xbar * cW2 / sqrt(2) << "   coeff of gamma term = " << alpEM / v <<endl; 
#endif  

}

void Sommerfeld::zero(){

  // Set all arrays to initial values
  
  for (int i = 0; i < Nsize; i++){
    for (int j = 0; j < Nsize; j++){
      N[i][j] = (0.0,0.0); // initial value
      Nin[i][j] = (0.0,0.0);
      A[i][j] = 0.0;
      Ain[i][j] = 0.0;
      for (int ik = 0; ik < 5; ik++){
	k[ik][i][j] = (0.0,0.0);
	ka[ik][i][j] = (0.0,0.0);
      }
      if(i == j) II[i][j] = 1.0;
      else II[i][j] = 0.0;
    }
  }

  N[0][0] = cplx(eps,0);
  N[1][1] = cplx(eps,0);

  A[0][0] = 1.0;
  A[1][1] = 1.0;
  
  return;
}



void Sommerfeld::getGammaP(){
  
  cplx Gamma[4][4], tGam[4][4];

  Gamma[0][0] = 1.0/2.0;
  Gamma[0][1] = 1.0/2.0;
  Gamma[1][0] = Gamma[0][1];
  Gamma[1][1] = 3.0/2.0;

  // Gamma for method 1
  // Gamma[0][0] = 1.0;
  // Gamma[0][1] = 1.0/sqrt(2.0);
  // Gamma[1][0] = Gamma[0][1];
  // Gamma[1][1] = 3.0/2.0;

  for (int i = 0; i < Nsize; i++){
    for (int j = 0; j < Nsize; j++) {
      tGam[i][j] = 0.0;
      GammaP[i][j] = 0.0;
      for (int l = 0; l < Nsize; l++) 
	tGam[i][j] += Gamma[i][l] * A[l][j];  
    }
  }

  for (int i = 0; i < Nsize; i++)
    for (int j = 0; j < Nsize; j++) 
      for (int l = 0; l < Nsize; l++) 
	GammaP[i][j] += conj(A[l][i]) * tGam[l][j];  
}
  
void Sommerfeld::printOutput(bool shortOP){
  
  if(shortOP) {
    cout << real(GammaP[0][0]) << endl;
    return;
  }
  for (int i = 0; i < Nsize; i++){
    for (int j = 0; j < Nsize; j++) {
      cout << "N[" << i << "][" << j << "]="<< N[i][j] << "  ";
    }
    cout << endl;
  }
  for (int i = 0; i < Nsize; i++){
    for (int j = 0; j < Nsize; j++) {
      cout << "A[" << i << "][" << j << "]="<< A[i][j] << "  ";
    }
    cout << endl;
  }

  cout << "Sommerfeld correction matrix = " << endl;
  for (int i = 0; i < Nsize; i++){
    for (int j = 0; j < Nsize; j++) {
      cout << "GammaP[" << i << "][" << j << "]="<< GammaP[i][j] << "  ";
    }
    cout << endl;
  }

  cout << "Correction factor = "<< real(GammaP[0][0]) <<endl;
  cout << "Maximum err: " << maxErr  << endl;
  cout << "Performed " << Ntot << " iterations "<< endl;
  
}
void Sommerfeld::eulerN(double xmax){
  
  double h = 0.0001;
  double x = 1.0e-7;
  zero();
  unsigned int n = 0;

  
  while (x < xmax) {

#ifdef DEBUG2
    cout << "Starting x =" << x << endl;
#endif
    
    for (int i = 0; i < Nsize; i++)
      for (int j = 0; j < Nsize; j++) {
	Nin[i][j] = N[i][j];
	Ain[i][j] = A[i][j];
      }
    setV(x);
    calculate(1);
    for (int i = 0; i < Nsize; i++)
      for (int j = 0; j < Nsize; j++) {
	N[i][j] += h * k[1][i][j];
	A[i][j] += h * ka[1][i][j];
      }

    x += h;
    n++;
  }

  getGammaP();
#ifdef DEBUG  
  cout << "Performed " << n << " iterations " << xmax << endl;
#endif  

}

void Sommerfeld::midpointN(double xmax){
  
  double h = 1.0e-4;
  double x = 1.0e-7;
  zero();

#ifdef DEBUG
  ofstream Nval;
  Nval.open("nval.txt");
#endif    

  unsigned int n = 0;
  while (x < xmax) {

#ifdef DEBUG
    Nval << x << "  " << abs(N[0][0]) << "  " << abs(N[0][1])
	 << "  " << abs(N[1][1]) << endl;
#endif

    
    for (int i = 0; i < Nsize; i++)
      for (int j = 0; j < Nsize; j++) {
	Nold[i][j] = N[i][j];
	Nin[i][j] = N[i][j];
	Ain[i][j] = A[i][j];
      }
    
    setV(x);
    calculate(1);
    
    for (int i = 0; i < Nsize; i++)
      for (int j = 0; j < Nsize; j++) {
	Nin[i][j] = N[i][j] + h * k[1][i][j] / 2.0;
	Ain[i][j] = A[i][j] + h * ka[1][i][j] / 2.0;
      }

    setV(x+h/Nsize);
    calculate(Nsize);
    
    for (int i = 0; i < Nsize; i++)
      for (int j = 0; j < Nsize; j++) {
	N[i][j] += h * k[2][i][j];
	A[i][j] += h * ka[2][i][j];
      }

    
    x += h;
    n++;
  }
  getGammaP();
  
#ifdef DEBUG  

#endif  
}

void Sommerfeld::rk4N(double xmax){
  
  double h = 1.0e-4;
  double x = 1.0e-7;
  zero();

  while (x < xmax) {
    bool accept = false;
    int n = 0;
    Ntot++;

    if (Ntot > 100000000) {
      cout << "Too many iterations" <<endl;
      return;
    }
    // Loop till step size enough for tolerance or more than 10 h-iterations for single x
    while (!accept && n < 10){
      n++;
	
      for (int i = 0; i < Nsize; i++)
	for (int j = 0; j < Nsize; j++) {
	  Nold[i][j] = N[i][j];
	  Nin[i][j] = N[i][j];
	  Ain[i][j] = A[i][j];
	}
    
      setV(x);
      calculate(1);

      for (int i = 0; i < Nsize; i++)
	for (int j = 0; j < Nsize; j++) {
	  Nin[i][j] = N[i][j] + h * k[1][i][j] / 2.0;
	  Ain[i][j] = A[i][j] + h * ka[1][i][j] / 2.0;
	}

      setV(x+h/2);
      calculate(2);

      for (int i = 0; i < Nsize; i++)
	for (int j = 0; j < Nsize; j++) {
	  Nin[i][j] = N[i][j] + h * k[2][i][j] / 2.0;
	  Ain[i][j] = A[i][j] + h * ka[2][i][j] / 2.0;
	}
      
      calculate(3);    

      for (int i = 0; i < Nsize; i++)
	for (int j = 0; j < Nsize; j++) {
	  Nin[i][j] = N[i][j] + h * k[3][i][j];
	  Ain[i][j] = A[i][j] + h * ka[3][i][j];
	}
    
      setV(x+h);
      calculate(4);
      for (int i = 0; i < Nsize; i++)
	for (int j = 0; j < Nsize; j++) {
	  N[i][j] += h * (k[1][i][j] + k[4][i][j]) / 6.0;
	  N[i][j] += h * (k[2][i][j] + k[3][i][j]) / 3.0; 
	  A[i][j] += h * (ka[1][i][j] + ka[4][i][j]) / 6.0;
	  A[i][j] += h * (ka[2][i][j] + ka[3][i][j]) / 3.0; 
	}

      // Check accuracy for adaptive stepsize
      // local: scale = rtol + |y| atol ~= tol * (1 + |y|)
      // global: scale = tol * h * dy/dx; exponent -> 1/4
      
      double delta = 10.0, scale = 10.0;
      bool local = true;
      
      for (int i = 0; i < Nsize; i++)
	for (int j = 0; j < Nsize; j++) {
	  delta = min(delta, abs(N[i][j] - Nold[i][j]));
	  if (local)
	    scale = 1 + max(max(abs(N[i][j]), abs(Nold[i][j])), scale);
	  else {
	    double mink = 1000.0;
	    for (int i = 0; i < 2; i++)
	      for (int j = 0; j < 2; j++) 
		if (abs(k[1][i][j]) < mink) mink = abs(k[1][i][j]);
	    scale = h * mink;
	  }
	}

      double err = delta / scale;
      // Is error within tolerance?
      if(err < tol || !doAdaptive) {
	maxErr = max(err, maxErr);
	accept = true;
	break;
      }
      // Change stepsize (Num.Rec.3rd; eqn 17.2.12)
      if (local)
	h = (0.95) * h / pow(err/tol , 1.0/5.0);
      else
	h = (0.95) * h / pow(err/tol , 1.0/4.0);
    }
    x += h;
  }
  getGammaP();
}


void Sommerfeld::calculate(int ik){

  cplx temp1[4][4], temp2[4][4], temp3[4][4], temp4[4][4];
  cplx temp1a[4][4], temp2a[4][4], temp3a[4][4];
  
  for (int i = 0; i < Nsize; i++){
    for (int j = 0; j < Nsize; j++){
      temp1[i][j] = cplx(0.0,0.0);
      temp2[i][j] = cplx(0.0,0.0);
      temp3[i][j] = cplx(0.0,0.0);
      temp4[i][j] = cplx(0.0,0.0);      

      temp1a[i][j] = cplx(0.0,0.0);
      temp2a[i][j] = cplx(0.0,0.0);
      temp3a[i][j] = 0.0;

      k[ik][i][j] = 0.0;
      ka[ik][i][j] = 0.0;

      for (int l = 0; l < Nsize; l++){
	temp1[i][j] += GG[i][l] * Nin[l][j]; 
	temp2[i][j] += Nin[i][l] * GG[l][j];
	temp3[i][j] += V[i][l] * Nin[l][j];
	temp2a[i][j] += V[i][l] * N[l][j];
      }
    }
  }
  for (int i = 0; i < Nsize; i++){
    for (int j = 0; j < Nsize; j++) {
      for (int l = 0; l <Nsize; l++) {
	temp4[i][j] += Nin[i][l] * temp3[l][j];
	temp3a[i][j] += A[i][l] * (GG[l][j] - temp2a[l][j]) ;
	      
      }
      k[ik][i][j] = II[i][j] + temp1[i][j] + temp2[i][j] - temp4[i][j];
      ka[ik][i][j] = temp3a[i][j];

    }
  }
  return;
}

void Sommerfeld::setV(double x){

  V[0][0] = 0.0;
  V[0][1] = -xbar / x * exp(-rW * x);
  V[1][0] = V[0][1];
  V[1][1] = - (alpEM/v + xbar * cW2 / sqrt(2) * exp(-rZ * x))/x;
  
#ifdef DEBUG2
  for (int i = 0; i < Nsize; i++)
    for (int j = 0; j < Nsize; j++)
      cout << "  V[" << i << "][" << j << "]="<< V[i][j] <<endl;
#endif

  return;
}

int main(int argc, char** argv){

  double M, dM, v, xmax;
  bool adapt = true;
  if (argc < 5) {
    cout << "Usage: ./somm.x <M> <dM> <v> <xmax> [doAdaptive = 0, 1] "<< endl;
    return -1;
  }
  else {
    M  = atof(argv[1]);
    dM = atof(argv[2]);
    v  = atof(argv[3]);
    xmax = atof(argv[4]);
    if(argc == 6) adapt = atoi(argv[5]);
  }
  

  Sommerfeld somm(M,dM,v,adapt);
  somm.tol = 1.0e-5;
  somm.rk4N(xmax);
  somm.printOutput(false);

  // Plot potential
  // ofstream vout;
  // vout.open("vout.txt");
  // for (int i = 1; i <=200; i++){
  //   double x = 0.01*i;
  //   somm.setV(x);
  //   vout << x << "  " <<somm.V[1][0] << "  " << somm.V[1][1] << endl;
  // }

  //To determine the dependence of S on xmax
  // ofstream output;
  // output.open("output.txt");
  // for (int i = 1; i <= 1000; i++){
  //   somm.rk4N(xmax/(i*1.0));
  //   output << xmax/(i*1.0) << "  " <<abs(somm.N[0][0]) << "  "
  // 	   << abs(somm.A[0][0]) << "  " <<abs(somm.GammaP[0][0]) << endl;
  // }

  // Determine dependence on deltam
  // for (int i = 1; i < 30; i++){
  //   double vin = v * i / 10.0;
  //   Sommerfeld somm(M,dM,vin,adapt);
  //   somm.tol = 1.0e-5;
  //   somm.rk4N(xmax);
  //   cout << vin << " ";
  //   somm.printOutput(true);
  // }


  return 0;
}


