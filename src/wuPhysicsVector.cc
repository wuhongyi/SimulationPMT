// wuPhysicsVector.cc --- 
// 
// Description: Coping from Geant4 G4PhysicsVector.cc
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 四 2月  5 20:21:05 2015 (+0800)
// Last-Updated: 四 2月  5 20:59:28 2015 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 3
// URL: http://wuhongyi.cn 

// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4PhysicsVector.cc
//
//  History:
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    03 Mar. 1996, K.Amako : Implemented the 1st version
//    01 Jul. 1996, K.Amako : Hidden bin from the user introduced
//    12 Nov. 1998, K.Amako : A bug in GetVectorLength() fixed
//    11 Nov. 2000, H.Kurashige : use STL vector for dataVector and binVector
//    18 Jan. 2001, H.Kurashige : removed ptrNextTable
//    09 Mar. 2001, H.Kurashige : added G4PhysicsVector type 
//    05 Sep. 2008, V.Ivanchenko : added protections for zero-length vector
//    11 May  2009, A.Bagulya : added new implementation of methods 
//            ComputeSecondDerivatives - first derivatives at edge points 
//                                       should be provided by a user
//            FillSecondDerivatives - default computation base on "not-a-knot"
//                                    algorithm
//    19 Jun. 2009, V.Ivanchenko : removed hidden bin 
//    17 Nov. 2009, H.Kurashige   : use pointer for DataVector
//    04 May  2010  H.Kurashige   : use G4PhyscisVectorCache
//    28 May  2010  H.Kurashige  : Stop using  pointers to G4PVDataVector
//    16 Aug. 2011  H.Kurashige  : Add dBin, baseBin and verboseLevel
// --------------------------------------------------------------

#include <iomanip>
#include "wuPhysicsVector.hh"

// --------------------------------------------------------------

wuPhysicsVector::wuPhysicsVector(bool)
 : type(T_wuPhysicsVector),
   edgeMin(0.), edgeMax(0.), numberOfNodes(0),
   useSpline(false), 
   dBin(0.), baseBin(0.),
   verboseLevel(0)
{

}

// --------------------------------------------------------------

wuPhysicsVector::~wuPhysicsVector() 
{
}

// --------------------------------------------------------------

wuPhysicsVector::wuPhysicsVector(const wuPhysicsVector& right)
{
  dBin         = right.dBin;
  baseBin      = right.baseBin;
  verboseLevel = right.verboseLevel;

  DeleteData();
  CopyData(right);
}

// --------------------------------------------------------------

wuPhysicsVector& wuPhysicsVector::operator=(const wuPhysicsVector& right)
{
  if (&right==this)  { return *this; }
  dBin         = right.dBin;
  baseBin      = right.baseBin;
  verboseLevel = right.verboseLevel;

  DeleteData();
  CopyData(right);
  return *this;
}

// --------------------------------------------------------------

int wuPhysicsVector::operator==(const wuPhysicsVector &right) const
{
  return (this == &right);
}

// --------------------------------------------------------------

int wuPhysicsVector::operator!=(const wuPhysicsVector &right) const
{
  return (this != &right);
}

// --------------------------------------------------------------

void wuPhysicsVector::DeleteData()
{
  useSpline = false;
  secDerivative.clear();
}

// --------------------------------------------------------------

void wuPhysicsVector::CopyData(const wuPhysicsVector& vec)
{
  type = vec.type;
  edgeMin = vec.edgeMin;
  edgeMax = vec.edgeMax;
  numberOfNodes = vec.numberOfNodes;
  useSpline = vec.useSpline;

  size_t i;
  dataVector.resize(numberOfNodes);
  for(i=0; i<numberOfNodes; ++i) { 
    dataVector[i] = (vec.dataVector)[i];
  }
  binVector.resize(numberOfNodes);
  for(i=0; i<numberOfNodes; ++i) { 
    binVector[i] = (vec.binVector)[i];
  }
  if(0 < (vec.secDerivative).size()) {
    secDerivative.resize(numberOfNodes);
    for(i=0; i<numberOfNodes; ++i){ 
      secDerivative[i] = (vec.secDerivative)[i];
    }
  }
}

// --------------------------------------------------------------

double wuPhysicsVector::GetLowEdgeEnergy(size_t binNumber) const
{
  return binVector[binNumber];
}

// --------------------------------------------------------------

bool wuPhysicsVector::Store(std::ofstream& fOut, bool ascii)
{
  // Ascii mode
  if (ascii)
  {
    fOut << *this;
    return true;
  } 
  // Binary Mode

  // binning
  fOut.write((char*)(&edgeMin), sizeof edgeMin);
  fOut.write((char*)(&edgeMax), sizeof edgeMax);
  fOut.write((char*)(&numberOfNodes), sizeof numberOfNodes);

  // contents
  size_t size = dataVector.size(); 
  fOut.write((char*)(&size), sizeof size);

  double* value = new double[2*size];
  for(size_t i = 0; i < size; ++i)
  {
    value[2*i]  =  binVector[i];
    value[2*i+1]=  dataVector[i];
  }
  fOut.write((char*)(value), 2*size*(sizeof (double)));
  delete [] value;

  return true;
}

// --------------------------------------------------------------

bool wuPhysicsVector::Retrieve(std::ifstream& fIn, bool ascii)
{
  // clear properties;
  dataVector.clear();
  binVector.clear();
  secDerivative.clear();

  // retrieve in ascii mode
  if (ascii){
    // binning
    fIn >> edgeMin >> edgeMax >> numberOfNodes; 
    if (fIn.fail())  { return false; }
    // contents
    int siz=0;
    fIn >> siz;
    if (fIn.fail())  { return false; }
    if (siz<=0)
    {
// #ifdef G4VERBOSE  
//       G4cerr << "G4PhysicsVector::Retrieve():";
//       G4cerr << " Invalid vector size: " << siz << G4endl;
// #endif
      std::cout << "wuPhysicsVector::Retrieve():";
      std::cout << " Invalid vector size: " << siz << std::endl;
      return false;
    }

    binVector.reserve(siz);
    dataVector.reserve(siz);
    double vBin, vData;

    for(int i = 0; i < siz ; i++)
    {
      vBin = 0.;
      vData= 0.;
      fIn >> vBin >> vData;
      if (fIn.fail())  { return false; }
      binVector.push_back(vBin);
      dataVector.push_back(vData);
    }

    // to remove any inconsistency 
    numberOfNodes = siz;
    edgeMin = binVector[0];
    edgeMax = binVector[numberOfNodes-1];
    return true ;
  }

  // retrieve in binary mode
  // binning
  fIn.read((char*)(&edgeMin), sizeof edgeMin);
  fIn.read((char*)(&edgeMax), sizeof edgeMax);
  fIn.read((char*)(&numberOfNodes), sizeof numberOfNodes ); 
 
  // contents
  size_t size;
  fIn.read((char*)(&size), sizeof size); 
 
  double* value = new double[2*size];
  fIn.read((char*)(value),  2*size*(sizeof(double)) );
  if (int(fIn.gcount()) != int(2*size*(sizeof(double))) )
  {
    delete [] value;
    return false;
  }

  binVector.reserve(size);
  dataVector.reserve(size);
  for(size_t i = 0; i < size; ++i)
  {
    binVector.push_back(value[2*i]);
    dataVector.push_back(value[2*i+1]);
  }
  delete [] value;

  // to remove any inconsistency 
  numberOfNodes = size;
  edgeMin = binVector[0];
  edgeMax = binVector[numberOfNodes-1];

  return true;
}

// --------------------------------------------------------------

void 
wuPhysicsVector::ScaleVector(double factorE, double factorV)
{
  size_t n = dataVector.size();
  size_t i;
  if(n > 0) { 
    for(i=0; i<n; ++i) {
      binVector[i]  *= factorE;
      dataVector[i] *= factorV;
    } 
  }
  //  n = secDerivative.size();
  // if(n > 0) { for(i=0; i<n; ++i) { secDerivative[i] *= factorV; } }
  secDerivative.clear();

  edgeMin *= factorE;
  edgeMax *= factorE;
}

// --------------------------------------------------------------

void 
wuPhysicsVector::ComputeSecondDerivatives(double firstPointDerivative, 
                                          double endPointDerivative)
  //  A standard method of computation of second derivatives 
  //  First derivatives at the first and the last point should be provided
  //  See for example W.H. Press et al. "Numerical recipes in C"
  //  Cambridge University Press, 1997.
{
  if(4 > numberOfNodes)   // cannot compute derivatives for less than 4 bins
  {
    ComputeSecDerivatives();
    return;
  }

  if(!SplinePossible()) { return; }

  useSpline = true;

  int n = numberOfNodes-1;

  double* u = new double [n];
  
  double p, sig, un;

  u[0] = (6.0/(binVector[1]-binVector[0]))
    * ((dataVector[1]-dataVector[0])/(binVector[1]-binVector[0])
       - firstPointDerivative);
 
  secDerivative[0] = - 0.5;

  // Decomposition loop for tridiagonal algorithm. secDerivative[i]
  // and u[i] are used for temporary storage of the decomposed factors.

  for(int i=1; i<n; ++i)
  {
    sig = (binVector[i]-binVector[i-1]) / (binVector[i+1]-binVector[i-1]);
    p = sig*(secDerivative[i-1]) + 2.0;
    secDerivative[i] = (sig - 1.0)/p;
    u[i] = (dataVector[i+1]-dataVector[i])/(binVector[i+1]-binVector[i])
         - (dataVector[i]-dataVector[i-1])/(binVector[i]-binVector[i-1]);
    u[i] = 6.0*u[i]/(binVector[i+1]-binVector[i-1]) - sig*u[i-1]/p;
  }

  sig = (binVector[n-1]-binVector[n-2]) / (binVector[n]-binVector[n-2]);
  p = sig*secDerivative[n-2] + 2.0;
  un = (6.0/(binVector[n]-binVector[n-1]))
    *(endPointDerivative - 
      (dataVector[n]-dataVector[n-1])/(binVector[n]-binVector[n-1])) - u[n-1]/p;
  secDerivative[n] = un/(secDerivative[n-1] + 2.0);

  // The back-substitution loop for the triagonal algorithm of solving
  // a linear system of equations.
   
  for(int k=n-1; k>0; --k)
  {
    secDerivative[k] *= 
      (secDerivative[k+1] - 
       u[k]*(binVector[k+1]-binVector[k-1])/(binVector[k+1]-binVector[k]));
  }
  secDerivative[0] = 0.5*(u[0] - secDerivative[1]);

  delete [] u;
}

// --------------------------------------------------------------

void wuPhysicsVector::FillSecondDerivatives()
  // Computation of second derivatives using "Not-a-knot" endpoint conditions
  // B.I. Kvasov "Methods of shape-preserving spline approximation"
  // World Scientific, 2000
{
  if(5 > numberOfNodes)  // cannot compute derivatives for less than 4 points
  {
    ComputeSecDerivatives();
    return;
  }

  if(!SplinePossible()) { return; }

  useSpline = true;
 
  int n = numberOfNodes-1;

  double* u = new double [n];
  
  double p, sig;

  u[1] = ((dataVector[2]-dataVector[1])/(binVector[2]-binVector[1]) -
          (dataVector[1]-dataVector[0])/(binVector[1]-binVector[0]));
  u[1] = 6.0*u[1]*(binVector[2]-binVector[1])
    / ((binVector[2]-binVector[0])*(binVector[2]-binVector[0]));
 
  // Decomposition loop for tridiagonal algorithm. secDerivative[i]
  // and u[i] are used for temporary storage of the decomposed factors.

  secDerivative[1] = (2.0*binVector[1]-binVector[0]-binVector[2])
    / (2.0*binVector[2]-binVector[0]-binVector[1]);

  for(int i=2; i<n-1; ++i)
  {
    sig = (binVector[i]-binVector[i-1]) / (binVector[i+1]-binVector[i-1]);
    p = sig*secDerivative[i-1] + 2.0;
    secDerivative[i] = (sig - 1.0)/p;
    u[i] = (dataVector[i+1]-dataVector[i])/(binVector[i+1]-binVector[i])
         - (dataVector[i]-dataVector[i-1])/(binVector[i]-binVector[i-1]);
    u[i] = (6.0*u[i]/(binVector[i+1]-binVector[i-1])) - sig*u[i-1]/p;
  }

  sig = (binVector[n-1]-binVector[n-2]) / (binVector[n]-binVector[n-2]);
  p = sig*secDerivative[n-3] + 2.0;
  u[n-1] = (dataVector[n]-dataVector[n-1])/(binVector[n]-binVector[n-1])
    - (dataVector[n-1]-dataVector[n-2])/(binVector[n-1]-binVector[n-2]);
  u[n-1] = 6.0*sig*u[n-1]/(binVector[n]-binVector[n-2])
    - (2.0*sig - 1.0)*u[n-2]/p;  

  p = (1.0+sig) + (2.0*sig-1.0)*secDerivative[n-2];
  secDerivative[n-1] = u[n-1]/p;

  // The back-substitution loop for the triagonal algorithm of solving
  // a linear system of equations.
   
  for(int k=n-2; k>1; --k)
  {
    secDerivative[k] *= 
      (secDerivative[k+1] - 
       u[k]*(binVector[k+1]-binVector[k-1])/(binVector[k+1]-binVector[k]));
  }
  secDerivative[n] = (secDerivative[n-1] - (1.0-sig)*secDerivative[n-2])/sig;
  sig = 1.0 - ((binVector[2]-binVector[1])/(binVector[2]-binVector[0]));
  secDerivative[1] *= (secDerivative[2] - u[1]/(1.0-sig));
  secDerivative[0]  = (secDerivative[1] - sig*secDerivative[2])/(1.0-sig);

  delete [] u;
}

// --------------------------------------------------------------

void 
wuPhysicsVector::ComputeSecDerivatives()
  //  A simplified method of computation of second derivatives 
{
  if(3 > numberOfNodes)  // cannot compute derivatives for less than 4 bins
  {
    useSpline = false;
    return;
  }

  if(!SplinePossible())  { return; }

  useSpline = true;

  size_t n = numberOfNodes-1;

  for(size_t i=1; i<n; ++i)
  {
    secDerivative[i] =
      3.0*((dataVector[i+1]-dataVector[i])/(binVector[i+1]-binVector[i]) -
           (dataVector[i]-dataVector[i-1])/(binVector[i]-binVector[i-1]))
      /(binVector[i+1]-binVector[i-1]);
  }
  secDerivative[n] = secDerivative[n-1];
  secDerivative[0] = secDerivative[1];
}

// --------------------------------------------------------------

bool wuPhysicsVector::SplinePossible()
  // Initialise second derivative array. If neighbor energy coincide 
  // or not ordered than spline cannot be applied
{
  bool result = true;
  for(size_t j=1; j<numberOfNodes; ++j)
  {
    if(binVector[j] <= binVector[j-1])  { 
      result = false; 
      useSpline = false;
      secDerivative.clear();
      break;
    }
  }  
  secDerivative.resize(numberOfNodes,0.0);
  return result;
}
   
// --------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const wuPhysicsVector& pv)
{
  // binning
  out << std::setprecision(12) << pv.edgeMin << " "
      << pv.edgeMax << " " << pv.numberOfNodes << std::endl; 

  // contents
  out << pv.dataVector.size() << std::endl; 
  for(size_t i = 0; i < pv.dataVector.size(); i++)
  {
    out << pv.binVector[i] << "  " << pv.dataVector[i] << std::endl;
  }
  out << std::setprecision(6);

  return out;
}

//---------------------------------------------------------------

double wuPhysicsVector::Value(double theEnergy, size_t& lastIdx) const
{
  double y;
  if(theEnergy <= edgeMin) {
    lastIdx = 0; 
    y = dataVector[0]; 
  } else if(theEnergy >= edgeMax) { 
    lastIdx = numberOfNodes-1; 
    y = dataVector[lastIdx]; 
  } else {
    lastIdx = FindBin(theEnergy, lastIdx);
    y = Interpolation(lastIdx, theEnergy);
  }
  return y;
}

//---------------------------------------------------------------

double wuPhysicsVector::FindLinearEnergy(double rand) const
{
  if(1 >= numberOfNodes) { return 0.0; }
  size_t n1 = 0;
  size_t n2 = numberOfNodes/2;
  size_t n3 = numberOfNodes - 1;
  double y = rand*dataVector[n3];
  while (n1 + 1 != n3)
    {
      if (y > dataVector[n2])
	{ n1 = n2; }
      else
	{ n3 = n2; }
      n2 = (n3 + n1 + 1)/2;
    }
  double res = binVector[n1];
  double del = dataVector[n3] - dataVector[n1];
  if(del > 0.0) { 
    res += (y - dataVector[n1])*(binVector[n3] - res)/del;  
  }
  return res;
}

//---------------------------------------------------------------

// 
// wuPhysicsVector.cc ends here
