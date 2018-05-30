// wuPhysicsVector.hh --- 
// 
// Description: Coping from Geant4 G4PhysicsVector.hh
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 四 2月  5 20:20:24 2015 (+0800)
// Last-Updated: 四 2月  5 20:55:07 2015 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 4
// URL: http://wuhongyi.cn 

//---------------------------------------------------------------
//      GEANT 4 class header file
//
//  G4PhysicsVector.hh
//
//  Class description:
//
//    A physics vector which has values of energy-loss, cross-section, 
//    and other physics values of a particle in matter in a given 
//    range of the energy, momentum, etc.
//    This class serves as the base class for a vector having various 
//    energy scale, for example like 'log', 'linear', 'free', etc.

//  History:
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    03 Mar. 1996, K.Amako : Implemented the 1st version
//    27 Apr. 1996, K.Amako : Cache mechanism added
//    01 Jul. 1996, K.Amako : Now GetValue not virtual
//    21 Sep. 1996, K.Amako : Added [] and () operators
//    11 Nov. 2000, H.Kurashige : Use STL vector for dataVector and binVector
//    09 Mar. 2001, H.Kurashige : Added G4PhysicsVectorType & Store/Retrieve()
//    02 Apr. 2008, A.Bagulya : Added SplineInterpolation() and SetSpline()
//    11 May  2009, V.Ivanchenko : Added ComputeSecondDerivatives
//    19 Jun. 2009, V.Ivanchenko : Removed hidden bin 
//    22 Dec. 2009  H.Kurashige  : Use pointers to G4PVDataVector
//    04 May. 2010  H.Kurashige  : Use G4PhysicsVectorCache
//    28 May  2010  H.Kurashige  : Stop using  pointers to G4PVDataVector
//    16 Aug. 2011  H.Kurashige  : Add dBin, baseBin and verboseLevel
//    02 Oct. 2013  V.Ivanchenko : FindBinLocation method become inlined;
//                                 instead of G4Pow G4Log is used
//---------------------------------------------------------------

#ifndef _WUPHYSICSVECTOR_H_
#define _WUPHYSICSVECTOR_H_

enum wuPhysicsVectorType
{
  T_wuPhysicsVector =0,
  T_wuPhysicsLinearVector,
  T_wuPhysicsLogVector,
  T_wuPhysicsLnVector,
  T_wuPhysicsFreeVector,
  T_wuPhysicsOrderedFreeVector,
  T_wuLPhysicsFreeVector
};

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

class wuPhysicsVector
{
public: // with description
  wuPhysicsVector(bool spline = false);
  // constructor  
  // This class is an abstract class with pure virtual method of
  // virtual size_t FindBinLocation(double theEnergy) const
  // So, default constructor is not supposed to be invoked explicitly

  wuPhysicsVector(const wuPhysicsVector&);
  wuPhysicsVector& operator=(const wuPhysicsVector&);
  // Copy constructor and assignment operator.

  virtual ~wuPhysicsVector();
  // destructor

    double Value(double theEnergy, size_t& lastidx) const; 
  // Get the cross-section/energy-loss value corresponding to the
  // given energy. An appropriate interpolation is used to calculate
  // the value. Consumer code got changed index and may reuse it
  // for the next call to save CPU for bin location. 


    inline double Value(double theEnergy) const; 
         // Get the cross-section/energy-loss value corresponding to the
         // given energy. An appropriate interpolation is used to calculate
         // the value. This method is kept for backward compatibility reason,
         // it should be used is bin location cannot be kept thread safe

    inline double GetValue(double theEnergy, bool& isOutRange) const;
         // Obsolete method to get value, isOutRange is not used anymore. 
         // This method is kept for the compatibility reason.

    int operator==(const wuPhysicsVector &right) const ;
    int operator!=(const wuPhysicsVector &right) const ;

    inline double operator[](const size_t binNumber) const ;
         // Returns simply the value in the bin specified by 'binNumber'
         // of the dataVector. The boundary check will not be done. 

    inline double operator()(const size_t binNumber) const ;
         // Returns simply the value in the bin specified by 'binNumber'
         // of the dataVector. The boundary check will not be Done. 

    inline void PutValue(size_t index, double theValue);
         // Put 'theValue' into the bin specified by 'binNumber'.
         // Take note that the 'index' starts from '0'.
         // To fill the vector, you have beforehand to construct a vector
         // by the constructor with Emin, Emax, Nbin. 'theValue' should
         // be the crosssection/energyloss value corresponding to the  
         // energy of the index. You can get this energy by the next method
         // or by the old method GetLowEdgeEnergy().

    virtual void ScaleVector(double factorE, double factorV);
         // Scale all values of the vector and second derivatives
         // by factorV, energies by vectorE. This method may be applied 
         // for example after Retrieve a vector from an external file to 
         // convert values into Geant4 units

    inline double Energy(size_t index) const;
         // Returns simply the value in the energy specified by 'index'
         // of the energy vector. The boundary check will not be done. 
         // Use this function when you fill physis vector by PutValue().

    inline double GetMaxEnergy() const;
         // Returns the energy of the last point of the vector

    double GetLowEdgeEnergy(size_t binNumber) const;
         // Obsolete method
         // Get the energy value at the low edge of the specified bin.
         // Take note that the 'binNumber' starts from '0'.
         // This value should be defined before the call.
         // The boundary check will not be done.

    inline size_t GetVectorLength() const;
         // Get the total length of the vector. 

    inline size_t FindBin(double energy, size_t idx) const;
         // find low edge index of a bin for given energy
         // min value 0, max value VectorLength-1
         // idx is suggested bin number from user code

    void FillSecondDerivatives();
        // Initialise second derivatives for spline keeping 
        // 3d derivative continues - default algorithm

    void ComputeSecDerivatives();
         // Initialise second derivatives for spline using algorithm 
         // which garantee only 1st derivative continues 
         // Warning: this method should be called when the vector 
         // is already filled

    void ComputeSecondDerivatives(double firstPointDerivative, 
                                  double endPointDerivative);
         // Initialise second derivatives for spline using 
         // user defined 1st derivatives at edge points
         // Warning: this method should be called when the vector 
         // is already filled

    double FindLinearEnergy(double rand) const;
         // Find energy using linear interpolation for vector
         // filled by cumulative probability function 
         // value of rand should be between 0 and 1

    inline bool IsFilledVectorExist() const;
         // Is non-empty physics vector already exist?

    inline wuPhysicsVectorType GetType() const;
         // Get physics vector type
  
    inline void SetSpline(bool);
         // Activate/deactivate Spline interpolation.

    virtual bool Store(std::ofstream& fOut, bool ascii=false);
    virtual bool Retrieve(std::ifstream& fIn, bool ascii=false);
         // To store/retrieve persistent data to/from file streams.

    friend std::ostream& operator<<(std::ostream&, const wuPhysicsVector&);

    inline void SetVerboseLevel(int value);
    inline int GetVerboseLevel(int);
         // Set/Get Verbose level

  protected:

    void DeleteData();
    void CopyData(const wuPhysicsVector& vec);
         // Internal methods for allowing copy of objects

  protected:

    wuPhysicsVectorType type;   // The type of PhysicsVector (enumerator)

    double edgeMin;           // Energy of first point
    double edgeMax;           // Energy of the last point

    size_t numberOfNodes;

    std::vector<double>  dataVector;    // Vector to keep the crossection/energyloss
    std::vector<double>  binVector;     // Vector to keep energy
    std::vector<double>  secDerivative; // Vector to keep second derivatives 

  private:

    bool SplinePossible();

    inline double LinearInterpolation(size_t idx, double energy) const;
         // Linear interpolation function
    inline double SplineInterpolation(size_t idx, double energy) const;
         // Spline interpolation function

    inline double Interpolation(size_t idx, double energy) const;

    inline size_t FindBinLocation(double theEnergy) const;
         // Find the bin# in which theEnergy belongs 

    bool     useSpline;

  protected:

    double dBin;          // Bin width - useful only for fixed binning
    double baseBin;       // Set this in constructor for performance

    int verboseLevel;

};


inline
 double wuPhysicsVector::operator[](const size_t binNumber) const
{
  return  dataVector[binNumber];
}

//---------------------------------------------------------------

inline
 double wuPhysicsVector::operator()(const size_t binNumber) const
{
  return dataVector[binNumber];
}

//---------------------------------------------------------------

inline
 double wuPhysicsVector::Energy(const size_t binNumber) const
{
  return binVector[binNumber];
}

//---------------------------------------------------------------

inline
 double wuPhysicsVector::GetMaxEnergy() const
{
  return edgeMax;
}

//---------------------------------------------------------------

inline 
 size_t wuPhysicsVector::GetVectorLength() const
{
  return numberOfNodes;
}

//---------------------------------------------------------------

inline
 double wuPhysicsVector::GetValue(double theEnergy, bool&) const 
{
  size_t idx=0;
  return Value(theEnergy, idx);
}

//------------------------------------------------

inline
 double wuPhysicsVector::LinearInterpolation(size_t idx, double e) const
{
  // Linear interpolation is used to get the value. Before this method
  // is called it is ensured that the energy is inside the bin
  // 0 < idx < numberOfNodes-1
  
  return dataVector[idx] +
         ( dataVector[idx + 1]-dataVector[idx] ) * (e - binVector[idx])
         /( binVector[idx + 1]-binVector[idx] );
}

//---------------------------------------------------------------

inline
 double wuPhysicsVector::SplineInterpolation(size_t idx, double e) const
{
  // Spline interpolation is used to get the value. Before this method
  // is called it is ensured that the energy is inside the bin
  // 0 < idx < numberOfNodes-1

  static const double onesixth = 1.0/6.0;

  // check bin value
  double x1 = binVector[idx];
  double x2 = binVector[idx + 1];
  double delta = x2 - x1;

  double a = (x2 - e)/delta;
  double b = (e - x1)/delta;
   
  // Final evaluation of cubic spline polynomial for return   
  double y1 = dataVector[idx];
  double y2 = dataVector[idx + 1];

  double res = a*y1 + b*y2 + 
        ( (a*a*a - a)*secDerivative[idx] +
          (b*b*b - b)*secDerivative[idx + 1] )*delta*delta*onesixth;

  return res;
}

//---------------------------------------------------------------

inline 
 double wuPhysicsVector::Interpolation(size_t idx, double e) const
{
  double res;
  if(useSpline) { res = SplineInterpolation(idx, e); }
  else          { res = LinearInterpolation(idx, e); }
  return res;
}

//---------------------------------------------------------------

inline 
 void wuPhysicsVector::PutValue(size_t binNumber, double theValue)
{
  dataVector[binNumber] = theValue;
}

//---------------------------------------------------------------

inline 
 bool wuPhysicsVector::IsFilledVectorExist() const
{
  bool status=false;

  if(numberOfNodes > 0) { status=true; }
  return status;
}

//---------------------------------------------------------------

inline 
 wuPhysicsVectorType wuPhysicsVector::GetType() const
{
  return type;
}

//---------------------------------------------------------------

// Flag useSpline is "true" only if second derivatives are filled 
inline 
 void wuPhysicsVector::SetSpline(bool val)
{
  if(val) {
    if(0 == secDerivative.size() && 0 < dataVector.size()) { 
      FillSecondDerivatives(); 
    }
  } else {
    useSpline = false;
    secDerivative.clear();
  }
}

//---------------------------------------------------------------

inline
void wuPhysicsVector::SetVerboseLevel(int value)
{
   verboseLevel = value;
}

//---------------------------------------------------------------

inline
int wuPhysicsVector::GetVerboseLevel(int)
{
   return verboseLevel;
}

//---------------------------------------------------------------

inline 
size_t wuPhysicsVector::FindBinLocation(double theEnergy) const
{
   size_t bin;
   if(type == T_wuPhysicsLogVector) {
     bin = size_t(std::log(theEnergy)/dBin - baseBin);
     if(bin + 2 > numberOfNodes) { bin = numberOfNodes - 2; }
     else if(bin > 0 && theEnergy < binVector[bin]) { --bin; }
     else if(bin + 2 < numberOfNodes && theEnergy > binVector[bin+1]) 
            { ++bin; }
   } else if(type == T_wuPhysicsLinearVector) {
     bin = size_t( theEnergy/dBin - baseBin ); 
     if(bin + 2 > numberOfNodes) { bin = numberOfNodes - 2; }   
     else if(bin > 0 && theEnergy < binVector[bin]) { --bin; }
     else if(bin + 2 < numberOfNodes && theEnergy > binVector[bin+1]) 
            { ++bin; }
   } else {
     bin = 0;
     size_t bin2;
     size_t bin3 = numberOfNodes - 1;
     while (bin != bin3 - 1) {
       bin2 = bin + (bin3 - bin + 1)/2;
       if (theEnergy > binVector[bin2]) { bin = bin2; }
       else { bin3 = bin2; }
     }
/*
// Bin location proposed by K.Genser (FNAL)
// V.I. Usage of this algorithm provides identical results
// no CPU advantage is observed in EM tests
// if new validation information will be known this code may be used
     PVDataVector::const_iterator it =
       std::lower_bound(binVector.begin(), binVector.end(), theEnergy);
     bin = it - binVector.begin() - 1;
*/
   }
   return bin;
}

//---------------------------------------------------------------

inline size_t wuPhysicsVector::FindBin(double e, size_t idx) const
{
  size_t id = idx;
  if(e < binVector[1]) { 
    id = 0; 
  } else if(e >= binVector[numberOfNodes-2]) { 
    id = numberOfNodes - 2; 
  } else if(idx >= numberOfNodes || e < binVector[idx] 
            || e > binVector[idx+1]) { 
    id = FindBinLocation(e); 
  }
  return id;
}

//---------------------------------------------------------------

inline
 double wuPhysicsVector::Value(double theEnergy) const
{
  size_t idx=0;
  return Value(theEnergy, idx);
}

//---------------------------------------------------------------






#endif /* _WUPHYSICSVECTOR_H_ */

// 
// wuPhysicsVector.hh ends here
