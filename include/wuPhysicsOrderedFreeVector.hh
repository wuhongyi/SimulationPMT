// wuPhysicsOrderedFreeVector.hh --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 四 2月  5 20:50:09 2015 (+0800)
// Last-Updated: 四 2月  5 20:58:47 2015 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 2
// URL: http://wuhongyi.cn 

////////////////////////////////////////////////////////////////////////
// PhysicsOrderedFreeVector Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:	G4PhysicsOrderedFreeVector.hh
// Created:     1996-08-13
// Author:      Juliet Armstrong
// Updated:     1997-03-25 by Peter Gumplinger
//		> cosmetics (only)
//              2000-11-11 by H.Kurashige
//              > use STL vector for dataVector and binVector
// mail:        gum@triumf.ca
//
// Class description:
//
//    A physics ordered free vector inherits from G4PhysicsVector which
//    has values of energy-loss, cross-section, and other physics values
//    of a particle in matter in a given range of the energy, momentum,
//    etc.). In addition, the ordered free vector provides a method for
//    the user to insert energy/value pairs in sequence.  Methods to
//    Retrieve the Max and Min energies and values from the vector are
//    also provided. 

////////////////////////////////////////////////////////////////////////

#ifndef _WUPHYSICSORDEREDFREEVECTOR_H_
#define _WUPHYSICSORDEREDFREEVECTOR_H_

/////////////
// Includes
/////////////

#include "wuPhysicsVector.hh"

/////////////////////
// Class Definition
/////////////////////

class wuPhysicsOrderedFreeVector : public wuPhysicsVector 
{

 public: // with description
	
  ////////////////////////////////
  // Constructors and Destructor
  ////////////////////////////////
  
  wuPhysicsOrderedFreeVector();
  wuPhysicsOrderedFreeVector(double* Energies,
			     double* Values,
			     size_t VectorLength);
  
  virtual ~wuPhysicsOrderedFreeVector();
  
  ////////////
  // Methods
  ////////////
  
  void InsertValues(double energy, double value); 

  //double GetLowEdgeEnergy(size_t binNumber) const;
  
  double GetMaxValue();
  
  double GetMinValue();
  
  double GetEnergy(double aValue);
  
  double GetMaxLowEdgeEnergy();
  
  double GetMinLowEdgeEnergy();
  
  void DumpValues();

 private:
  
  size_t FindValueBinLocation(double aValue);
  
  double LinearInterpolationOfEnergy(double aValue, size_t theLocBin);
};

////////////////////
// Inline methods
////////////////////

inline
double wuPhysicsOrderedFreeVector::GetMaxValue()
{
	return dataVector.back();
}

inline
double wuPhysicsOrderedFreeVector::GetMinValue()
{
	return dataVector.front();
}

inline
double wuPhysicsOrderedFreeVector::GetMaxLowEdgeEnergy()
{
	return binVector.back();
}

inline
double wuPhysicsOrderedFreeVector::GetMinLowEdgeEnergy()
{
	return binVector.front();
}

inline
void wuPhysicsOrderedFreeVector::DumpValues()
{
   for (size_t i = 0; i < numberOfNodes; i++)
   {
     std::cout << binVector[i] << "\t" << dataVector[i] << std::endl;
   }
}

#endif /* _WUPHYSICSORDEREDFREEVECTOR_H_ */

// 
// wuPhysicsOrderedFreeVector.hh ends here
