// wuPhysicsOrderedFreeVector.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 四 2月  5 20:50:21 2015 (+0800)
// Last-Updated: 四 2月  5 20:58:46 2015 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 2
// URL: http://wuhongyi.cn 

////////////////////////////////////////////////////////////////////////
// PhysicsOrderedFreeVector Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4PhysicsOrderedFreeVector.cc
// Version:     2.0
// Created:     1996-08-13
// Author:      Juliet Armstrong
// Updated:     1997-03-25 by Peter Gumplinger
//              > cosmetics (only)
//              1998-11-11 by Peter Gumplinger
//              > initialize all data members of the base class in 
//                derived class constructors
//              2000-11-11 by H.Kurashige
//              > use STL vector for dataVector and binVector
//              2009-06-19 by V.Ivanchenko 
//              > removed hidden bin 
//              2013-10-02 by V.Ivanchenko removed FindBinLocation   
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "wuPhysicsOrderedFreeVector.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        /////////////////
        // Constructors
        /////////////////

wuPhysicsOrderedFreeVector::wuPhysicsOrderedFreeVector(double *Energies,
                                                       double *Values,
                                                       size_t VectorLength)
  : wuPhysicsVector()
{
  type = T_wuPhysicsOrderedFreeVector;

  for (size_t i = 0 ; i < VectorLength ; i++)
    {
      InsertValues(Energies[i], Values[i]);
    }
}

wuPhysicsOrderedFreeVector::wuPhysicsOrderedFreeVector()
  : wuPhysicsVector()
{
  type = T_wuPhysicsOrderedFreeVector;
}

        ////////////////
        // Destructors
        ////////////////

wuPhysicsOrderedFreeVector::~wuPhysicsOrderedFreeVector() {}

        ////////////
        // Methods
        ////////////
  
void wuPhysicsOrderedFreeVector::InsertValues(double energy, double value)
{
        std::vector<double>::iterator binLoc =
                 std::lower_bound(binVector.begin(), binVector.end(), energy);

        size_t binIdx = binLoc - binVector.begin();	// Iterator difference!

        std::vector<double>::iterator dataLoc = dataVector.begin() + binIdx;

        binVector.insert(binLoc, energy);
        dataVector.insert(dataLoc, value);

        numberOfNodes++;
        edgeMin = binVector.front();
        edgeMax = binVector.back();
}

double  wuPhysicsOrderedFreeVector::GetEnergy(double aValue)
{

        if (aValue <= GetMinValue()) {
                return GetMinLowEdgeEnergy();
        } else if (aValue >= GetMaxValue()) {
                return GetMaxLowEdgeEnergy();
        } else { 
        size_t closestBin = FindValueBinLocation(aValue);
        double theEnergy = LinearInterpolationOfEnergy(aValue, closestBin);

        return theEnergy;
        }
}

size_t wuPhysicsOrderedFreeVector::FindValueBinLocation(double aValue)
{
   int n1 = 0;
   int n2 = numberOfNodes/2;
   int n3 = numberOfNodes - 1;
   while (n1 != n3 - 1) {
      if (aValue > dataVector[n2])
         { n1 = n2; }
      else
         { n3 = n2; }
      n2 = n1 + (n3 - n1 + 1)/2;
   }
   return (size_t)n1;
}

double wuPhysicsOrderedFreeVector::LinearInterpolationOfEnergy(double aValue,
								 size_t theLocBin)
{
  double intplFactor = (aValue-dataVector[theLocBin])
     / (dataVector[theLocBin+1]-dataVector[theLocBin]); // Interpolation factor

  return binVector[theLocBin] +
         ( binVector[theLocBin+1]-binVector[theLocBin] ) * intplFactor;
}


// 
// wuPhysicsOrderedFreeVector.cc ends here
