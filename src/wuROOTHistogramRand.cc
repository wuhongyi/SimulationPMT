// wuROOTHistogramRand.cc --- 
// 
// Description: 保持与wuHistogramRand.cc一致
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 六 2月 20 13:39:59 2016 (+0800)
// Last-Updated: 六 2月 20 13:44:15 2016 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 2
// URL: http://wuhongyi.cn 

#include "wuROOTHistogramRand.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
wuROOTHistogramRand::wuROOTHistogramRand()
{

}

wuROOTHistogramRand::~wuROOTHistogramRand()
{

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void wuROOTHistogramRand::SetHistogram(const std::vector<double> *x,const std::vector<double> *y, const int n)
{
  N=n;
  X.clear();
  Y.clear();
  for (int i = 0; i < N; ++i)
    {
      X.push_back(x->at(i));
      Y.push_back(y->at(i));
    }

  if(N==1) return;
  double count=0;
  for (int i = 0; i < N-1; ++i)
    {
      count+=((Y[i+1]+Y[i])*(X[i+1]-X[i])/2);
    }
  sumweight.clear();
  sumweight.push_back(0);
  for (int i = 0; i < N-1; ++i)
    {
      sumweight.push_back(((Y[i+1]+Y[i])*(X[i+1]-X[i])/2)/count);
    }
  double temp=0;
  for (int i = 0; i < N; ++i)
    {
      temp=sumweight.at(i)+temp;
      sumweight.at(i)=temp;
    }
}

double wuROOTHistogramRand::GetRand()
{
  // if(N==1) return HistogramWeight[0].first;//单能
  if(N==1) return X[0];//单能
  double rand = gRandom->Rndm();
  for (int i = 0; i < N; ++i)
    {
      if(rand <= sumweight[i])
	{
	  return GetTrapezoidX(X[i-1],Y[i-1],sumweight[i-1],X[i],Y[i],sumweight[i],rand);
	}
    }

}

double wuROOTHistogramRand::GetTrapezoidX(double x1,double y1,double eta1,double x2,double y2,double eta2,double eta)
{
  //能谱抽样中，对每相邻的两个点构成的梯形，eta2-eta1为其在x1到x2之间的概率。返回概率eta所对应的x
  if(std::abs(y2-y1)<1.0e-6)// ==
    {
      return x1+(eta-eta1)*(x2-x1)/(eta2-eta1);
    }
  else// !=
    {
      double k=(y2-y1)/(x2-x1);
      double b=(eta-eta1)*(y2+y1)*(x2-x1)/(eta2-eta1);
      return x1+(-y1+std::sqrt(y1*y1+k*b))/k;
    }
}


// 
// wuROOTHistogramRand.cc ends here










