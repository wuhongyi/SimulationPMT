// wuROOTHistogramRand.hh --- 
// 
// Description: 保持与wuHistogramRand.hh同步
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 六 2月 20 13:35:45 2016 (+0800)
// Last-Updated: 六 2月 20 13:39:01 2016 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 1
// URL: http://wuhongyi.cn 

//ROOT使用
//采用ROOT随机数

#ifndef _WUROOTHISTOGRAMRAND_H_
#define _WUROOTHISTOGRAMRAND_H_

#include <map>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>

#include "TRandom.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class wuROOTHistogramRand
{
public:
  wuROOTHistogramRand();
  virtual ~wuROOTHistogramRand();

public:
  void SetHistogram(const std::vector<double> *x,const std::vector<double> *y, const int n);
  double GetRand();
  double GetTrapezoidX(double x1,double y1,double eta1,double x2,double y2,double eta2,double eta);

private:
  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> sumweight; 
  // std::map<double,double>::iterator it;
  int N;//==1 时为单能，>1 时为谱分布
};



#endif /* _WUROOTHISTOGRAMRAND_H_ */

// 
// wuROOTHistogramRand.hh ends here
