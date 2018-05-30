// wuSimulationPMT.hh --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 五 2月 19 21:23:43 2016 (+0800)
// Last-Updated: 二 3月 29 19:56:27 2016 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 53
// URL: http://wuhongyi.cn 

// 模型参考文献
// Nuclear Instruments and Methods in Physics Research A 491(2002)54-68
// CPC(HEP&NP),2009,33(10):860-865
// 北京大学博士后研究工作报告，2009，宋玉收
// Nuclear Instruments and Methods in Physics Research A 624(2010)583-590

#ifndef _WUSIMULATIONPMT_H_
#define _WUSIMULATIONPMT_H_

#include "wuPhysicsOrderedFreeVector.hh"
#include "wuROOTHistogramRand.hh"

#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TObject.h"
#include "TGraph.h"


#include <cstring>
#include <string>
#include <cmath>

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class wuSimulationPMT
{
public:
  wuSimulationPMT(string name,string PMTTytle = "XP2020",string fn = "PMT.txt");//必须保证每个类的name不一样（name内存中用到）
  virtual ~wuSimulationPMT();

  void ResetDigitizer();//每个Event前都得调用
  void AddPhoton(double time,double energy);//添加光子信息，一个Event中多少个光子到达PMT就得调用多少次

  double GetTDC();//获得时间
  double GetQDC();//获得能量
  void GetTDCQDC(double *tdc,double *qdc);

  TH1D* GetTH1Signal(const char* name);//获得该Event的输出
  TGraph* GetGraphSignal(const char* name);//获得该Event的输出

  TH1D* GetTH1SinglePhoton(const char* name);//获得单光子的等效输出（不含增益）
  TGraph* GetGraphSinglePhoton(const char* name);//获得单光子的等效输出（不含增益）

private:
  void InitializationPMT();
  void InitializationDigitizer();

  void ResetPMT();

  double GetGain();
  void GetSinglePhotonSignal();//不包含Gain参数，避免每次重复计算

private:


private:
  double PhotonArrivedEnergy;
  double PhotonArrivedTime;

  TF1 *singlephoton;
  TH1D *SinglePhoton;
  TH1D *OutputSignal;

  double SinglePhotonTimeRange;//ns
  double SignalTimeLeft;//ns
  double SignalTimeRight;//ns
  double DigitizerTimeResolution;//ps

  double IntegralTimePreTriger;//ns
  double IntegralTimePostTriger;//ns

private:
  string FileNameIn;
  string PMTName;
  string NameFlag;

  wuPhysicsOrderedFreeVector FactorCollection_PMT;//能量与收集因子
  wuPhysicsOrderedFreeVector QuantumEff_PMT;//能量与量子效率


  int FlagGain;//0 Const; 1 HV~G; 2 Gaussian; 3 Polya distribution.
  double GainConst;
  double GainSigma;
  wuROOTHistogramRand RandHV;//高压谱  GetRand()
  wuPhysicsOrderedFreeVector HVGain;//高压与增益
  double GainM;//玻利亚分布
  TF1 *GainPolya;
  
  double RiseTime_PMT;//上升时间


  int FlagOutputAlgorithm;//单光子输出电流算法 1-2002;2-2010
  //2002参数
  double FactorQ2V_PMT;
  double CCouple_PMT;
  //2010参数
  double TTS_FHWM;//ns
  double SigmaTTS;
  
  int FlagDiscriminantAlgorithm;//1-前沿甄别；2-恒比定时甄别
  double LEDThreshold;//前沿甄别
  double CFDThreshold;//恒比定时甄别

};

#endif /* _WUSIMULATIONPMT_H_ */
// 
// wuSimulationPMT.hh ends here

