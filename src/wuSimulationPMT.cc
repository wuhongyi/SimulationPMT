// wuSimulationPMT.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 五 2月 19 21:24:11 2016 (+0800)
// Last-Updated: 二 3月 29 22:07:49 2016 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 145
// URL: http://wuhongyi.cn 

#include "wuSimulationPMT.hh"
#include "wuReadData.hh"

#include "TRandom.h"
#include "TMath.h"
#include <sstream>
#include <vector>
#include <iostream>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double wuSimulationPMTSinglePhotonelectronSignal1(Double_t *x, Double_t *par) 
{
  return x[0]*x[0]*TMath::Exp(-(x[0]*x[0])/(par[0]*par[0]));
}

double wuSimulationPMTSinglePhotonelectronSignal2(Double_t *x, Double_t *par) 
{
  if(x[0] > 0) 
    {
      return exp(-0.5*(pow((x[0]/par[0]),0.85)+exp(-x[0]/par[0])))/sqrt(2*TMath::Pi());
    }
  else
    {
      return exp(-0.5*(x[0]/par[0]+exp(-x[0]/par[0])))/sqrt(2*TMath::Pi());
    }		 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
wuSimulationPMT::wuSimulationPMT(string name,string PMTTytle,string fn)
{
  FileNameIn = fn;
  PMTName = PMTTytle;
  NameFlag = name;

  InitializationPMT();

  // cout<<FactorCollection_PMT.Value(3.0)<<endl;
  // cout<<QuantumEff_PMT.Value(3.0)<<endl;
}

wuSimulationPMT::~wuSimulationPMT()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void wuSimulationPMT::InitializationPMT()
{
  vector<double> ePhotonVec;
  int NL=wuReadData::ReadVector(PMTName+"CathodePhotonSpectrum",FileNameIn, &ePhotonVec);
  // cout<<endl;
  double *photonenergyrange=new double[NL];
  for(int i=0; i<NL; i++)
    {   
      photonenergyrange[i]=ePhotonVec[i];
      photonenergyrange[i]=6.62606957e-34*2.99792456e8/(photonenergyrange[i]*1.0e-9)/1.60217657e-19; //transfer from wavelength(nm) to energy (eV)
      // cout<<photonenergyrange[i]<<"  ";
    }

  int narr = 0;
  double *FactorCollection = new double[NL];  
  double *QuantumEff = new double[NL];

  narr = wuReadData::ReadArray(PMTName+"FactorCollection_PMT",FileNameIn,FactorCollection);
  if(narr<NL) for(int i=narr; i<NL; i++) FactorCollection[i] = FactorCollection[narr-1];

  narr = wuReadData::ReadArray(PMTName+"QuantumEff_PMT",FileNameIn,QuantumEff);
  if(narr<NL) for(int i=narr; i<NL; i++) QuantumEff[i] = QuantumEff[narr-1];

  for (int i = 0; i < NL; ++i)
    {
      FactorCollection_PMT.InsertValues(photonenergyrange[i],FactorCollection[i]);    
      QuantumEff_PMT.InsertValues(photonenergyrange[i],QuantumEff[i]);
    }
  delete [] FactorCollection;
  delete [] QuantumEff;
  delete [] photonenergyrange;

  FlagGain = wuReadData::ReadValue<int>(PMTName+"FlagGain",FileNameIn);
  if(FlagGain == 0 || FlagGain ==2 || FlagGain ==3) 
    {
      GainConst = wuReadData::ReadValue<double>(PMTName+"GainConst",FileNameIn);
      // cout<<FlagGain<<"  "<<GainConst<<endl;
    }
  if(FlagGain == 1)
    {
      vector<double> hv;
      vector<double> hvw;
      vector<double> hvg;
      int nhv=0;int nhvw=0;int nhvg=0;
      nhv=wuReadData::ReadVector(PMTName+"HighVoltage", FileNameIn, &hv);
      nhvw=wuReadData::ReadVector(PMTName+"HVWeight", FileNameIn, &hvw);
      RandHV.SetHistogram(&hv, &hvw, nhv);
      nhvg=wuReadData::ReadVector(PMTName+"HVGain", FileNameIn, &hvg);
      for (int i = 0; i < nhv; ++i)
	{
	  HVGain.InsertValues(hv[i],hvg[i]);
	}
    }
   if(FlagGain == 2)
     {
       GainSigma = wuReadData::ReadValue<double>(PMTName+"GainSigma",FileNameIn);
      // cout<<FlagGain<<"  "<<GainSigma<<endl;
     }
   if(FlagGain == 3)
     {
      GainM = wuReadData::ReadValue<double>(PMTName+"GainM",FileNameIn);
      
      stringstream ss;
      string temp;
      string tempm;
      string tempm1;
      string tempgainconst;
      ss.clear();
      ss<<GainM;
      ss>>tempm;
      ss.clear();
      ss<<(GainM-1);
      ss>>tempm1;
      ss.clear();
      ss<<GainConst;
      ss>>tempgainconst;
      temp="TMath::Power("+tempm+"*x/"+tempgainconst+","+tempm1+")*TMath::Exp(-"+tempm+"*x/"+tempgainconst+")*"+tempm+"/("+tempgainconst+"*TMath::Gamma("+tempm+"))";
      // cout<<temp<<endl;
      GainPolya = new TF1("GainPolya",temp.c_str(),0.5*GainConst,GainConst+0.5*GainConst);
      GainPolya -> SetNpx(10000);
     }


   RiseTime_PMT = wuReadData::ReadValue<double>(PMTName+"RiseTime_PMT",FileNameIn);
   
   InitializationDigitizer();//

}

void wuSimulationPMT::InitializationDigitizer()
{
  SinglePhotonTimeRange = wuReadData::ReadValue<double>("SinglePhotonTimeRange",FileNameIn);
  SignalTimeLeft = wuReadData::ReadValue<double>("SignalTimeLeft",FileNameIn);
  SignalTimeRight = wuReadData::ReadValue<double>("SignalTimeRight",FileNameIn);
  DigitizerTimeResolution = wuReadData::ReadValue<double>("DigitizerTimeResolution",FileNameIn);
  
  IntegralTimePreTriger = wuReadData::ReadValue<double>("IntegralTimePreTriger",FileNameIn);
    IntegralTimePostTriger = wuReadData::ReadValue<double>("IntegralTimePostTriger",FileNameIn);

  FlagOutputAlgorithm = wuReadData::ReadValue<int>("FlagOutputAlgorithm",FileNameIn);

  string tempname = NameFlag+"OutputSignal";
  OutputSignal = new TH1D(tempname.c_str(),"",int((SignalTimeRight-SignalTimeLeft)*1000.0/DigitizerTimeResolution),SignalTimeLeft,SignalTimeRight);

  if(FlagOutputAlgorithm == 1)
    {
      tempname = NameFlag+"SinglePhoton";
      SinglePhoton = new TH1D(tempname.c_str(),"",int(1000.0*SinglePhotonTimeRange/DigitizerTimeResolution),0,SinglePhotonTimeRange);
      // cout<<SinglePhotonTimeRange<<"  "<<DigitizerTimeResolution<<endl;

      tempname = NameFlag+"singlephoton";
      singlephoton  = new TF1(tempname.c_str(),wuSimulationPMTSinglePhotonelectronSignal1,0,SinglePhotonTimeRange,1);
      singlephoton->SetParameter(0,RiseTime_PMT);
      // cout<<singlephoton->Integral(0,3)<<"  "<<singlephoton->Integral(0,6)<<endl;

      FactorQ2V_PMT = wuReadData::ReadValue<double>("FactorQ2V_PMT",FileNameIn);
      CCouple_PMT = wuReadData::ReadValue<double>("CCouple_PMT",FileNameIn);
    }
  else
    {
      tempname = NameFlag+"SinglePhoton";
      SinglePhoton = new TH1D(tempname.c_str(),"",int(1000.0*1.25*SinglePhotonTimeRange/DigitizerTimeResolution),-0.25*SinglePhotonTimeRange,SinglePhotonTimeRange);

      TTS_FHWM = wuReadData::ReadValue<double>(PMTName+"TTS_FHWM",FileNameIn);
      SigmaTTS = 0.586*TTS_FHWM/(2*sqrt(2*TMath::Log(2)));
      tempname = NameFlag+"singlephoton";
      singlephoton  = new TF1(tempname.c_str(),wuSimulationPMTSinglePhotonelectronSignal2,-0.25*SinglePhotonTimeRange,SinglePhotonTimeRange,1);
      singlephoton->SetParameter(0,SigmaTTS);//参数TTS
    }


  //甄别方法
  FlagDiscriminantAlgorithm = wuReadData::ReadValue<int>("FlagDiscriminantAlgorithm",FileNameIn);
  if(FlagDiscriminantAlgorithm == 1)
    {
      LEDThreshold = wuReadData::ReadValue<double>("LEDThreshold",FileNameIn);
    }
  else
    {
      CFDThreshold = wuReadData::ReadValue<double>("CFDThreshold",FileNameIn);
    }

  GetSinglePhotonSignal();//

}

void wuSimulationPMT::GetSinglePhotonSignal()
{
  SinglePhoton->Reset("ICES");//每次使用前清空
  double tempt = 0;
  double temp;

  if(FlagOutputAlgorithm == 1)
    {
      for (int i = 1; i <= SinglePhoton->GetNbinsX(); ++i)
	{
	  tempt = SinglePhoton->GetBinCenter(i);
	  temp = tempt*tempt*TMath::Exp(-((tempt*tempt)/(RiseTime_PMT*RiseTime_PMT)));
	  SinglePhoton->SetBinContent(i, 1.6e-19*FactorQ2V_PMT*temp/(CCouple_PMT*singlephoton->Integral(0,SinglePhotonTimeRange)));//刨除增益项 GetGain()
	}
    }
  else
    {
      for (int i = 1; i <= SinglePhoton->GetNbinsX(); ++i)
	{
	  tempt = SinglePhoton->GetBinCenter(i);
	  temp = wuSimulationPMTSinglePhotonelectronSignal2(&tempt,&SigmaTTS);//参数TTS
	  SinglePhoton->SetBinContent(i, temp);//刨除增益项 GetGain()
	}
    }

  // TFile *f = new TFile("hhh.root","RECREATE");
  // SinglePhoton->Write();
  // f->Close();
  // SinglePhoton->Draw();
}

double wuSimulationPMT::GetGain()
{
  if(FlagGain == 0)
    {
      return GainConst;
    }
  else
    {
      if(FlagGain == 1)
	{
	  return HVGain.Value(RandHV.GetRand());
	}
      else
	{
	  if(FlagGain == 2)
	    {
	      return gRandom->Gaus(GainConst,GainSigma);
	    }
	  else
	    {
	      if(FlagGain == 3)
		{
		  return  GainPolya->GetRandom();
		}
	      else
		{
		  cout<<"FlagGain 不在指定范围！"<<endl;
		  return 0;
		}
	    }
	}
    }
}

void wuSimulationPMT::AddPhoton(double time,double energy)
{
  PhotonArrivedEnergy = 1.0e6*energy;
  PhotonArrivedTime = time;

  //每添加一个光子就画在总图上
  if(gRandom->Rndm() > FactorCollection_PMT.Value(PhotonArrivedEnergy)) return;
  if(gRandom->Rndm() > QuantumEff_PMT.Value(PhotonArrivedEnergy)) return;

  if(FlagOutputAlgorithm == 1)
    {
      for (int i = 1; i <= SinglePhoton->GetNbinsX(); ++i)
  	{
  	  OutputSignal->Fill(PhotonArrivedTime+(i-1)*DigitizerTimeResolution/1000.0,GetGain()*SinglePhoton->GetBinContent(i));//每个光子的输出增益在这里体现，还没考虑渡越时间的展宽！！！
  	}
    }
  else
    {
      for (int i = 1; i <= SinglePhoton->GetNbinsX(); ++i)
  	{
	  //SinglePhoton 从-0.25Range开始而不是从零开始
  	  OutputSignal->Fill(PhotonArrivedTime-0.25*SinglePhotonTimeRange+(i-1)*DigitizerTimeResolution/1000.0,GetGain()*SinglePhoton->GetBinContent(i));
  	}
    }

}

void wuSimulationPMT::ResetDigitizer()
{
  OutputSignal->Reset("ICES");//每次使用前清空
}


TH1D* wuSimulationPMT::GetTH1Signal(const char* name)
{
  return (TH1D *)OutputSignal->Clone(name);
}

TGraph* wuSimulationPMT::GetGraphSignal(const char* name)
{
  TGraph *graph = new TGraph();
  graph->SetName(name);
    for (int i = 1; i <= OutputSignal->GetNbinsX(); ++i)
	{
	  graph->SetPoint(i-1,OutputSignal->GetBinCenter(i) ,OutputSignal->GetBinContent(i));
	}
  return graph;
}

TH1D* wuSimulationPMT::GetTH1SinglePhoton(const char* name)
{
  return (TH1D *)SinglePhoton->Clone(name);
}

TGraph* wuSimulationPMT::GetGraphSinglePhoton(const char* name)
{
  TGraph *graph = new TGraph();
  graph->SetName(name);
    for (int i = 1; i <= SinglePhoton->GetNbinsX(); ++i)
	{
	  graph->SetPoint(i-1,SinglePhoton->GetBinCenter(i) ,SinglePhoton->GetBinContent(i));
	}
  return graph;
}

double wuSimulationPMT::GetTDC()
{
  int triger;
  if(FlagDiscriminantAlgorithm == 1)//前沿甄别
    {
      triger = OutputSignal->FindFirstBinAbove(LEDThreshold);
    }
  else
    {
      triger = OutputSignal->FindFirstBinAbove(CFDThreshold*OutputSignal->GetBinContent(OutputSignal->GetMaximumBin()));
    }
  return OutputSignal->GetBinCenter(triger);

}

double wuSimulationPMT::GetQDC()
{
  int triger;
  if(FlagDiscriminantAlgorithm == 1)
    {
      triger = OutputSignal->FindFirstBinAbove(LEDThreshold);
    }
  else
    {
      triger = OutputSignal->FindFirstBinAbove(CFDThreshold*OutputSignal->GetBinContent(OutputSignal->GetMaximumBin()));
    }

  if(triger-IntegralTimePreTriger*1000/DigitizerTimeResolution <1) 
    {
      // cout<<"参数 IntegralTimePreTriger 设置过大。"<<endl;
      return -1;
    }
  if(triger+IntegralTimePostTriger*1000/DigitizerTimeResolution > OutputSignal->GetNbinsX()) 
    {
      // cout<<"参数 IntegralTimePostTriger 设置过大。"<<endl;
      return -1;
    }
  return OutputSignal->Integral(triger-IntegralTimePreTriger*1000/DigitizerTimeResolution,triger+IntegralTimePostTriger*1000/DigitizerTimeResolution);

}

void wuSimulationPMT::GetTDCQDC(double *tdc,double *qdc)
{
  int triger;
  if(FlagDiscriminantAlgorithm == 1)
    {
      triger = OutputSignal->FindFirstBinAbove(LEDThreshold);
    }
  else
    {
      triger = OutputSignal->FindFirstBinAbove(CFDThreshold*OutputSignal->GetBinContent(OutputSignal->GetMaximumBin()));
    }

      *tdc = OutputSignal->GetBinCenter(triger);

      if(triger-IntegralTimePreTriger*1000/DigitizerTimeResolution <1) 
	{
	  // cout<<"参数 IntegralTimePreTriger 设置过大。"<<endl;
	  *qdc = -1;return;
	}
      if(triger+IntegralTimePostTriger*1000/DigitizerTimeResolution > OutputSignal->GetNbinsX()) 
	{
	  // cout<<"参数 IntegralTimePostTriger 设置过大。"<<endl;
	  *qdc = -1;return;
	}
      *qdc = OutputSignal->Integral(triger-IntegralTimePreTriger*1000/DigitizerTimeResolution,triger+IntegralTimePostTriger*1000/DigitizerTimeResolution);
      
}


// 
// wuSimulationPMT.cc ends here
