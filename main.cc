// main.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 六 4月 19 09:14:41 2014 (+0800)
// Last-Updated: 四 5月 31 05:39:30 2018 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 288
// URL: http://wuhongyi.cn 

#include "TRint.h"
#include "TObject.h"

#include <iostream>
#include "wuTimes.hh"
#include "wuSimulationPMT.hh"
#include "wuReadData.hh"

#include "TBenchmark.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TFormula.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TLeaf.h"
#include "TLeafC.h"
#include "TLeafD.h"
#include "TLeafElement.h"
#include "TLeafF.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TRandom.h"
#include "TSelector.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char *argv[])
{
  // Create an interactive ROOT application
  TRint *theApp = new TRint("Rint", &argc, argv);

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  //初始参数
  TString TreeName="t";//这里为要处理的文件中 tree name
  // create first the chain with all the files
  TChain *fChain=new TChain(TreeName);
  // std::cout << "Creating the chain" << std::endl;
  fChain->SetCacheSize(20*1024*1024);
  TString dir = gSystem->DirName(__FILE__);//获取当前文件main.cc所在路径 
  dir.ReplaceAll("/./","/");
  // std::cout<<dir<<std::endl;
  //=======================================================
  //以下两个选一个： 手动填写root文件所在路径 或 者直接使用当前文件所在路径
  // gSystem->Setenv("Dir","/home/wuhongyi");//手动填写路径
  string InputFile = wuReadData::ReadValue<string>("InputFileName","ReadData.txt");
  gSystem->Setenv("Dir",dir);//当前文件路径
  //=======================================================
  //将要处理的文件放在这里，支持tree名相同的多个结构相同的文件。特别适合用于Geant4多线程模拟的输出文件处理。
  string InputFileName = "$Dir"+InputFile;
  fChain->Add(InputFileName.c_str());
  // chain->Print();

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  // Declaration of leaf types
  Int_t           EventID;
  Int_t           TrackID;
  Char_t          PName[32];
  Char_t          CreatorProcess[32];
  Char_t          SDname[32];
  Double_t        GlobalTimePost;
  Double_t        EkPost;
  Double_t        xPost;
  Double_t        yPost;
  Double_t        zPost;
  Int_t           Information;
  // List of branches
  TBranch        *b_EventID;   //!
  TBranch        *b_TrackID;   //!
  TBranch        *b_PName;   //!
  TBranch        *b_CreatorProcess;   //!
  TBranch        *b_SDname;   //!
  TBranch        *b_GlobalTimePost;   //!
  TBranch        *b_EkPost;   //!
  TBranch        *b_xPost;   //!
  TBranch        *b_yPost;   //!
  TBranch        *b_zPost;   //!
  TBranch        *b_Information;   //!
  fChain->SetBranchAddress("EventID", &EventID, &b_EventID);
  fChain->SetBranchAddress("TrackID", &TrackID, &b_TrackID);
  fChain->SetBranchAddress("PName", PName, &b_PName);
  fChain->SetBranchAddress("CreatorProcess", CreatorProcess, &b_CreatorProcess);
  fChain->SetBranchAddress("SDname", SDname, &b_SDname);
  fChain->SetBranchAddress("GlobalTimePost", &GlobalTimePost, &b_GlobalTimePost);
  fChain->SetBranchAddress("EkPost", &EkPost, &b_EkPost);
  fChain->SetBranchAddress("xPost", &xPost, &b_xPost);
  fChain->SetBranchAddress("yPost", &yPost, &b_yPost);
  fChain->SetBranchAddress("zPost", &zPost, &b_zPost);
  fChain->SetBranchAddress("Information", &Information, &b_Information);

  std::cout <<std::endl<< "=== Running Hongyi Wu Analysis ===" << std::endl;

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  gBenchmark->Start("tree");//计时开始
  Long64_t TotalEntry = fChain->GetEntries();//拿到TChain中总entry行数
  wuTimes *wutimes = new wuTimes();//显示进度条

  Int_t ID = -1;
  string OutputFileName = wuReadData::ReadValue<string>("OutputFileName","ReadData.txt");
  TFile *newfile = new TFile(OutputFileName.c_str(),"RECREATE");
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  wuSimulationPMT leftPMT("left","XP2020");
  wuSimulationPMT rightPMT("right","XP2020");


  bool flag = false;

  TH1D *th1signal;

  TH1D *TDC = new TH1D("TDC","",2000,0,200);
  TH2D *TDCQDC = new TH2D("TDCQDC","",1000,1.0e10,1.0e12,1000,0,200);
  double tdc1,tdc2,qdc1,qdc2;
  TGraph *graphsignalleft;
  TGraph *graphsignalright;
  TGraph *singlephoton;
  singlephoton = rightPMT.GetGraphSinglePhoton("rightSinglePhoton");
  singlephoton ->Write();
  for (Long64_t entry = 0; entry < TotalEntry; ++entry)
    {//循环处理从这里开始

      fChain->GetEvent(entry);//这个是重点，拿到TChain中第entry行数据
      if(entry % 1000 == 0) wutimes->ProgressBarWithTime(entry ,TotalEntry);//必须从0开始才能触发计时起点

      if(EventID != ID)
	{
	  if(ID == 1)
	    {
	      // th1signal = leftPMT.GetTH1Signal("signal");
	      graphsignalleft = leftPMT.GetGraphSignal("leftgraph");
	      graphsignalright = rightPMT.GetGraphSignal("rightgraph");
	      newfile->cd();
	      graphsignalleft->Write();
	      graphsignalright->Write();
	    }

	  ID = EventID;
	  
	  if(flag)
	    {
	      leftPMT.GetTDCQDC(&tdc1,&qdc1);
	      rightPMT.GetTDCQDC(&tdc2,&qdc2);
	      TDCQDC->Fill((qdc1+qdc2)/2,(tdc1+tdc2)/2);
	      TDC->Fill((tdc1+tdc2)/2);
	    }
	  flag = false;

	  leftPMT.ResetDigitizer();
	  rightPMT.ResetDigitizer();
	}

      if(PName[0] != 'o')//跳过非optivalphoton
	{
	  continue;
	}

      flag = true;
      if(SDname[7] == '1')
	{
	  leftPMT.AddPhoton(GlobalTimePost,EkPost);
	}
      else
	{
	  rightPMT.AddPhoton(GlobalTimePost,EkPost);
	}

    }//循环处理到这里结束
  cout<<endl;


  // th1signal->Draw();
  // graphsignal2->Draw();
  // TDCQDC->Draw();
  newfile->cd();
  TDC->Write();
  TDCQDC->Write();

  newfile->Close();

  cout<<endl;
  gBenchmark->Show("tree");//计时结束并输出时间

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  // and enter the event loop...
  theApp->Run();

  delete theApp;

  return 0;

}




// 
// main.cc ends here
