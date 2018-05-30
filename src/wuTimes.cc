// wuTimes.cc --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 六 2月 22 19:38:00 2014 (+0800)
// Last-Updated: 日 10月  9 21:54:30 2016 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 18
// URL: http://wuhongyi.cn 

#include"wuTimes.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

wuTimes::wuTimes()
{
  stime=0;
  etime=0;
  ProgressBarStart=0;
  ProgressBarStop=0;
  //std::cout<<"creating wuTimes..."<<std::endl;

  Percent=0;PercentOld=100;
}

wuTimes::~wuTimes()
{
  //std::cout<<std::endl<<"deleting wuTimes..."<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void wuTimes::CountTimeStart()
{
  time( &stime );
}

void wuTimes::CountTimeEnd()
{
  time( &etime );
}

void wuTimes::CoutTime()
{
  std::cout<<std::endl;
  std::cout<<"---------------TIME---------------------"<<std::endl;
  std::cout<<"The running time： "<<(etime-stime)/3600<<"hour "<<((etime-stime)%3600)/60<<"min "<<((etime-stime)%3600)%60<<"sec。"<<std::endl;
  std::cout<<"----------------------------------------"<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void wuTimes::ProgressBar(unsigned int i,unsigned int total)
{
  Percent=(int)(100.0*(double)i/(double)total);
  if(Percent!=PercentOld)
    {
      cout<<"\r"<<Percent<<"% [";

      int temp = (int)((double)Percent/(double)100 * 50);
      for (int i = 0; i < temp-1; ++i)
	{
	  cout << "=";
	}
      cout<<">";
      for (int j = 0; j <50 - temp-1; ++j)
	{
	  cout<<" ";
	}
      cout<<"]  ";
    }
  else
    {
      cout << "\r" ;
    }
  std::cout << flush;
  PercentOld=Percent;
}

void wuTimes::ProgressBarWithTime(unsigned int i,unsigned int total)
{
  if (i==0){ time(&ProgressBarStart);}

  Percent=(int)(100.0*(double)i/(double)total);
  if(Percent!=PercentOld)
    {

      cout<<"\r"<<Percent<<"% [";
      time(&ProgressBarStop);
      int temp = (int)((double)Percent/(double)100 * 50);
      for (int i = 0; i < temp-1; ++i)
	{
	  cout << "=";
	}
      cout<<">";
      for (int j = 0; j <50 - temp-1; ++j)
	{
	  cout<<" ";
	}
      cout<<"]   running： "<<(ProgressBarStop-ProgressBarStart)/3600<<"h "<<((ProgressBarStop-ProgressBarStart)%3600)/60<<"m "<<((ProgressBarStop-ProgressBarStart)%3600)%60<<"s。";
    }
  else
    {
      cout << "\r" ;
    }
  std::cout << flush;
  PercentOld=Percent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
string wuTimes::getTimeString()
{
  time_t lt = time(NULL);
  tm* current = localtime( &lt );
  char str[100];
  strftime( str , 100 , "%Y%m%d", current);

  return std::string(str);
}

string wuTimes::getTimeStringAll()
{
  time_t lt = time(NULL);
  tm* current = localtime( &lt );
  char str[100];
  strftime( str , 100 , "%Y%m%d%H%M%S", current);

  return std::string(str);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 
// wuTimes.cc ends here
