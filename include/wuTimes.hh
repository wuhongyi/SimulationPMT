// wuTimes.hh --- 
// 
// Description: 
// Author: Hongyi Wu(吴鸿毅)
// Email: wuhongyi@qq.com 
// Created: 六 2月 22 19:37:41 2014 (+0800)
// Last-Updated: 三 4月 27 12:40:04 2016 (+0800)
//           By: Hongyi Wu(吴鸿毅)
//     Update #: 12
// URL: http://wuhongyi.cn 

#ifndef WUTIMES_H
#define WUTIMES_H

#include <ctime>
#include <iostream>
#include <string>
#include <cstring>
using namespace std;
class wuTimes
{
private:
  time_t stime;
  time_t etime;
  int Percent;
  int PercentOld;  
public:
  wuTimes();
  ~wuTimes();
  void CountTimeStart();//开始计时
  void CountTimeEnd();//结束计数
  void CoutTime();//输出运行时间

  void ProgressBar(unsigned int ,unsigned int );//进度条
  void ProgressBarWithTime(unsigned int ,unsigned int );//带时间的进度条,必须从0开始

  string getTimeString();//返回年月日
  string getTimeStringAll();//返回年月日时分秒

private:
  time_t ProgressBarStart;
  time_t ProgressBarStop;
};
#endif

// 
// wuTimes.hh ends here
