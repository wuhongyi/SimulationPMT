#CMakeLists.txt --- 
# 
# Description: CMake编译、运行程序模板。
# Author: Hongyi Wu(吴鸿毅)
# Email: wuhongyi@qq.com 
# Created: 六 2月 22 19:39:08 2014 (+0800)
# Last-Updated: 四 5月 31 05:35:33 2018 (+0800)
#           By: Hongyi Wu(吴鸿毅)
#     Update #: 112
# URL: http://wuhongyi.cn 

#cmake编译、运行程序模板！
#使用方法：将所有点hh文件放置到include文件夹，点cc文件放置到src，动态、静态链接库放置到lib
#将主函数放在与本CMakeLists.txt并列的目录中。
#在当前目录打开终端，输入        cmake .[Enter]   make
#或在build目录下打开终端，输入cmake ..[Enter]   make
#----------------------------------------------------------------------------
#判断cmake版本是否符合，对工程命名
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CODE_DIR /home/wuhongyi/CodeProject )   #设置工程文件路径
#set(CMAKE_C_COMPILER "/usr/local/bin/gcc")  #调用自己安装在local的gcc
#set(CMAKE_CXX_COMPILER "/usr/local/bin/g++") #调用自己安装在local的g++
set(PROJ_NAME wu)
project(${PROJ_NAME})
set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated") #支持最新c++11
#c++0x是新的C++意思。而-Wno-deprecated是表示 no-deprecated这类警告忽略。
SET(CMAKE_BUILD_TYPE Debug) #编译调试版本
#----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CODE_DIR}/cmake) #设置 FindROOT.cmake 文件路径
# Find ROOT (required package)
find_package(ROOT REQUIRED)
# find_package( ROOT COMPONENTS ASImage EG FFTW FitPanel Ged Geom Gui GuiBld GuiHtml GX11 GX11TTF Hbook Hist Html MathMore Minuit Minuit2 MLP Proof RGL RooFit RooFitCore Smatrix Spectrum Table TMVA TreePlayer TreeViewer )#COMPONENTS可选 ，REQUIRED COMPONENTS必须有
#----------------------------------------------------------------------------
#指出所需hh，cc，so，a文件位置路径，并用变量命名
include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/* )
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*)
file(GLOB libs ${PROJECT_SOURCE_DIR}/lib/* )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build) #生成可执行文件放在bin文件夹中


#----------------------------------------------------------------------------
#add_executable()后面的是一个参数列表，可带多个要编译的文件名，中间以空格或回车等隔开。生成一个可执行文件
#add_executable(${PROJ_NAME} main.cc wuDict.cxx ${sources} ${headers}) #main.cc
add_executable(${PROJ_NAME} main.cc ${sources} ${headers}) #main.cc
#链接库，后面补上名称
target_link_libraries(${PROJ_NAME} ${libs} ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )

#----------------------------------------------------------------------------
#将所需数据文件复制到bin文件夹
set(CopyFile
    ReadData.txt rootlogon.C rootlogoff.C PMT.txt
  )

foreach(copyfile ${CopyFile})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${copyfile}
    ${PROJECT_SOURCE_DIR}/build/${copyfile}
    COPYONLY
    )
endforeach()

###################################
