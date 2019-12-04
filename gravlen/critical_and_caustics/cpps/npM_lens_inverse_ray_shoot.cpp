#include <iostream>
#include "argparse.h"
#include <math.h>
#include <cstdlib>
#include <ctime>
#include<cmath>
#include <unistd.h>
// #include "ppl.h"//Concurrency::parallel_for
// https://raw.githubusercontent.com/icestudent/vc-19-changes/master/ppl.h

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>

// #include<string>
// #include <array>

#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

//https://www.cnblogs.com/lonelypinky/p/11579911.html
//ubuntu下安装Opencv 4.1.1


//g++ npM_lens_inverse_ray_shoot.cpp `pkg-config opencv4 --libs --cflags` -o mainn
//./mainn --masN 10 --imgS 1024 --rayN 3000 --rayL 3.5 --pmrl 0 --pmru 1 --pxyr 1.5


//error: ISO C++ forbids declaration of ‘inverse_ray_shooting’ with no type [-fpermissive]
//https://stackoverflow.com/questions/7929477/explain-the-error-iso-c-forbids-declaration-of-personlist-with-no-type

//https://segmentfault.com/a/1190000015653101
//C++下的OpenCV中Mat类型存储的图像格式

// https://docs.opencv.org/2.4/doc/tutorials/highgui/video-write/video-write.html
// creat video using opencv

/*

https://bbs.csdn.net/topics/390901939?list=lz
用openmp解决了。


https://blog.csdn.net/lianshaohua/article/details/82390026
C++高性能的并行数据结构

https://justyy.com/archives/2257
C++ 使用 Parallel For 多线程 计算 圆周率 – Monte Carlo

https://www.cnblogs.com/allenlooplee/archive/2012/07/25/2608820.html
遇见C++ PPL：C++ 的并行和异步

https://blog.csdn.net/qq_28314759/article/details/52604726
简单了解并行计算库PPL

https://stackoverflow.com/questions/11564749/using-parallel-patterns-library-ppl-h
Using Parallel Patterns Library (ppl.h)

https://docs.microsoft.com/en-us/cpp/parallel/concrt/parallel-patterns-library-ppl?redirectedfrom=MSDN&view=vs-2019
Parallel Patterns Library (PPL), Microsoft


https://stackoverflow.com/questions/36246300/parallel-loops-in-c/36246386
Parallel Loops in C++

Guide into OpenMP: Easy multithreading programming for C++
https://bisqwit.iki.fi/story/howto/openmp/
*/



const float PI = std::atan(1.0)*4;
const float arcsec2rad = PI/180/3600;
const float rad2arcsec = 180*3600/PI;

// 生成有n个元素的随机数组,每个元素的随机范围为[rangeL, rangeR] 闭区间
void generateRandomArray(float *arr, int n, float rangeL, float rangeR, int srandnum) {
    assert(rangeL <= rangeR);
    
    // float arr = new float[n]; // 创建一个 n个元素的数组
    // float arr[n];

    srand(srandnum); // 随机种子, time(NULL)
    for (int i = 0; i < n; i++)
        arr[i] = (rand()/double(RAND_MAX)) * (rangeR - rangeL) + rangeL;
    // rand() / double(RAND_MAX)
    //要取得0～1之间的浮点数，可以使用rand() / double(RAND_MAX)。
    // return arr;
}

class Onelens
{
public:
    float mass, *pos, dis;
    Onelens(float mass, float *pos, float dis): mass(mass), pos(pos), dis(dis){}    
};


class Nlenses
{
public:
    float *masses, *xs, *ys;
    int pointsNum;
    // float massratio, mu1, mu2, *betax, *betay, *mag, m1, m2, beta;
public:
    Nlenses(float *masses,float *xs,float *ys,int pointsNum);
    void inverse_ray_shooting(float thetax, float thetay, float *betaxy);
    float comp_mag_samez(float thetax, float thetay);
    float ray_shoot_comp_mag_in_one(float thetax, float thetay, float *betaxy);
    void get_imgs_lessmem(Mat srcplaneIMG, Mat srcplaneIMG_withoutlens, Mat imgplaneIMG, unsigned short int *ImgSize, float *xlim, float *ylim, unsigned short int raynum_1side);
};

Nlenses::Nlenses(float *masses,float *xs,float *ys,int pointsNum)
{
float msum = 0;
// int len = sizeof( masses ) / sizeof( *masses );

// scale total mass to 1

for(int i=0;i<pointsNum;i++){
    msum+=masses[i];}
for(int i=0;i<pointsNum;i++){
    masses[i]/=msum;}

this -> masses = masses;
this -> xs = xs;
this -> ys = ys;
this -> pointsNum = pointsNum;
}

void Nlenses::inverse_ray_shooting(float thetax, float thetay, float *betaxy)
{
    // float r1_2 = pow(thetax - this->lens1->pos[0],2) + pow(thetay - this->lens1->pos[1],2);
    // float r2_2 = pow(thetax - this->lens2->pos[0],2) + pow(thetay - this->lens2->pos[1],2);

    // betaxy[0] = thetax - this->mu1*(thetax - this->lens1->pos[0])/r1_2 - this->mu2*(thetax - this->lens2->pos[0])/r2_2;
    // betaxy[1] = thetay - this->mu1*(thetay - this->lens1->pos[1])/r1_2 - this->mu2*(thetay - this->lens2->pos[1])/r2_2;
    
    betaxy[0] = thetax;
    betaxy[1] = thetay;
    float Ri2;
    for (int i=0;i<pointsNum;i++){
    Ri2 = pow(thetax - this->xs[i],2) + pow(thetay - this->ys[i],2);
    betaxy[0] -= this->masses[i]*(thetax - this->xs[i])/Ri2;
    betaxy[1] -= this->masses[i]*(thetay - this->ys[i])/Ri2;
    }
}

float Nlenses::comp_mag_samez(float thetax, float thetay)
{
    float mag,A11=1,A12=0,A22=1,detA,Ri2;


    // R1 = sqrt(pow(thetax-x1,2)+pow(thetay-y1,2));
    // R2 = sqrt(pow(thetax-x2,2)+pow(thetay-y2,2));
for (int i=0;i<pointsNum;i++){
    // A11 = 1 - this->mu1*(1/pow(R1, 2) - 2*pow(thetax-x1,2)/pow(R1,4)) - this->mu2*(1/pow(R2, 2) - 2*pow(thetax-x2,2)/pow(R2,4));
    // A12 = -2*this->mu1*(thetax-x1)*(thetay-y1)/pow(R1,4) - 2*this->mu2*(thetax-x2)*(thetay-y2)/pow(R2,4);
    // A22 = 1 - this->mu1*(1/pow(R1, 2) - 2*pow(thetay-y1,2)/pow(R1,4)) - this->mu2*(1/pow(R2, 2) - 2*pow(thetay-y2,2)/pow(R2,4));
    Ri2 = pow(thetax-this->xs[i],2)+pow(thetay-this->ys[i],2);
    A11 -= this->masses[i]*(1/Ri2 - 2*pow(thetax-this->xs[i],2)/pow(Ri2,2));
    A22 -= this->masses[i]*(1/Ri2 - 2*pow(thetay-this->ys[i],2)/pow(Ri2,2));
    A12 -= 2*this->masses[i]*(thetax-this->xs[i])*(thetay-this->ys[i])/pow(Ri2,2);

}
    detA = (A11*A22 - pow(A12,2));
    mag = 1/abs(detA);
    return mag;
}

float Nlenses::ray_shoot_comp_mag_in_one(float thetax, float thetay, float *betaxy)
{
    float mag,A11=1,A12=0,A22=1,detA,Ri2;
    betaxy[0] = thetax;
    betaxy[1] = thetay;

for (int i=0;i<pointsNum;i++){
    Ri2 = pow(thetax-this->xs[i],2)+pow(thetay-this->ys[i],2);

    betaxy[0] -= this->masses[i]*(thetax - this->xs[i])/Ri2;
    betaxy[1] -= this->masses[i]*(thetay - this->ys[i])/Ri2;
    
    A11 -= this->masses[i]*(1/Ri2 - 2*pow(thetax-this->xs[i],2)/pow(Ri2,2));
    A22 -= this->masses[i]*(1/Ri2 - 2*pow(thetay-this->ys[i],2)/pow(Ri2,2));
    A12 -= 2*this->masses[i]*(thetax-this->xs[i])*(thetay-this->ys[i])/pow(Ri2,2);

}
    detA = (A11*A22 - pow(A12,2));
    mag = 1/abs(detA);
    return mag;
}



void Nlenses::get_imgs_lessmem(Mat srcplaneIMG, Mat srcplaneIMG_withoutlens, Mat imgplaneIMG, unsigned short int *ImgSize, float *xlim, float *ylim, unsigned short int raynum_1side)
{
    float ratiox = (ImgSize[1]-1)/(xlim[1]-xlim[0]);
    float ratioy = (ImgSize[0]-1)/(ylim[1]-ylim[0]);
    float incx = (xlim[1]-xlim[0])/raynum_1side; // increment of x 
    float incy = (ylim[1]-ylim[0])/raynum_1side; // increment of y

    unsigned long int cnt = 1;
    float betaxy[2] = {0, 0};
    float thetay, thetax;
    // https://justyy.com/archives/2257

    // float thetax = xlim[0];
    // while(thetax<xlim[1]){

    // for (float thetax = xlim[0]; thetax<xlim[1]; thetax=thetax+incx){
    for (int ix = 0; ix<raynum_1side; ix++){
    // Concurrency::parallel_for(int(0); raynum_1side; [&](int i){

    thetax = xlim[0]+ix*incx;

    thetay = ylim[0];
    while(thetay<ylim[1]){
    // for (float thetax = xlim[0]; thetax<=xlim[1]; thetax=thetax+incx){
    // for (float thetay = ylim[0]; thetax<=ylim[1]; thetay=thetay+incy){

    std::cout.width(3);//i的输出为3位宽
    std::cout << cnt*100/(raynum_1side*raynum_1side) << "%";
    std::cout << "\b\b\b\b";//回删三个字符，使数字在原地变化

    // std::cout << cnt<<", " <<thetax <<", " << thetay << "\n";
    //运行时，对内存操作有误，常报Segmentation fault，Core Dump
    cnt = cnt + 1;

    //http://www.cplusplus.com/forum/beginner/106352/
    //c++ class calling one member function from another
    

    // float betaxy[2] = {0, 0};
    // betax, betay = this -> inverse_ray_shooting(thetax, thetay, betaxy)
    // *betaxy = this -> inverse_ray_shooting(thetax, thetay, betaxy);

  
    // this -> inverse_ray_shooting(thetax, thetay, betaxy);
    // float mag = this -> comp_mag_samez(thetax, thetay);
    
    // 
    float mag = this -> ray_shoot_comp_mag_in_one(thetax, thetay, betaxy);

    float xdata, ydata, betaxdata, betaydata;
    xdata = thetax - xlim[0];// thetax
    ydata = thetay - ylim[0];// thetax
    betaxdata = betaxy[0] - xlim[0];
    betaydata = betaxy[1] - ylim[0];

    // mat_CV_8UC1.at<uchar>(nrow,ncol);   
    // https://blog.csdn.net/csdn_zhishui/article/details/83039207
    // OpenCV::Mat 修改特定位置像素的值

    /**
    核心：
    *   *(img.data + img.step[0] * i + img.step[1] * j + img.elemSize1() * c)=new_val
    *   这行代码是解析出来image的第i行第j列(即坐标为[j,i])第c通道，然后就可以对它进行赋值了。
    */

    long int bi, bj, i, j; //65441,
    bj = round(betaydata*ratioy);
    bi = round(betaxdata*ratiox);
    j = round(ydata*ratioy);
    i = round(xdata*ratiox);
    
    // std::cout <<  bi<<", " <<bj <<", " << i <<", " << j << "\n\n";

    // try{
    //     // srcplaneIMG[round(betaydata*ratioy), round(betaxdata*ratiox)]+=mag
    //     // srcplaneIMG[bj, bi]+=mag
    //     *(srcplaneIMG.data + srcplaneIMG.step[0] * bi + srcplaneIMG.step[1] * bj + srcplaneIMG.elemSize1() * 0) += mag;
    // }
    // catch (...) {
    //     cout << "catch (...)" << endl;
    // }

    // cout << "\nmag: " << mag;

    if ((bi>=0 && bi<ImgSize[1])&&(bj>=0&&bj<ImgSize[0])){
    // *(srcplaneIMG.data + srcplaneIMG.step[0] * bi + srcplaneIMG.step[1] * bj + srcplaneIMG.elemSize1() * 0) += mag;
    srcplaneIMG.at<float>(bj,bi) += mag;

}
      
    // srcplaneIMG_withoutlens[j, i]+=1
    if ((i>=0 && i<ImgSize[1])&&(j>=0&&j<ImgSize[0])){
    // *(srcplaneIMG_withoutlens.data + srcplaneIMG_withoutlens.step[0] * i + srcplaneIMG_withoutlens.step[1] * j + srcplaneIMG_withoutlens.elemSize1() * 0) += 1;
    srcplaneIMG_withoutlens.at<float>(j,i) += 1.f;
    // imgplaneIMG[j, i]+=mag    
    // *(imgplaneIMG.data + imgplaneIMG.step[0] * i + imgplaneIMG.step[1] * j + imgplaneIMG.elemSize1() * 0) += mag;
    imgplaneIMG.at<float>(j,i) += mag;
}

    thetay += incy;
    }
    // thetax += incx;
}

// srcplaneIMG /= srcplaneIMG_withoutlens;
// Mat result;
// result = srcplaneIMG.clone();
    // int rowNumber = result.rows;
    // int colNumber = result.cols;

    // computing real magnification in source plane to get the Caustics

    // for (int i = 0; i < ImgSize[0]; i++)
    // {
    //     for (int j = 0; j < ImgSize[1]; j++){
    //     srcplaneIMG.at<float>(j,i) = srcplaneIMG.at<float>(j,i) / srcplaneIMG_withoutlens.at<float>(j,i);
    //     }
    // }


// return void;
}

//https://github.com/jamolnng/argparse
//A simple C++ header only command line argument parser
//


int main(int argc, char* argv[])
{

// float mass1, posX;
// unsigned short int raynum_1side;
// unsigned short int imgsize;
// unsigned short int pointsNum;

unsigned short int raynum_1side, imgsize, pointsNum;
float raylim, masses_rangel, masses_rangeh, pmxlim, pmylim;

// float raylim = 2.5;

// raynum_1side = 1000;
// imgsize = 1024;

ArgumentParser parser("Argument parser");
parser.add_argument("--masN", "an int, number of point masses, e.g. 10");
parser.add_argument("--imgS", "an int, image size, e.g. 1024");
parser.add_argument("--rayN", "an int, ray number in one side, e.g. 3000");
parser.add_argument("--rayL", "an positive float, light ray range in lens plane, e.g. 2.5");
parser.add_argument("--pmrl", "an float, low mass limit of point masses, e.g. 0");
parser.add_argument("--pmru", "an float, up mass limit of point masses, e.g. 1");
parser.add_argument("--pxyr", "an positive float, limit of x,y position of point masses, e.g. 1");
// parser.add_argument("--pyr", "an positive float, limit of y position of point masses, e.g. 1");


// parser.add_argument("-p", "an float number, position X of one of the point masses");
  try {
    parser.parse(argc, argv);
  } catch (const ArgumentParser::ArgumentNotFound& ex) {
    std::cout << ex.what() << std::endl;
    return 0;
  }


if (parser.exists("h")){
cout <<"N point mass lenses, ploting the Critical lines and Caustics"<<endl;
// cout<<"\b \b";
// cin.get();
return 0;
}

pointsNum = parser.get<unsigned short int>("masN");
raynum_1side = parser.get<unsigned short int>("rayN");
imgsize = parser.get<unsigned short int>("imgS");
raylim = parser.get<float>("rayL");
masses_rangel= parser.get<float>("pmrl");
masses_rangeh= parser.get<float>("pmru");
pmxlim = parser.get<float>("pxyr");
// pmylim = parser.get<float>("pxyr");
pmylim = pmxlim;

// cout<<pointsNum<<","<<pointsNum<<","<<
// raynum_1side<<","<<
// imgsize<<","<<
// raylim<<","<<
// masses_rangel<<","<<
// masses_rangeh<<","<<
// pmxlim<<","<<
// pmylim<<endl;

// pointsNum = 2;
// raynum_1side = 1000;
// imgsize = 1024;
// raylim = 2.5;
// masses_rangel = 0;
// masses_rangeh = 1;
// pmxlim = 1;
// pmylim = 1;

float masses_range[2] = {masses_rangel,masses_rangeh};
float xs_range[2] = {-pmxlim,pmxlim};
float ys_range[2] = {-pmylim,pmylim};

// mass1 = parser.get<float>("m");
// posX = parser.get<float>("p");
// cout<< mass1 << "\n"<<raynum_1side << "\n"<<imgsize << "\n";


// using namespace PersonNamespace;
// float pos1[2] = {-posX, 0};
// float mass1 = 1;


// cout<< argv[1]<< ","<<*argv[1]<< ","<< argc<< "\n";
// char *temp = argv[1];
// stringstream ss;
// ss<<argv[1];
// ss>>mass1;
// mass1 = (float) (*temp);
// cout<< mass1 << "\n";



// Onelens lens1(mass1, pos1 , 0);
// float pos2[2] = {posX, 0};
// Onelens lens2(1.0, pos2 , 0);

// float massratio = lens1.mass/lens2.mass;

// if(abs(posX)>=1.5){lim = abs(posX)+2.5;}
// if((massratio)<-0.5){lim = abs(posX)+2.5;}
float xlim[2] = {-raylim,raylim};
float ylim[2] = {-raylim,raylim};
// float xlim[2] = {-2.5,2.5};
// float ylim[2] = {-2.5,2.5};

unsigned short int imshowsize = 640;
unsigned short int ImgSize[2] = {imgsize, imgsize}; // raw, colume
//unsigned short int    2 个字节   0 到 65,535



// THETAX, THETAY = genxy(xlim=xlim,ylim=ylim,num=num)
// def genxy(xlim=(-1,1),ylim=(-1,1),num=100):
//     x = np.linspace(xlim[0],xlim[1],num)
//     y = np.linspace(ylim[0],ylim[1],num)
//     X,Y = np.meshgrid(x,y)
//     # return X.reshape(1,-1), Y.reshape(1,-1) # shape: (1,10)
//     return X.reshape(1,-1)[0], Y.reshape(1,-1)[0] # shape: (3,)


// Nlenses twolens( &lens1,  &lens2);


float masses[pointsNum];
float xs[pointsNum];
float ys[pointsNum];

// int srandnum = time(NULL);
generateRandomArray(masses, pointsNum, masses_range[0],masses_range[1],time(NULL));
sleep(1); //延时1秒 
generateRandomArray(xs, pointsNum, xs_range[0],xs_range[1],time(NULL));
sleep(1); //延时1秒 
generateRandomArray(ys, pointsNum, ys_range[0],ys_range[1],time(NULL));


// float masses[pointsNum]= {1,1};
// float xs[pointsNum]= {-0.3, 0.3};
// float ys[pointsNum]={0,0};


// unsigned short int pointsNum2 = 9;
// float scalexy = 0.8;
// float scalem = -0.2;//-0.2;
// float masses2[pointsNum2]= {1,1*scalem,1,1*scalem,1,1*scalem,1,1*scalem,1};
// float xs2[pointsNum2]= {-1,0,1,-1,0,1,-1,0,1};
// float ys2[pointsNum2]={1,1,1,0,0,0,-1,-1,-1};
// for( int i=0;i<pointsNum2;i++ ){
//     xs2[i]*=scalexy;
//     ys2[i]*=scalexy;
// }

// Nlenses nlens(masses2, xs2, ys2, pointsNum2);


Nlenses nlens(masses, xs, ys, pointsNum);

cout <<"masses: ";
for( int i=0;i<pointsNum;i++ )
    cout << masses[i] <<" " ;
cout <<endl ;
cout <<"xs: ";
for( int i=0;i<pointsNum;i++ )
    cout  << xs[i] <<" " ;
cout <<endl ;
cout <<"ys: ";
for( int i=0;i<pointsNum;i++ )
    cout  << ys[i] <<" " ;
cout <<endl ;

//Mat Matrix_name(行数，列数，存储元素的数据类型，每个矩阵点的通道数)
// Mat srcplaneIMG(ImgSize[0], ImgSize[1], CV_8UC1, 1);
// Mat srcplaneIMG_withoutlens(ImgSize[0], ImgSize[1], CV_8UC1, 1);
// Mat imgplaneIMG(ImgSize[0], ImgSize[1], CV_8UC1, 1);

Mat srcplaneIMG(ImgSize[0], ImgSize[1], CV_32F, 1);
Mat srcplaneIMG_withoutlens(ImgSize[0], ImgSize[1], CV_32F, 1);
Mat imgplaneIMG(ImgSize[0], ImgSize[1], CV_32F, 1);

nlens.get_imgs_lessmem(srcplaneIMG, srcplaneIMG_withoutlens, imgplaneIMG, ImgSize, xlim, ylim, raynum_1side);
// ./mainn --masN 10 --imgS 2048 --rayN 2000 --rayL 2 --pmrl 0 --pmru 1 --pxyr   12.53s user 0.23s system 83% cpu 15.260 total

//cv::Vec3b vec3b      = img.at<cv::Vec3b>(0,0);
// uchar     vec3b0     = img.at<cv::Vec3b>(0,0)[0];
//https://blog.csdn.net/xiaowei_cqu/article/details/19839019 OpenCV】访问Mat中每个像素的值（新）


Mat srcplaneIMG_8U, imgplaneIMG_8U;
srcplaneIMG.convertTo( srcplaneIMG_8U, CV_8U);
imgplaneIMG.convertTo( imgplaneIMG_8U, CV_8U);

// cerr << srcplaneIMG << endl;
// cerr << imgplaneIMG << endl;

// namedWindow("1.picture one");//窗口显示图片
// imshow("srcplaneIMG", srcplaneIMG_8U);

// // namedWindow("2.picture two");
// imshow("imgplaneIMG", imgplaneIMG_8U);


// https://docs.opencv.org/3.4/d3/d50/group__imgproc__colormap.html
// Holds the colormap version of the image:
Mat srcplaneIMG_color, imgplaneIMG_color;
// Apply the colormap:
// https://docs.opencv.org/master/d3/d50/group__imgproc__colormap.html
// https://github.com/yuki-koyama/tinycolormap
applyColorMap(srcplaneIMG_8U, srcplaneIMG_color, COLORMAP_VIRIDIS);//COLORMAP_VIRIDIS
applyColorMap(imgplaneIMG_8U, imgplaneIMG_color, COLORMAP_VIRIDIS);
// Show the result:

// Mat srcplaneIMG_color_resize, imgplaneIMG_color_resize;
// Size size(imshowsize,imshowsize);//the dst image size,e.g.100x100
// resize(srcplaneIMG_color, srcplaneIMG_color_resize, size);
// resize(imgplaneIMG_color,imgplaneIMG_color_resize,size);

//https://blog.csdn.net/m0_37303351/article/details/78944904
//调整OpenCV弹出窗口大小
//https://blog.csdn.net/wuyoy520/article/details/47111295

// namedWindow("srcplaneIMG",WINDOW_NORMAL | WINDOW_KEEPRATIO | WINDOW_GUI_EXPANDED);
namedWindow("srcplaneIMG-Caustics", 0);
resizeWindow("srcplaneIMG-Caustics", imshowsize, imshowsize);
moveWindow("srcplaneIMG-Caustics", imshowsize+30, 0);
// resize(srcplaneIMG_color, outImage,  Size(256,256),INTER_NEAREST);
imshow("srcplaneIMG-Caustics", srcplaneIMG_8U);
// cvResizeWindow("srcplaneIMG", 320, 320);

// namedWindow("imgplaneIMG",WINDOW_NORMAL | WINDOW_KEEPRATIO | WINDOW_GUI_EXPANDED);
namedWindow("imgplaneIMG-Critical lines", 0);
resizeWindow("imgplaneIMG-Critical lines", imshowsize, imshowsize);
moveWindow("imgplaneIMG-Critical lines", 0, 0);
// resize(_inputImage, outImage,  Size(256,256),INTER_NEAREST);
imshow("imgplaneIMG-Critical lines", imgplaneIMG_color);
// cvResizeWindow("imgplaneIMG", 320, 320);

// https://docs.opencv.org/3.2.0/d7/dfc/group__highgui.html
// cout << "OpenCV version : " << CV_VERSION << endl;
// cout << "Major version : " << CV_MAJOR_VERSION << endl;
// cout << "Minor version : " << CV_MINOR_VERSION << endl;
// cout << "Subminor version : " << CV_SUBMINOR_VERSION << endl;
// // OpenCV version : 3.2.0
// Major version : 3
// Minor version : 2
// Subminor version : 0

// OpenCV version : 4.1.1
// Major version : 4
// Minor version : 1
// Subminor version : 1


// $ pkg-config --modversion opencv
// 3.2.0

// SetWindowPos(yourhwnd,0,0,0,newWidth,newHeight,SWP_NOMOVE|SWP_NOZORDER|SWP_NOACTIVATE);



// cout << "\nmassratio: " << massratio;

// cout << "mass: " << lens1.mass << "\npos: " << lens1.pos[0] <<","<< lens1.pos[1] << "\ndis: "<<lens1.dis;
// cout << "\n\nmass: " << lens2.mass << "\npos: " << lens2.pos[0] <<","<< lens2.pos[1] << "\ndis: "<<lens2.dis;
// cout << "\nmassratio: " << twolens.massratio;
// cout << "\nmu1: " << twolens.mu1;
// cout << "\nmu2: " << twolens.mu2;
// cout << "\nlens1.mass: " << twolens.lens1->mass;
// cout << "\nlens2.dis: " << twolens.lens2->dis;

// waitKey(5000);

// c++ saved file name create
// https://stackoverflow.com/questions/13108973/creating-file-names-automatically-c
// file name with system time
// https://stackoverflow.com/questions/22318389/pass-system-date-and-time-as-a-filename-in-c

        // Get the pressed value
        int key = (waitKey(0) & 0xFF);

        if (key == 's'){ //a
            // String name = format("img%04d.png", i++); // NEW !
            // imwrite(name, frame);
  String suffix;
  cout << "Please enter image suffix >>: ";
  cin >> suffix;
imwrite("savedimgs/srcplaneIMG_8U_"+suffix+".png", srcplaneIMG_8U); // A JPG FILE IS BEING SAVED
imwrite("savedimgs/imgplaneIMG_8U_"+suffix+".png", imgplaneIMG_8U); // A JPG FILE IS BEING SAVED
imwrite("savedimgs/srcplaneIMG_color_"+suffix+".png", srcplaneIMG_color); // A JPG FILE IS BEING SAVED
imwrite("savedimgs/imgplaneIMG_color_"+suffix+".png", imgplaneIMG_color); // A JPG FILE IS BEING SAVED
        }
        else if (key == 'q') return 0; // stop capturing by pressing q
        else {
            // Pressed an invalid key... continue with next frame
        }


// waitKey(0);

}