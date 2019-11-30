#include <iostream>
#include <math.h>
#include<cmath>
// #include<string>
// #include <array>

#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

//g++ -std=c++11 two_pointM_lens_inverse_ray_shoot.cpp -o two_pointM_lens_inverse_ray_shoot.out && ./two_pointM_lens_inverse_ray_shoot.out

//error: ISO C++ forbids declaration of ‘inverse_ray_shooting’ with no type [-fpermissive]
//https://stackoverflow.com/questions/7929477/explain-the-error-iso-c-forbids-declaration-of-personlist-with-no-type


const float PI = std::atan(1.0)*4;
const float arcsec2rad = PI/180/3600;
const float rad2arcsec = 180*3600/PI;


class Onelens
{
public:
    float mass, *pos, dis;
    Onelens(float mass, float *pos, float dis): mass(mass), pos(pos), dis(dis){}    
};

class Twolenses
{
public:
    Onelens *lens1 , *lens2;
    float massratio, mu1, mu2, *betax, *betay, *mag, m1, m2, beta;
public:
    Twolenses(Onelens*, Onelens*);
    void inverse_ray_shooting(float thetax, float thetay, float *betaxy);
    float comp_mag_samez(float thetax, float thetay);
    void get_imgs_lessmem(Mat srcplaneIMG, Mat srcplaneIMG_withoutlens, Mat imgplaneIMG, unsigned short int *ImgSize, float *xlim, float *ylim, unsigned short int raynum_1side);
};

Twolenses::Twolenses(Onelens *lens1, Onelens *lens2)
{
this -> lens1 = lens1;
this -> lens2  = lens2;
this -> massratio = lens1->mass/lens2->mass;
this -> mu1 = this -> massratio/(1+this -> massratio);
this -> mu2 = 1/(1+this -> massratio);
if (this -> lens2->dis != 0){
    float M1 = this -> mu1;
    float M2 = this -> mu2;
    float D1s = 1 - this -> lens1->dis;
    float D2 = this -> lens2->dis;
    float D2s = 1 - this -> lens2->dis;
    float D1 = this -> lens1->dis;
    this -> m1 = M1/(M1+D2s*D1/(D1s*D2)*M2);
    this -> m2 = 1 - this -> m1;
    this -> beta = (D2 - D1)/(D1s*D2); // D12*Ds/D1s/D2 from Erdl & Schneider 1993
}
}

void Twolenses::inverse_ray_shooting(float thetax, float thetay, float *betaxy)
{
    float r1_2 = pow(thetax - this->lens1->pos[0],2) + pow(thetay - this->lens1->pos[1],2);
    float r2_2 = pow(thetax - this->lens2->pos[0],2) + pow(thetay - this->lens2->pos[1],2);

    betaxy[0] = thetax - this->mu1*(thetax - this->lens1->pos[0])/r1_2 - this->mu2*(thetax - this->lens2->pos[0])/r2_2;
    betaxy[1] = thetay - this->mu1*(thetay - this->lens1->pos[1])/r1_2 - this->mu2*(thetay - this->lens2->pos[1])/r2_2;
    // return betaxy;
}

float Twolenses::comp_mag_samez(float thetax, float thetay)
{
    float mag,x1,y1,x2,y2,A11,A12,A22,detA,R1,R2;
    x1 = this->lens1->pos[0];
    y1 = this->lens1->pos[1];
    x2 = this->lens2->pos[0];
    y2 = this->lens2->pos[1];
    R1 = sqrt(pow(thetax-x1,2)+pow(thetay-y1,2));
    R2 = sqrt(pow(thetax-x2,2)+pow(thetay-y2,2));

    A11 = 1 - this->mu1*(1/pow(R1, 2) - 2*pow(thetax-x1,2)/pow(R1,4)) - this->mu2*(1/pow(R2, 2) - 2*pow(thetax-x2,2)/pow(R2,4));
    A12 = -2*this->mu1*(thetax-x1)*(thetay-y1)/pow(R1,4) - 2*this->mu2*(thetax-x2)*(thetay-y2)/pow(R2,4);
    A22 = 1 - this->mu1*(1/pow(R1, 2) - 2*pow(thetay-y1,2)/pow(R1,4)) - this->mu2*(1/pow(R2, 2) - 2*pow(thetay-y2,2)/pow(R2,4));
    detA = (A11*A22 - pow(A12,2));
    mag = 1/abs(detA);
    return mag;
}

void Twolenses::get_imgs_lessmem(Mat srcplaneIMG, Mat srcplaneIMG_withoutlens, Mat imgplaneIMG, unsigned short int *ImgSize, float *xlim, float *ylim, unsigned short int raynum_1side)
{
    float ratiox = (ImgSize[1]-1)/(xlim[1]-xlim[0]);
    float ratioy = (ImgSize[0]-1)/(ylim[1]-ylim[0]);
    float incx = (xlim[1]-xlim[0])/raynum_1side; // increment of x 
    float incy = (ylim[1]-ylim[0])/raynum_1side; // increment of y

    unsigned long int cnt = 1;

    float thetax = xlim[0];
    while(thetax<xlim[1]){

    float thetay = ylim[0];
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
    float betaxy[2] = {0, 0};
    // betax, betay = this -> inverse_ray_shooting(thetax, thetay, betaxy)
    // *betaxy = this -> inverse_ray_shooting(thetax, thetay, betaxy);
    this -> inverse_ray_shooting(thetax, thetay, betaxy);
    float mag = this -> comp_mag_samez(thetax, thetay);

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
    srcplaneIMG.at<uchar>(bj,bi) += mag;

}
      
    // srcplaneIMG_withoutlens[j, i]+=1
    if ((i>=0 && i<ImgSize[1])&&(j>=0&&j<ImgSize[0])){
    // *(srcplaneIMG_withoutlens.data + srcplaneIMG_withoutlens.step[0] * i + srcplaneIMG_withoutlens.step[1] * j + srcplaneIMG_withoutlens.elemSize1() * 0) += 1;
    srcplaneIMG_withoutlens.at<uchar>(j,i) += 1;
    // imgplaneIMG[j, i]+=mag    
    // *(imgplaneIMG.data + imgplaneIMG.step[0] * i + imgplaneIMG.step[1] * j + imgplaneIMG.elemSize1() * 0) += mag;
    imgplaneIMG.at<uchar>(j,i) += mag;
}

    thetay += incy;
    }
    thetax += incx;
}

// // srcplaneIMG /= srcplaneIMG_withoutlens;
// Mat result;
// result = srcplaneIMG.clone();
//     // int rowNumber = result.rows;
//     // int colNumber = result.cols;
//     for (int i = 0; i < ImgSize[0]; i++)
//     {
//         for (int j = 0; j < ImgSize[1]; j++){
//         srcplaneIMG.at<uchar>(j,i) = result.at<uchar>(j,i) / srcplaneIMG_withoutlens.at<uchar>(j,i);
//         }

//     }


// return void;
}


int main()
{
// using namespace PersonNamespace;
float pos1[2] = {-0.3, 0};
float mass1 = 1;
Onelens lens1(mass1, pos1 , 0);

float pos2[2] = {0.3, 0};
Onelens lens2(1.0, pos2 , 0);

float massratio = lens1.mass/lens2.mass;
float xlim[2] = {-2.5,2.5};
float ylim[2] = {-2.5,2.5};
unsigned short int ImgSize[2] = {512, 512}; // raw, colume
//unsigned short int    2 个字节   0 到 65,535

unsigned short int raynum_1side = 10000;

// THETAX, THETAY = genxy(xlim=xlim,ylim=ylim,num=num)
// def genxy(xlim=(-1,1),ylim=(-1,1),num=100):
//     x = np.linspace(xlim[0],xlim[1],num)
//     y = np.linspace(ylim[0],ylim[1],num)
//     X,Y = np.meshgrid(x,y)
//     # return X.reshape(1,-1), Y.reshape(1,-1) # shape: (1,10)
//     return X.reshape(1,-1)[0], Y.reshape(1,-1)[0] # shape: (3,)

Twolenses twolens( &lens1,  &lens2);

//Mat Matrix_name(行数，列数，存储元素的数据类型，每个矩阵点的通道数)
Mat srcplaneIMG(ImgSize[0], ImgSize[1], CV_16U, 1);
Mat srcplaneIMG_withoutlens(ImgSize[0], ImgSize[1], CV_16U, 1);
Mat imgplaneIMG(ImgSize[0], ImgSize[1], CV_16U, 1);

twolens.get_imgs_lessmem(srcplaneIMG, srcplaneIMG_withoutlens, imgplaneIMG, ImgSize, xlim, ylim, raynum_1side);

//cv::Vec3b vec3b      = img.at<cv::Vec3b>(0,0);
// uchar     vec3b0     = img.at<cv::Vec3b>(0,0)[0];
//https://blog.csdn.net/xiaowei_cqu/article/details/19839019 OpenCV】访问Mat中每个像素的值（新）

// namedWindow("1.picture one");//窗口显示图片
imshow("srcplaneIMG", srcplaneIMG);

// namedWindow("2.picture two");
imshow("imgplaneIMG", imgplaneIMG);


// cout << "\nmassratio: " << massratio;

// cout << "mass: " << lens1.mass << "\npos: " << lens1.pos[0] <<","<< lens1.pos[1] << "\ndis: "<<lens1.dis;
// cout << "\n\nmass: " << lens2.mass << "\npos: " << lens2.pos[0] <<","<< lens2.pos[1] << "\ndis: "<<lens2.dis;
// cout << "\nmassratio: " << twolens.massratio;
// cout << "\nmu1: " << twolens.mu1;
// cout << "\nmu2: " << twolens.mu2;
// cout << "\nlens1.mass: " << twolens.lens1->mass;
// cout << "\nlens2.dis: " << twolens.lens2->dis;

// waitKey(5000);
waitKey();
return 0;
}