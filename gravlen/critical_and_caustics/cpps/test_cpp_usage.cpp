#include <iostream>
#include <math.h>
// #include<string>
// #include <array>
using namespace std;

//g++ -std=c++11 f1_first.cpp -o f1_first.out && ./f1_first.out
//https://blog.csdn.net/liitdar/article/details/81488329
//https://zhuanlan.zhihu.com/p/27306140

const float PI = std::atan(1.0)*4;
const float arcsec2rad = PI/180/3600;
const float rad2arcsec = 180*3600/PI;

// class Onelens
// {
// public:
//     float mass;   //
//     float *pos;  //
//     float dis;   //
// public:
//     // Onelens(float, float [2], float);
//     Onelens(float, float *, float);
// };

// Onelens::Onelens(float mass, float *pos, float dis)
// {
// this -> mass = mass;
// this -> pos = pos;
// this -> dis = dis;
// }

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
    float test;
    // Twolenses();
    // float massratio, mu1, mu2, *betax, *betay, *mag, m1, m2, beta;
public:
    Twolenses(Onelens*, Onelens*);
    // Twolenses(Onelens *lens1, Onelens *lens2);
};

Twolenses::Twolenses(Onelens *lens1, Onelens *lens2)
{
this -> lens1 = lens1;
this -> lens2  = lens2;
}

//a类创建的对象传递给b类的方法参数
//https://blog.csdn.net/yhblog/article/details/100777826

// c++ 类对象作为函数参数    
//https://blog.csdn.net/zjwson/article/details/56301449

//C++ 成员变量是别的类的对象
//https://blog.csdn.net/love9099/article/details/43059235
//https://blog.csdn.net/qq_35721743/article/details/83592415
// https://blog.csdn.net/Mr_xiao_1/article/details/80219162
// https://blog.csdn.net/HighKit/article/details/8752180
// https://blog.csdn.net/qq_43312665/article/details/88202001

// Twolenses::Twolenses(Onelens lens1, Onelens lens2)
// {
// this -> lens1 = lens1;
// this -> lens2  = lens2;
// // this -> massratio = lens1.mass/lens2.mass;
// // this -> mu1 = this -> massratio/(1+this -> massratio);
// // this -> mu2 = 1/(1+this -> massratio);
// // if (this -> lens2.dis != 0){
// //     float M1 = this -> mu1;
// //     float M2 = this -> mu2;
//     // float D1s = 1 - this -> lens1.dis;
//     // float D2 = this -> lens2.dis;
//     // float D2s = 1 - this -> lens2.dis;
//     // float D1 = this -> lens1.dis;
//     // this -> m1 = M1/(M1+D2s*D1/(D1s*D2)*M2);
//     // this -> m2 = 1 - this -> m1;
//     // this -> beta = (D2 - D1)/(D1s*D2); // D12*Ds/D1s/D2 from Erdl & Schneider 1993
// // }
// }



int main()
{

float pos1[2] = {-0.3, 0};
Onelens lens1(1.0, pos1 , 0);

float pos2[2] = {0.3, 0};
Onelens lens2(1.0, pos2 , 0);

// Twolenses twolens(*lens1, *lens2);
// Twolenses twolens;
// twolens.test = 1;

// Onelens *lens1 = new Onelens

// Twolenses twolens;
// twolens.test = 1;
// twolens.lens1 = &lens1;
// twolens.lens2 = &lens2;

Twolenses twolens( &lens1,  &lens2);

twolens.test = 1;


cout << "mass: " << lens1.mass << "\npos: " << lens1.pos[0] <<","<< lens1.pos[1] << "\ndis: "<<lens1.dis;
cout << "\n\nmass: " << lens2.mass << "\npos: " << lens2.pos[0] <<","<< lens2.pos[1] << "\ndis: "<<lens2.dis;
// cout << "massratio: " << twolens.massratio << "\n mu1: " << twolens.mu1 << "\n mu2: "<<twolens.mu2;

cout << "\nlens1: " << twolens.test;
cout << "\nlens1.mass: " << twolens.lens1->mass;
cout << "\nlens2.dis: " << twolens.lens2->dis;
return 0;
}