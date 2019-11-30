#include <iostream>
#include <math.h>
// #include<string>
// #include <array>
using namespace std;

//g++ -std=c++11 two_pointM_lens_inverse_ray_shoot.cpp -o two_pointM_lens_inverse_ray_shoot.out && ./two_pointM_lens_inverse_ray_shoot.out

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



int main()
{

float pos1[2] = {-0.3, 0};
Onelens lens1(1.0, pos1 , 0);

float pos2[2] = {0.3, 0};
Onelens lens2(1.0, pos2 , 0);

Twolenses twolens( &lens1,  &lens2);

cout << "mass: " << lens1.mass << "\npos: " << lens1.pos[0] <<","<< lens1.pos[1] << "\ndis: "<<lens1.dis;
cout << "\n\nmass: " << lens2.mass << "\npos: " << lens2.pos[0] <<","<< lens2.pos[1] << "\ndis: "<<lens2.dis;
cout << "\nmassratio: " << twolens.massratio;
cout << "\nmu1: " << twolens.mu1;
cout << "\nmu2: " << twolens.mu2;
cout << "\nlens1.mass: " << twolens.lens1->mass;
cout << "\nlens2.dis: " << twolens.lens2->dis;
return 0;
}