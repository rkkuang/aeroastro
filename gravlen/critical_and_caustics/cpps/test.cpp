#include<iostream>
using namespace std;
int main(){
    int a=123;
    //&a表示a在内存中的地址，也就是123在内存中的地址
    cout<<"a: "<<a<<endl<<"a's address："<<&a<<endl;
    //此时p是一个指针，指向a所在的位置
    int *p=&a;
    cout<<"p: "<<p<<endl;
    //声明p之后，在p之前添加*表示p指向内存的值
    cout<<"p's value: "<<*p<<endl;
    //同时p也是 一个变量，在内存中也有一个地址储存它，但其地址不是a的地址
    cout<<"p's address: "<<&p<<endl;
    //试试*&组合使用是什么效果
    cout<<"*&p: "<<*&p<<endl;
    //&p是一个内存地址，*&p表示&p指向地址内存空间的值,在这里表示a的地址
    cout<<"**&p: "<<**&p<<endl;
    //刚才我们已经知道*&p是a的地址，那么**&p就表示a的值
return 0;}