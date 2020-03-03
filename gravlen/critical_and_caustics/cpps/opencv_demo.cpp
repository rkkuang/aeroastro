#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

//https://blog.csdn.net/luanpeng825485697/article/details/81181825
// 链接: https://pan.baidu.com/s/1NSiWERRIM0Wbd_oD1tvQqg  密码: uglw
// install opencv c++
// https://docs.opencv.org/master/d7/d9f/tutorial_linux_install.html

// sudo apt-get install libopencv-dev  -- not work

//https://docs.opencv.org/4.1.2/d7/d9f/tutorial_linux_install.html
// //opencv 依赖包lib查询 pkg-config --libs opencv

//https://www.jianshu.com/p/487a40f1add9

//https://www.cnblogs.com/woshijpf/p/3840530.html
//在linux环境下编译运行OpenCV程序的两种方法 *****

// gcc opencv_demo.cpp -o DisplayImage `pkg-config --cflags --libs opencv` -lstdc++
// g++ opencv_demo.cpp `pkg-config opencv --libs --cflags` -o DisplayImage


//https://www.cnblogs.com/blfshiye/p/4656902.html
//OpenCV（C++接口)学习笔记1-图像读取、显示、保存

// https://www.cnblogs.com/HL-space/p/10546604.html
// C++ Opencv Mat类型使用的几个注意事项及自写函数实现Laplace图像锐化


int main(void)
{
    Mat src_image_ = imread("bg.jpg");
    imshow("src_image_", src_image_);
    waitKey(5000);
    return 0;
}