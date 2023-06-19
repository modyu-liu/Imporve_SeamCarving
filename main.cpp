#include<bits/stdc++.h>
#include<opencv2/opencv.hpp>
#include "SeamCarving.h"
using namespace std;
using namespace cv;

int main(){
    Mat img ;
    img = imread("/Users/liusiyu/CLionProjects/第二轮考核项目/RectPanoramic/data/1_input.jpg");
    Mat mask(img.rows , img.cols , CV_8UC1);
    for(int i = 0 ;i < img.rows ; i++){
        for(int j = 0; j < img.cols ; j++){
            if(img.at<Vec3b>(i , j) == Vec3b(255 , 255 , 255)){
                mask.at<uchar>(i , j) = BG;
            }
            else {
                mask.at<uchar>(i , j) = FG;
            }
        }
    }
    SeamCarving seam(img , mask);
    return 0 ;
}