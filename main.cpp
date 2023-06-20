#include<bits/stdc++.h>
#include<opencv2/opencv.hpp>
#include "SeamCarving.h"

using namespace std;
using namespace cv;

int main(){
    Mat img ;
    //img = imread("2.png");
    img = imread("10a_input.jpg");
    //clean_img(img);



    auto prepare = [&](){
        int up = 0 , bottom = 0 , left = 0 , right = 0;
        bool ok = 1;
        while(ok) {
            for (int i = 0; i < img.cols; i++) {
                if (img.at<Vec3b>(up, i) != WHITE) {
                    ok = 0;
                    break;
                }
            }
            if(ok)up++;
        }
        ok = 1;
        while(ok) {
            for (int i = 0; i < img.cols; i++) {
                if (img.at<Vec3b>(img.rows - 1 - bottom , i) != WHITE) {
                    ok = 0;
                    break;
                }
            }
            if(ok)bottom++;
        }
        ok = 1;
        while(ok) {
            for (int i = 0; i < img.rows; i++) {
                if (img.at<Vec3b>(i, left) != WHITE) {
                    ok = 0;
                    break;
                }
            }
            if(ok)left++;
        }
        ok = 1;
        while(ok) {
            for (int i = 0; i < img.rows; i++) {
                if (img.at<Vec3b>(i, img.cols - 1 - right) != WHITE) {
                    ok = 0;
                    break;
                }
            }
            if(ok)right++;
        }
        //imshow("img" , img);
        //waitKey(0);
        //Rect roi(1000 , 0 , 800 , 300);

        Rect roi(left , up , img.cols - left - right , img.rows - up - bottom);
        img = img(roi).clone();

        //imshow("img" , subimg);
        //waitKey(0);

    };

    prepare();
    Mat mask(img.rows , img.cols , CV_8UC1);
    for(int i = 0; i < img.rows ; i ++){
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