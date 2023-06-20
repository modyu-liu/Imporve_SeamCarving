#include<bits/stdc++.h>
#include<opencv2/opencv.hpp>
#include "SeamCarving.h"

using namespace std;
using namespace cv;
void fillHole(Mat src , Mat &dst){
    cv::Size m_Size = src.size();
    Mat Temp = Mat::zeros(m_Size.height + 2, m_Size.width + 2, src.type());
    src.copyTo(Temp(cv::Range(1, m_Size.height + 1), cv::Range(1, m_Size.width + 1)));

    Rect ccomp;
    cv::floodFill(Temp, cv::Point(0, 0), cv::Scalar(255) , &ccomp ,cv::Scalar(5));
    Mat cutImg;
    Temp(cv::Range(1, m_Size.height + 1), cv::Range(1, m_Size.width + 1)).copyTo(cutImg);
    dst = src | (~cutImg);
}
Mat clean(Mat src ){
    Mat bw;
    cv::cvtColor(src, src, cv::COLOR_BGR2GRAY);
    uchar thr = 252;
    Mat mask = Mat::zeros(src.size(), CV_8UC1);
    for (int row = 0; row < src.rows; row++) {
        for (int col = 0; col < src.cols; col++) {
            if (src.at<uchar>(row, col)<thr) {
                mask.at<uchar>(row, col) = 255;
            }
        }
    }
    fillHole(mask, bw);
    bw = ~bw;

    Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(5, 5));
    Mat dilate_out;
    cv::dilate(bw, dilate_out, element);
    cv::dilate(dilate_out, dilate_out, element);
    cv::dilate(dilate_out, dilate_out, element);
    Mat erode_out;
    erode(dilate_out, erode_out, element);
    return erode_out;
}
int main(){
    Mat img ;
    //img = imread("2.png");
    img = imread("1.jpg");
    //img = imread("10a_input.jpg");
    //clean_img(img);
    Mat mask = clean(img);
    for(int i = 0; i < mask.rows ; i++){
        for(int j = 0 ; j < mask.cols ; j++){
            if(mask.at<uchar>(i , j) == 255){
                mask.at<uchar>(i , j) = BG;
            }
            else {
                mask.at<uchar>(i , j) = FG;
            }

        }
    }
    auto prepare = [&](){
        int up = 0 , bottom = 0 , left = 0 , right = 0;
        bool ok = 1;
        while(ok) {
            for (int i = 0; i < img.cols; i++) {
                if (mask.at<uchar>(up , i) != BG) {
                    ok = 0;
                    break;
                }
            }
            if(ok)up++;
        }
        ok = 1;
        while(ok) {
            for (int i = 0; i < img.cols; i++) {
                if (mask.at<uchar>(img.rows - 1 - bottom , i) != BG) {
                    ok = 0;
                    break;
                }
            }
            if(ok)bottom++;
        }
        ok = 1;
        while(ok) {
            for (int i = 0; i < img.rows; i++) {
                if (mask.at<uchar>(i, left) != BG) {
                    ok = 0;
                    break;
                }
            }
            if(ok)left++;
        }
        ok = 1;
        while(ok) {
            for (int i = 0; i < img.rows; i++) {
                if (mask.at<uchar>(i, img.cols - 1 - right) != BG) {
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
        mask = mask(roi).clone();
        //imshow("img" , subimg);
        //waitKey(0);

    };

    prepare();

    SeamCarving seam(img , mask);
    return 0 ;
}