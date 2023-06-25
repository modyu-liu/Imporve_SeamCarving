//
// Created by 刘思语 on 2023/6/18.
//
#include<opencv2/opencv.hpp>
#include<bits/stdc++.h>
#include"Color.h"
#include<sys/time.h>
using namespace std;
using namespace cv;

const int INF = 1e8;


class SeamCarving{
private:
    Mat img , mask , emap , grayimg;
    pair<int , int>end;
    vector<vector<double>>dp ;
    vector<vector<int>>fa;

    vector<vector<Point>>dis , ordinate;
    int side;
    vector<int>pos;

public:
    clock_t time1 , time2;

    double time = 0;
    SeamCarving(Mat &img , Mat &mask);
    void get_long_boundary();
    bool add_seam();
    void calc_seam();
    void calc_e1();
    void Rect_wrap();
    void irrect_wrap();
    vector<vector<Point>> get_ordinate();

};