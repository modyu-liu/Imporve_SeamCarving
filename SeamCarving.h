//
// Created by 刘思语 on 2023/6/18.
//
#include<opencv2/opencv.hpp>
#include<bits/stdc++.h>
#include"Color.h"
using namespace std;
using namespace cv;

const int INF = 1e8;


class SeamCarving{
private:
    Mat img , mask , emap , subimg , submask ;
    pair<int , int>end;
    vector<vector<double>>dp ;
    vector<vector<int>>fa;
    vector<vector<Point>>dis , ordinate;
    int side;
    vector<int>pos;


public:
    SeamCarving(Mat &img , Mat &mask);
    void get_long_boundary();
    bool add_seam();
    void calc_seam();
    void calc_e1();
    void Rect_wrap();
    void irrect_wrap();
    void show();
    double get_diff(int x1 , int y1 , int x2 , int y2);
    vector<vector<Point>> get_ordinate();

};