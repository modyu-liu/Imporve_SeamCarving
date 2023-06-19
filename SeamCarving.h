//
// Created by 刘思语 on 2023/6/18.
//
#include<opencv2/opencv.hpp>
#include<bits/stdc++.h>
#include"Color.h"
#include "graph.h"
using namespace std;
using namespace cv;

const int INF = 1e8;
typedef Graph<double , double ,double > GraphType;
class SeamCarving{
private:
    Mat img;
    Mat mask;
    pair<int , int>end;
    int side;
    GraphType *g ;


public:
    SeamCarving(Mat &img , Mat &mask);
    void get_long_boundary();
    bool add_seam();
    void build_graph();
    void show();


};