//
// Created by 刘思语 on 2023/6/21.
//
#include <bits/stdc++.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

//#include<bits/stdc++.h>
#include<opencv2/opencv.hpp>

#include "lsd.h"
using namespace cv;
using namespace std;
using namespace Eigen;


const double PI = acos(-1);

class globalwarp{
private:
    Mat img ;
    vector<vector<Point>>mordinate;
    vector<pair<Point , Point>>line;
    vector<vector<vector<pair<Point , Point>>>>seg_line;
    vector<vector<vector<int>>>allocate;
    vector<double>rotate;
    int seg_num;


public:
    globalwarp(Mat img , vector<vector<Point>>mordinate);
    void get_line();
    bool inmesh(Point s , int x ,int y);
    Point get_intersection(Point a , Point b , Point c , Point d);
    bool is_intersection(Point a , Point b , Point c , Point d);
    void mseg_line();
    void init_rotate();
    MatrixXd get_single_shape_preservation(int x ,int y);
    SparseMatrix<double> get_shape_preservation();
    MatrixXd inv_biliner(Point P ,  int x , int y);
    SparseMatrix<double> get_line_preservation();
    pair<SparseMatrix<double> , VectorXd> get_boundary_constraints();
    SparseMatrix<double> get_position_information();
    void update_rotate(VectorXd &V);
    SparseMatrix<double> Connect_mat(SparseMatrix<double> m1 , SparseMatrix<double> m2);
    void start_learn();
    bool in_line(Vector2d &p , Vector2d &a , Vector2d &b);
    void show();

};