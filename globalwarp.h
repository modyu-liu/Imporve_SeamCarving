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
struct Pointd{
    double x , y;
    Pointd(double x , double y):x(x) , y(y){}
    Pointd(int x , int y):x(x) , y(y){}
    Pointd(Point p) : x(p.x) , y(p.y){}
    Pointd(){}
    Pointd &operator+=(const Pointd &p) {
        x += p.x, y += p.y;
        return *this;
    }
    Pointd &operator-=(const Pointd &p) {
        x -= p.x, y -= p.y;
        return *this;
    }
    friend Pointd operator-(const Pointd &p) {
        return Pointd(-p.x, -p.y);
    }
    friend Pointd operator+(Pointd lhs, const Pointd &rhs) {
        return lhs += rhs;
    }
    friend Pointd operator-(Pointd lhs, const Pointd &rhs) {
        return lhs -= rhs;
    }
};
class globalwarp{
private:
    Mat img ;
    vector<vector<Point>>mordinate;
    vector<pair<Pointd , Pointd>>line;
    vector<vector<vector<pair<Pointd , Pointd>>>>seg_line;
    vector<vector<vector<int>>>allocate;
    vector<double>rotate;
    vector<vector<vector<pair<MatrixXd , MatrixXd>>>>seg_line_w;
    vector<vector<Point>>ver;

    int seg_num;


public:
    globalwarp(Mat img , vector<vector<Point>>mordinate);
    void get_line();
    bool inmesh(Pointd& s , int x ,int y);
    Pointd get_intersection(Pointd& a , Pointd& b , Pointd& c , Pointd& d);
    bool is_intersection(Pointd& a , Pointd& b , Pointd& c , Pointd& d);
    void mseg_line();
    Pointd modify_line(Pointd& s , int i , int j);

    void init_rotate();
    MatrixXd get_single_shape_preservation(int x ,int y);
    SparseMatrix<double> get_shape_preservation();
    MatrixXd inv_biliner(Pointd& P ,  int x , int y);
    SparseMatrix<double> get_line_preservation();
    pair<SparseMatrix<double> , VectorXd> get_boundary_constraints();
    SparseMatrix<double> get_position_information();
    void update_rotate();
    SparseMatrix<double> Connect_mat(SparseMatrix<double>& m1 , SparseMatrix<double>& m2);
    void start_learn();
    bool in_line(Vector2d &p , Vector2d &a , Vector2d &b);
    vector<vector<Point>> get_ordinate();
};