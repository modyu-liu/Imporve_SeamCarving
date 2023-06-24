//
// Created by 刘思语 on 2023/6/18.
//
#include <opencv2/opencv.hpp>
using namespace cv;

#ifndef RECTPANORAMIC_COLOR_H
#define RECTPANORAMIC_COLOR_H

#endif //RECTPANORAMIC_COLOR_H
const Vec3b BLUE = Vec3b(255,0,0); //Background
const Vec3b GREEN = Vec3b(0,255,0); //Foreground
const Vec3b LIGHTBLUE = Vec3b(255,255,160); //ProbBackground
const Vec3b PINK = Vec3b(230,130,255); //ProbBackground
const Vec3b RED = Vec3b(0,0,255); //color of Rectangle
const Vec3b WHITE = Vec3b(255,255,255);
const Vec3b BLACK = Vec3b(0 , 0 , 0);

enum Side{
    TOP = 0 ,
    BOTTOM = 1,
    LEFT = 2,
    RIGHT = 3,
};
enum State{
    BG = 0 ,
    FG = 1 ,
    ADT = 2,
};