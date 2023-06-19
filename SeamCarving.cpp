//
// Created by 刘思语 on 2023/6/18.
//

#include "SeamCarving.h"

SeamCarving::SeamCarving(Mat &img , Mat &mask) {
    this->img = img;
    this->mask = mask;
    get_long_boundary();
}

void SeamCarving::get_long_boundary() {
    Side side ;
    pair<int , int>end;
    end.first = -1;end.second = -1;
    int res1 = 0 , res2 = 0;
    int maxn = 0 ;
    for(int i = 0 ; i < img.rows ; i++){
        if(mask.at<uchar>(i , 0) == BG){
            res1++;
        }
        else {
            if(maxn < res1){
                maxn = res1;
                side = LEFT;
                end.first = i - maxn;
                end.second = i - 1;
            }
            res1 = 0;

        }
        if(mask.at<uchar>(i , mask.cols - 1) == BG){
            res2++;
        }
        else {
            if(maxn < res2){
                maxn = res2 ;
                side = RIGHT;

                end.first = i - maxn;
                end.second = i - 1;
            }
            res2 = 0;

        }
    }
    if(maxn < res1){
        maxn = res1;
        side = LEFT;

        end.first = mask.rows - maxn;
        end.second = mask.rows - 1;
    }
    if(maxn < res2){
        maxn = res2;
        side = RIGHT;
        end.first = mask.rows - maxn;
        end.second = mask.rows - 1;

    }
    res1 = 0 , res2 = 0;

    for(int i = 0 ; i < img.cols ; i++){
        if(mask.at<uchar>(0 , i) == BG){
            res1++;
        }
        else {
            if(maxn < res1){
                maxn = res1;
                side = TOP;

                end.first = i - maxn;
                end.second = i - 1;
            }
            res1 = 0;

        }
        if(mask.at<uchar>(mask.rows - 1 , i) == BG){
            res2++;
        }
        else {
            //cout<<"find::::"<<i<<'\n';
            //cout<<"ok??"<<' ' << (int)mask.at<uchar>(mask.rows - 1 , i)<<'\n';
            //cout<<(int)img.at<Vec3b>(img.rows - 1, i)[0]<<' '<<(int)img.at<Vec3b>(img.rows - 1, i)[1]<<' ' << (int)img.at<Vec3b>(img.rows - 1, i)[2]<<'\n';


            if(maxn < res2){
                maxn = res2 ;
                side = BOTTOM;
                end.first = i - maxn;
                end.second = i - 1;
            }
            res2 = 0;

        }
    }
    if(maxn < res1){
        maxn = res1;
        side = TOP;

        end.first = mask.cols - maxn;
        end.second = mask.cols - 1;
    }
    if(maxn < res2){
        maxn = res2;
        side = BOTTOM;
        end.first = mask.cols - maxn;
        end.second = mask.cols - 1;
    }
    this->end = end;
    this->side = side;
    //show();
}

//g -> add_tweights( 0,   /* capacities */  1, 5 );
//g -> add_tweights( 1,   /* capacities */  2, 6 );
//g -> add_edge( 0, 1,    /* capacities */  3, 4 );

void SeamCarving::build_graph() {
    if(this->side == LEFT || this->side == RIGHT){
        int node_num = (this->end.second - this->end.first + 1) * img.cols;
        g = new GraphType(node_num , node_num * 4);
        for(int i = this->end.first ; i <= this->end.second ; i++){
            for(int j = 0; j < img.cols ; j++){
                g->add_node();
            }
        }
        for(int i = this->end.first ; i <= this->end.second ; i++){
            g->add_tweights((i - end.first) * img.cols , INF , 0);
            g->add_tweights((i - end.first) * img.cols + img.cols - 1 , 0 , INF);
        }
        for(int i = this->end.first ; i <= this->end.second ; i++ ){
            for(int j = 0 ; j < img.cols ; j++){

            }
        }


    }
}

bool SeamCarving::add_seam(){
    get_long_boundary();
    if(this->end.first == -1)return true;
    build_graph();
}
void SeamCarving::show(){
    Mat pre = img.clone();
    cout<<"111"<<' '<<end.first<<' '<<end.second<<' ' << this->side<< '\n';
    //cout<<(int)pre.at<Vec3b>(img.rows - 1, 4182)[0]<<' '<<(int)pre.at<Vec3b>(img.rows - 1, 4182)[1]<<' ' << (int)pre.at<Vec3b>(img.rows - 1, 4182)[2]<<'\n';
    //cout<<(int)img.at<Vec3b>(img.rows - 1, 4182)[0]<<' '<<(int)img.at<Vec3b>(img.rows - 1, 4182)[1]<<' ' << (int)img.at<Vec3b>(img.rows - 1, 4182)[2]<<'\n';
    //cout<<(int)mask.at<uchar>(img.rows - 1 , 4182)<<'\n';

    if(this->side == LEFT){

        for(int i = end.first ; i <= end.second ; i++){
            pre.at<Vec3b>(i , 0) = BLUE;
        }
    }
    else if(this->side == RIGHT){

        for(int i = end.first ; i <= end.second ; i++){
            pre.at<Vec3b>(i , pre.cols - 1) = BLUE;
        }
    }
    else if(this->side == TOP){

        for(int i = end.first ; i <= end.second ; i++){
            pre.at<Vec3b>(0 , i) = BLUE;
        }
    }
    else if(this->side == BOTTOM){

        for(int i = end.first ; i <= end.second ; i++){
            pre.at<Vec3b>(pre.rows - 1 , i) = BLUE;
        }
    }
    cout<<"inthis"<<'\n';

    imshow("image" , pre);
    waitKey(0);

}