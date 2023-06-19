//
// Created by 刘思语 on 2023/6/18.
//

#include "SeamCarving.h"

SeamCarving::SeamCarving(Mat &img , Mat &mask) {
    this->img = img;
    this->mask = mask;
    while(!add_seam());
    imshow("img" , img);
    waitKey(0);

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
            int now = (i - end.first) * img.cols;
            if(j < img.cols - 1){
                if(j == img.cols - 2){
                    g->add_edge(now , now + 1 , 0 , INF);
                }
                else {
                    Vec3d diff = ((Vec3d)img.at<Vec3b>(i , j) - (Vec3d)img.at<Vec3b>(i , j + 2));
                    double res = diff.dot(diff);
                    g->add_edge(now , now + 1 , res , INF);
                }
            }
            if(i < this->end.second){
                int idx = now + img.cols;
                if(j == 0) {
                    g->add_edge(now, idx, 0, 0);
                }
                else {
                    Vec3d diff = ((Vec3d)img.at<Vec3b>(i + 1, j) - (Vec3d)img.at<Vec3b>(i , j - 1));
                    double cnt1 = diff.dot(diff);
                    diff = ((Vec3d)img.at<Vec3b>(i , j) - (Vec3d)img.at<Vec3b>(i + 1 , j - 1));
                    double cnt2 = diff.dot(diff);
                    g->add_edge(now , idx , cnt1 , cnt2);
                }
            }
            {
                if(j == img.cols - 1 || i == end.second)continue;
                int idx1 = now + 1;
                int idx2 = now + img.cols;
                g->add_edge(idx1 , idx2 , INF , 0);
            }
            {
                if(j == img.cols - 1 || i == end.second)continue;
                int idx1 = now;
                int idx2 = now + img.cols + 1;
                g->add_edge(idx2 , idx1 , INF , 0);

            }
        }
    }

}
void SeamCarving::segment(){
    g->maxflow();
    this->pos.clear();
    for(int i = end.first; i <= end.second ; i++){
        for(int j = 0; j < img.cols - 1 ; j++){
            int now = (i - end.first) * img.cols + j;
            if(g->what_segment(now) != g->what_segment(now + 1)){
                this.pos.emplace_back(j);
                break;
            }
        }
    }
}
bool SeamCarving::add_seam(){
    get_long_boundary();
    if(this->end.first == -1)return true;
    if(this->side == TOP || this->side == BOTTOM){
        transpose(img, img);
        transpose(mask , mask);
    }
    build_graph();
    for(int i = end.first ; i <= end.second ; i++){
        int p = i - end.first;

        if(this->side == TOP || this->side == LEFT){
            for(int j = 0 ; j < p ; j ++){
                img.at<Vec3b>(i , j) = img.at<Vec3b>(i , j + 1);
                mask.at<uchar>(i , j) = mask.at<uchar>(i , j + 1);
            }

        }
        else {
            for(int j = img.cols - 1 ; j > p ; j --){
                img.at<Vec3b>(i , j) = img.at<Vec3b>(i , j - 1);
                mask.at<uchar>(i , j) = mask.at<uchar>(i , j - 1);
            }
        }
        int cnt =0 ;
        if(p == 0 || p == img.cols - 1){
            if(p == 0 ){
                img.at<Vec3b>(i , p) = img.at<Vec3b>(i , p + 1);
            }
            else {
                img.at<Vec3b>(i , p) = img.at<Vec3b>(i , p - 1);

            }
        }
        else {
            img.at<Vec3b>(i , p) = (img.at<Vec3b>(i , p - 1) + img.at<Vec3b>(i , p + 1)) / 2;
        }
        mask.at<uchar>(i , p) = ADT;

    }
    return false;

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