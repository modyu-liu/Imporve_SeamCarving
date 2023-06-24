//
// Created by 刘思语 on 2023/6/18.
//

#include "SeamCarving.h"

SeamCarving::SeamCarving(Mat &img , Mat &mask) {
    this->img = img;
    this->mask = mask;
    this->dp = vector<vector<double>>(max(img.rows , img.cols) , vector<double>(max(img.rows , img.cols) , 0));
    this->fa = vector<vector<int>>(max(img.rows , img.cols) , vector<int>(max(img.rows , img.cols) , 0));
    this->dis = vector<vector<Point>>(img.rows , vector<Point>(img.cols , Point(0 , 0)));
    add_seam();
}
//使用e1能量函数
void SeamCarving::calc_e1() {
    Mat gx, gy;
    Sobel(this->subimg, gx, CV_32F, 1, 0);
    Sobel(this->subimg, gy, CV_32F, 0, 1);
    magnitude(gx, gy, this->emap);  // 只计算合梯度的幅值
    for(int i = 0; i < subimg.rows ; i++){
        for(int j = 0; j < subimg.cols ; j++){
            if(submask.at<uchar>(i , j) == BG){
                this->emap.at<float>(i , j) = INF;
            }
        }
    }


}

void SeamCarving::get_long_boundary() {
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
}

double SeamCarving::get_diff(int x1 , int y1 , int x2 , int y2){
    double R = (int)subimg.at<Vec3b>(x1 , y1)[0] - (int)subimg.at<Vec3b>(x2 , y2)[0];
    double G = (int)subimg.at<Vec3b>(x1 , y1)[1] - (int)subimg.at<Vec3b>(x2 , y2)[1];
    double B = (int)subimg.at<Vec3b>(x1 , y1)[2] - (int)subimg.at<Vec3b>(x2 , y2)[2];
    return sqrt(R * R + G * G + B * B);
}

void SeamCarving::calc_seam(){
    //Mat grayimg ;
    //cvtColor(subimg, grayimg, COLOR_BGR2GRAY);
    for(int i = 0 ; i < subimg.rows ; i++){
        if(i == 0){
            for(int j = 0; j < subimg.cols ;j++){
                dp[i][j] = emap.at<float>(i , j);
            }
        }
        else {
            for (int j = 0; j < subimg.cols; j++) {
                dp[i][j] = 0;
                dp[i][j] += emap.at<float>(i, j);
                if(j == 0){
                    double CU = 0 , CR = 0 ;
                    CU = dp[i - 1][j];
                    CR = dp[i - 1][j + 1] + get_diff(i - 1, j , i , j + 1);
                    if(CU < CR){
                        dp[i][j] += CU;
                        fa[i][j] = 0;
                    }
                    else {
                        dp[i][j] += CR;
                        fa[i][j] = 1;
                    }
                }
                else if(j == subimg.cols - 1){
                    double CL = 0 , CU = 0 ;
                    CL = dp[i - 1][j - 1] + get_diff(i - 1 , j , i , j - 1);
                    CU = dp[i - 1][j];
                    if(CL < CU){
                        dp[i][j] += CL;
                        fa[i][j] = -1;
                    }
                    else {
                        dp[i][j] += CU;
                        fa[i][j] = 0;
                    }
                }
                else {
                    double CL = 0, CU = 0, CR = 0;
                    CL += dp[i - 1][j - 1];
                    CL += get_diff(i , j - 1 , i , j + 1);
                    CL += get_diff(i - 1 , j , i , j - 1);
                    CU += dp[i - 1][j];
                    CU += get_diff(i , j - 1 , i , j + 1);
                    CR += dp[i - 1][j + 1];
                    CR += get_diff(i , j - 1 , i , j + 1);
                    CR += get_diff(i - 1 , j , i , j + 1);
                    double minn = min({CL , CU , CR});
//                    cout<<"find::"<<' '<<CL<<' '<<CU<<' '<<CR<<'\n';
                    dp[i][j] += minn;
                    if(minn == CL )fa[i][j] = -1;
                    else if(minn == CU) fa[i][j] = 0;
                    else fa[i][j] = 1;
                }

            }
        }

    }
    int col = 0;
    double minn = 1e18;
    for(int j = 0; j < subimg.cols ; j++){
        if(dp[subimg.rows - 1][j] < minn){
            minn = dp[subimg.rows - 1][j];
            col = j;
        }
    }
    this->pos.clear();
    int row = subimg.rows - 1;
    bool inf = 0;
    if(dp[row][col] >= INF){
        inf = 1;
        cout<<"check::"<<dp[row][col]<<' ' << col<<'\n';
        cout<<"Error!"<<'\n';
    }
    bool ok = 1;
    while(row){
        this->pos.emplace_back(col);
        if(submask.at<uchar>(row , col) == BG){
            ok = 0;
        }
        col += fa[row][col];
        row--;
        if(col >= subimg.cols || col < 0){
            cout<<"dp error!"<<'\n';
            cout<<"col out image!"<<'\n';
            cout<<"img size "<< subimg.rows <<' '<<subimg.cols<<'\n';

            cout<<"ok??? "<<row<<' '<<col<<'\n';


            exit(0);
        }
    }
    this->pos.emplace_back(col);

    reverse(this->pos.begin() , this->pos.end());


}

bool SeamCarving::add_seam(){

    while(1) {
        get_long_boundary();
        if (this->end.first == -1 || this->end.second == this->end.first )break;
        if (this->side == TOP || this->side == BOTTOM) {
            transpose(img, img);
            transpose(mask, mask);
        }
        Rect roiRect(0, end.first, img.cols, end.second - end.first + 1);
        this->subimg = this->img(roiRect).clone();
        this->submask = this->mask(roiRect).clone();
        calc_e1();
        calc_seam();
        this->subimg.copyTo(this->img(roiRect));
        //continue;
        for (int i = end.first; i <= end.second; i++) {
            int p = i - end.first;
            if (this->side == TOP || this->side == LEFT) {
                for (int j = 0; j < this->pos[p]; j++) {
                    img.at<Vec3b>(i, j) = img.at<Vec3b>(i, j + 1);
                    mask.at<uchar>(i, j) = mask.at<uchar>(i, j + 1);
                }
            } else {
                for (int j = img.cols - 1; j > this->pos[p]; j--) {
                    img.at<Vec3b>(i, j) = img.at<Vec3b>(i, j - 1);
                    mask.at<uchar>(i, j) = mask.at<uchar>(i, j - 1);
                }
            }
            int cnt = 0;
            if (this->pos[p] == 0 || this->pos[p] == img.cols - 1) {
                if (this->pos[p] == 0) {
                    img.at<Vec3b>(i, this->pos[p]) = img.at<Vec3b>(i, this->pos[p] + 1);
                } else {
                    img.at<Vec3b>(i, this->pos[p]) = img.at<Vec3b>(i, this->pos[p] - 1);
                }
            } else {
                img.at<Vec3b>(i, this->pos[p]) = img.at<Vec3b>(i, this->pos[p] - 1) / 2 + img.at<Vec3b>(i, this->pos[p] + 1) / 2;
                Vec3b rev(255 - img.at<Vec3b>(i, this->pos[p])[0] , 255 - img.at<Vec3b>(i, this->pos[p])[1] , 255 - img.at<Vec3b>(i, this->pos[p])[2]);
                img.at<Vec3b>(i, this->pos[p]) = rev;
            }
            if(mask.at<uchar>(i , this->pos[p]) != BG) {
                mask.at<uchar>(i, this->pos[p]) = ADT;
            }
        }
        //cout << "point3" << '\n';
        if (this->side == TOP || this->side == BOTTOM) {
            transpose(img, img);
            transpose(mask, mask);
        }
        //imshow("img" , img);
        //waitKey(1);
       if(side == LEFT){
           for(int i = end.first ; i <= end.second ; i++){
               for(int j = 0 ; j < this->pos[i - end.first]; j ++){
                   dis[i][j].y = dis[i][j + 1].y - 1;
                   dis[i][j].x = dis[i][j + 1].x;
               }
           }
       }
       else if(side == RIGHT){
           for(int i = end.first ; i <= end.second ; i++){
               for(int j = img.cols - 1 ; j > this->pos[i - end.first]; j --){
                   dis[i][j].y = dis[i][j - 1].y + 1;
                   dis[i][j].x = dis[i][j - 1].x;
               }
           }
       }
       else if(side == TOP){
           for(int i = end.first ; i <= end.second ; i++){
               for(int j = 0 ; j < this->pos[i - end.first]; j ++){
                   dis[j][i].y = dis[j + 1][i].y ;
                   dis[j][i].x = dis[j + 1][i].x - 1;
               }
           }
       }
       else if(side = BOTTOM){
           for(int i = end.first ; i <= end.second ; i++){
               for(int j = img.rows - 1 ; j > this->pos[i - end.first]; j --){
                   dis[j][i].y = dis[j - 1][i].y ;
                   dis[j][i].x = dis[j - 1][i].x + 1;
               }
           }
       }
    }

    Rect_wrap();

}
void SeamCarving::Rect_wrap(){
    /* a = x / y
     * xy = 400
     * ay^2 = 400
     * y = sqrt(400/a)
     * */
    double a = (double)img.rows / img.cols;
    //cout<<"check::"<<a<<'\n';

    int cols = sqrt(400.0 / a);

    int rows = 400 / cols;
    //cout<<"find::"<<rows<<' '<<cols<<'\n';

    int row_num = img.rows / rows;
    int col_num = img.cols / cols;

    //cout<<"ppp::"<<row_num<<' '<<col_num<<'\n';
    vector<vector<Point>>res;

    for(int i = 0 ;i < img.rows ; i += row_num){
        vector<Point>pre;
        for(int j = 0; j < img.cols ; j += col_num ){
            pre.emplace_back(Point(i , j));
        }
        if(pre.back().y != img.cols - 1){
            pre.emplace_back(Point(i , img.cols - 1));
        }
        res.emplace_back(pre);
    }

    if(res.back().back().x != img.rows - 1){
        vector<Point>pre;
        for(int j = 0; j < img.cols ; j += col_num ){
            pre.emplace_back(Point(img.rows - 1 , j));
        }
        if(pre.back().y != img.cols - 1){
            pre.emplace_back(Point(img.rows - 1 , img.cols - 1));
        }
        res.emplace_back(pre);
    }
    /*
    for(int i = 0 ; i < res.size() ;i ++){
        for(auto it : res[i]){
            circle(img, it, 2, Scalar(0, 255, 0), 2);
        }
    }
    imshow("img" , img);
    waitKey(0);
    */
    this->ordinate = res;
    irrect_wrap();

}
void SeamCarving::irrect_wrap(){
    for(int i = 0; i < this->ordinate.size() ; i++){
        for(int j = 0; j < this->ordinate[i].size() ; j++){
            int x = this->ordinate[i][j].x;
            int y = this->ordinate[i][j].y;
            int new_x = x - dis[x][y].x;
            int new_y = y - dis[x][y].y;
            this->ordinate[i][j] = Point(new_x , new_y);
        }
    }
    /*
    Mat pre = Mat(img.rows , img.cols , CV_8UC3  ,Scalar(255,255,255));
    for(int i = 0; i < img.rows ; i++){
        for(int j = 0; j < img.cols ; j++){
            if(mask.at<uchar>(i , j) == BG || mask.at<uchar>(i , j) == ADT){
                continue;
            }
            else {
                int x = i - dis[i][j].x;
                int y = j - dis[i][j].y;
                pre.at<Vec3b>(x , y) = img.at<Vec3b>(i , j);
            }
        }
    }
    for(int i = 0; i < this->ordinate.size() ; i++){
        for(int j = 0; j < this->ordinate[i].size() ; j++){
            Point now = Point(this->ordinate[i][j].y , this->ordinate[i][j].x);
            circle(pre , now , 2 , Scalar(0, 255, 0) , 2);
        }
    }
    imshow("img" , pre);
    waitKey(0);
    */
}
vector<vector<Point>> SeamCarving::get_ordinate() {
    return this->ordinate;
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