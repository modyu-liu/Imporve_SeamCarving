//
// Created by 刘思语 on 2023/6/18.
//

#include "SeamCarving.h"

SeamCarving::SeamCarving(Mat &img , Mat &mask) {
    this->img = img;
    cvtColor(this->img , this->img , COLOR_RGB2GRAY);
    this->mask = mask;
    this->dp = vector<vector<double>>(max(img.rows , img.cols) , vector<double>(max(img.rows , img.cols) , 0));
    this->fa = vector<vector<int>>(max(img.rows , img.cols) , vector<int>(max(img.rows , img.cols) , 0));
    this->dis = vector<vector<Point>>(img.rows , vector<Point>(img.cols , Point(0 , 0)));

    add_seam();


}
//使用e1能量函数
void SeamCarving::calc_e1() {
    Mat gx, gy;
    Sobel(this->img, gx, CV_32F, 1, 0);
    Sobel(this->img, gy, CV_32F, 0, 1);
    magnitude(gx, gy, this->emap);  // 只计算合梯度的幅值
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


void SeamCarving::calc_seam(){

    if(side == TOP || side == BOTTOM){
        for(int i = end.first ; i <= end.second ; i++){
            if(i == end.first){
                for(int j = 0; j < img.rows ;j++){
                    if(mask.at<uchar>(j , i) == BG)dp[j][i] = INF;
                    else if(mask.at<uchar>(j , i) == ADT)dp[j][i] = 3000;
                    else {
                        dp[j][i] = emap.at<float>(j, i);
                    }
                }
            }
            else {
                for (int j = 0; j < img.rows; j++) {
                    dp[j][i] = 0;
                    if(mask.at<uchar>(j , i) == BG)dp[j][i] = INF;
                    else if(mask.at<uchar>(j , i) == ADT)dp[j][i] = 3000;
                    else dp[j][i] = emap.at<float>(j, i);
                    if(j == 0){
                        double CU = 0 , CR = 0 ;
                        CU = dp[j][i - 1];
                        CR = dp[j + 1][i - 1] + abs((int)img.at<uchar>(j, i - 1) - (int)img.at<uchar>(j + 1 , i));
                        if(CU < CR){
                            dp[j][i] += CU;
                            fa[j][i] = 0;
                        }
                        else {
                            dp[j][i] += CR;
                            fa[j][i] = 1;
                        }
                    }
                    else if(j == img.rows - 1){
                        double CL = 0 , CU = 0 ;
                        CL = dp[j - 1][i - 1] + abs((int)img.at<uchar>(j , i - 1) - (int)img.at<uchar>(j - 1 , i));
                        CU = dp[j][i - 1];
                        if(CL < CU){
                            dp[j][i] += CL;
                            fa[j][i] = -1;
                        }
                        else {
                            dp[j][i] += CU;
                            fa[j][i] = 0;
                        }
                    }
                    else {
                        double CL = 0, CU = 0, CR = 0;
                        CL += dp[j - 1][i - 1];
                        CL += abs((int)img.at<uchar>(j - 1 , i) - (int)img.at<uchar>(j + 1 , i));
                        CL += abs((int)img.at<uchar>(j , i - 1) - (int)img.at<uchar>(j - 1 , i));
                        CU += dp[j][i - 1];
                        CU += abs((int)img.at<uchar>(j - 1 , i) - (int)img.at<uchar>(j + 1 , i));
                        CR += dp[j + 1][i - 1];
                        CR += abs((int)img.at<uchar>(j - 1 , i) - (int)img.at<uchar>(j + 1 , i));
                        CR += abs((int)img.at<uchar>(j , i - 1) - (int)img.at<uchar>(j + 1 , i));
                        double minn = min({CL , CU , CR});
                        dp[j][i] += minn;
                        if(minn == CL )fa[j][i] = -1;
                        else if(minn == CU) fa[j][i] = 0;
                        else fa[j][i] = 1;
                    }

                }
            }

        }
        int row = 0;
        double minn = 1e18;
        for(int j = 0; j < img.rows ; j++){
            if(dp[j][end.second] < minn){
                minn = dp[j][end.second];
                row = j;
            }
        }
        this->pos.clear();
        int col = end.second;
        bool inf = 0;
        if(dp[row][col] >= INF){
            inf = 1;
            cout<<"check::"<<dp[row][col]<<' ' << col<<'\n';
            cout<<"Error!"<<'\n';
        }
        bool ok = 1;
        while(col > end.first){
            this->pos.emplace_back(row);
            if(mask.at<uchar>(row , col) == BG){
                ok = 0;
            }
            row += fa[row][col];
            col--;
            if(row >= img.rows || row < 0){
                cout<<"dp error!"<<'\n';
                cout<<"col out image!"<<'\n';
                exit(0);
            }
        }
        this->pos.emplace_back(row);

        reverse(this->pos.begin() , this->pos.end());

    }
    else {
        for (int i = end.first; i <= end.second; i++) {
            if (i == end.first) {
                for (int j = 0; j < img.cols; j++) {
                    if(mask.at<uchar>(i , j) == BG)dp[i][j] = INF;
                    else if(mask.at<uchar>(i , j) == ADT)dp[i][j] = 3000;
                    else dp[i][j] = emap.at<float>(i, j);
                }
            } else {
                for (int j = 0; j < img.cols; j++) {
                    dp[i][j] = 0;
                    if(mask.at<uchar>(i , j) == BG)dp[i][j] = INF;
                    else if(mask.at<uchar>(i , j) == ADT)dp[i][j] = 3000;
                    else dp[i][j] = emap.at<float>(i, j);
                    if (j == 0) {
                        double CU = 0, CR = 0;
                        CU = dp[i - 1][j];
                        CR = dp[i - 1][j + 1] +
                             abs((int) img.at<uchar>(i - 1, j) - (int) img.at<uchar>(i, j + 1));

                        if (CU < CR) {
                            dp[i][j] += CU;
                            fa[i][j] = 0;
                        } else {
                            dp[i][j] += CR;
                            fa[i][j] = 1;
                        }
                    } else if (j == img.cols - 1) {
                        double CL = 0, CU = 0;
                        CL = dp[i - 1][j - 1] +
                             abs((int) img.at<uchar>(i - 1, j) - (int) img.at<uchar>(i, j - 1));
                        CU = dp[i - 1][j];
                        if (CL < CU) {
                            dp[i][j] += CL;
                            fa[i][j] = -1;
                        } else {
                            dp[i][j] += CU;
                            fa[i][j] = 0;
                        }
                    } else {
                        double CL = 0, CU = 0, CR = 0;
                        CL += dp[i - 1][j - 1];
                        CL += abs((int) img.at<uchar>(i, j - 1) - (int) img.at<uchar>(i, j + 1));
                        CL += abs((int) img.at<uchar>(i - 1, j) - (int) img.at<uchar>(i, j - 1));
                        CU += dp[i - 1][j];
                        CU += abs((int) img.at<uchar>(i, j - 1) - (int) img.at<uchar>(i, j + 1));
                        CR += dp[i - 1][j + 1];
                        CR += abs((int) img.at<uchar>(i, j - 1) - (int) img.at<uchar>(i, j + 1));
                        CR += abs((int) img.at<uchar>(i - 1, j) - (int) img.at<uchar>(i, j + 1));
                        double minn = min({CL, CU, CR});
                        dp[i][j] += minn;
                        if (minn == CL)fa[i][j] = -1;
                        else if (minn == CU) fa[i][j] = 0;
                        else fa[i][j] = 1;
                    }

                }
            }

        }
        int col = 0;
        double minn = 1e18;
        for (int j = 0; j < img.cols; j++) {
            if (dp[end.second][j] < minn) {
                minn = dp[end.second][j];
                col = j;
            }
        }
        this->pos.clear();
        int row = end.second;
        bool inf = 0;
        if (dp[row][col] >= INF) {
            inf = 1;
            cout << "check::" << dp[row][col] << ' ' << col << '\n';
            cout << "Error!" << '\n';
        }
        bool ok = 1;
        while (row > end.first) {
            this->pos.emplace_back(col);
            if (mask.at<uchar>(row, col) == BG) {
                ok = 0;
            }
            col += fa[row][col];
            row--;
            if (col >= img.cols || col < 0) {
                cout << "dp error!" << '\n';
                cout << "col out image!" << '\n';
                exit(0);
            }
        }
        this->pos.emplace_back(col);

        reverse(this->pos.begin(), this->pos.end());
    }

}

bool SeamCarving::add_seam(){

    while(1) {
        get_long_boundary();
        if (this->end.first == -1 || this->end.second == this->end.first )break;

        calc_e1();
        calc_seam();

        if(side == TOP || side == BOTTOM){
            if(side == TOP){
                for(int i = end.first ; i <= end.second ; i++){
                    int p = i - end.first;
                    for (int j = 0; j < this->pos[p]; j++) {
                        img.at<uchar>(j, i) = img.at<uchar>(j + 1, i);
                        mask.at<uchar>(j, i) = mask.at<uchar>(j + 1, i);
                        dis[j][i].y = dis[j + 1][i].y ;
                        dis[j][i].x = dis[j + 1][i].x - 1;
                    }
                    if (this->pos[p] == 0 || this->pos[p] == img.rows - 1) {
                        if (this->pos[p] == 0) {
                            img.at<uchar>(this->pos[p] , i ) = img.at<uchar>(this->pos[p] + 1 , i);
                        } else {
                            img.at<uchar>(this->pos[p], i) = img.at<uchar>(this->pos[p] - 1, i);
                        }
                    } else {
                        img.at<uchar>(this->pos[p], i) = img.at<uchar>(this->pos[p] - 1, i) / 2 + img.at<uchar>(this->pos[p] + 1, i) / 2 ;

                    }
                    img.at<uchar>(this->pos[p], i) = ((int)img.at<uchar>(this->pos[p], i) + 127) % 255;

                    if(mask.at<uchar>(this->pos[p] , i) != BG) {
                        mask.at<uchar>( this->pos[p] , i ) = ADT;
                    }
                }
            }
            else {
                for(int i = end.first ; i <= end.second ; i++){
                    int p = i - end.first;
                    for (int j = img.rows - 1; j > this->pos[p]; j--) {
                        img.at<uchar>(j, i) = img.at<uchar>(j - 1, i);
                        mask.at<uchar>(j, i) = mask.at<uchar>(j - 1, i);
                        dis[j][i].y = dis[j - 1][i].y ;
                        dis[j][i].x = dis[j - 1][i].x + 1;
                    }
                    if (this->pos[p] == 0 || this->pos[p] == img.rows - 1) {
                        if (this->pos[p] == 0) {
                            img.at<uchar>(this->pos[p] , i ) = img.at<uchar>(this->pos[p] + 1 , i);
                        } else {
                            img.at<uchar>(this->pos[p], i) = img.at<uchar>(this->pos[p] - 1, i);
                        }
                    } else {
                        img.at<uchar>(this->pos[p], i) = img.at<uchar>(this->pos[p] - 1, i) / 2 + img.at<uchar>(this->pos[p] + 1, i) / 2;
                    }
                    img.at<uchar>(this->pos[p], i) = ((int)img.at<uchar>(this->pos[p], i) + 127) % 255;

                    if(mask.at<uchar>(this->pos[p] , i) != BG) {
                        mask.at<uchar>( this->pos[p] , i ) = ADT;
                    }
                }
            }
        }
        else {
            if(side == LEFT){
                for(int i = end.first ; i <= end.second ; i++){
                    int p = i - end.first;
                    for (int j = 0; j < this->pos[p]; j++) {
                        img.at<uchar>(i, j) = img.at<uchar>(i, j + 1);
                        mask.at<uchar>(i, j) = mask.at<uchar>(i, j + 1);
                        //emap.at<float>(i ,j) = emap.at<float >(i , j + 1);

                        //grayimg.at<uchar>(i , j) = grayimg.at<uchar>(i , j + 1);

                        dis[i][j].y = dis[i][j + 1].y - 1;
                        dis[i][j].x = dis[i][j + 1].x;
                    }
                    if (this->pos[p] == 0 || this->pos[p] == img.cols - 1) {
                        if (this->pos[p] == 0) {
                            img.at<uchar>(i, this->pos[p]) = img.at<uchar>(i, this->pos[p] + 1);
                            //grayimg.at<uchar>(i , this->pos[p]) = grayimg.at<uchar>(i , this->pos[p] + 1);
                        } else {
                            img.at<uchar>(i, this->pos[p]) = img.at<uchar>(i, this->pos[p] - 1);
                            //grayimg.at<uchar>(i , this->pos[p]) = grayimg.at<uchar>(i , this->pos[p] - 1);
                        }
                    } else {
                        img.at<uchar>(i, this->pos[p]) = img.at<uchar>(i, this->pos[p] - 1) / 2 + img.at<uchar>(i, this->pos[p] + 1) / 2;
                    }
                    img.at<uchar>(i, this->pos[p]) = ((int)img.at<uchar>(i, this->pos[p]) + 127) % 255;

                    if(mask.at<uchar>(i , this->pos[p]) != BG) {
                        mask.at<uchar>(i, this->pos[p]) = ADT;
                    }

                }
            }
            else {
                for(int i = end.first ; i <= end.second ; i++){
                    int p = i - end.first;
                    for (int j = img.cols - 1; j > this->pos[p]; j--) {
                        img.at<uchar>(i, j) = img.at<uchar>(i , j - 1);
                        mask.at<uchar>(i, j) = mask.at<uchar>(i, j - 1);
                        dis[i][j].y = dis[i][j - 1].y + 1;
                        dis[i][j].x = dis[i][j - 1].x;
                    }

                    if (this->pos[p] == 0 || this->pos[p] == img.cols - 1) {
                        if (this->pos[p] == 0) {
                            img.at<uchar>(i, this->pos[p]) = img.at<uchar>(i, this->pos[p] + 1);
                        } else {
                            img.at<uchar>(i, this->pos[p]) = img.at<uchar>(i, this->pos[p] - 1);
                        }
                    } else {
                        img.at<uchar>(i, this->pos[p]) = img.at<uchar>(i, this->pos[p] - 1) / 2 + img.at<uchar>(i, this->pos[p] + 1) / 2;
                    }
                    img.at<uchar>(i, this->pos[p]) = ((int)img.at<uchar>(i, this->pos[p]) + 127) % 255;

                    if(mask.at<uchar>(i , this->pos[p]) != BG) {
                        mask.at<uchar>(i, this->pos[p]) = ADT;
                    }
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
    int cols = sqrt(400.0 / a);

    int rows = 400 / cols;
    int row_num = img.rows / rows;
    int col_num = img.cols / cols;

    vector<vector<Point>>res;

    for(int i = 0 ;i < img.rows ; i += row_num){
        vector<Point>pre;
        for(int j = 0; j < img.cols ; j += col_num ){
            pre.emplace_back(Point(i , j));
        }
        if(pre.back().y != img.cols - 1 && img.cols - pre.back().y >= 5){
            pre.emplace_back(Point(i , img.cols - 1));
        }
        res.emplace_back(pre);
    }

    if(res.back().back().x != img.rows - 1 && img.rows - res.back().back().x >= 5){
        vector<Point>pre;
        for(int j = 0; j < img.cols ; j += col_num ){
            pre.emplace_back(Point(img.rows - 1 , j));
        }
        if(pre.back().y != img.cols - 1){
            pre.emplace_back(Point(img.rows - 1 , img.cols - 1));
        }
        res.emplace_back(pre);
    }

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
}
vector<vector<Point>> SeamCarving::get_ordinate() {
    return this->ordinate;
}
