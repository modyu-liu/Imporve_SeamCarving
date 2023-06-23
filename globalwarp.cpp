//
// Created by 刘思语 on 2023/6/21.
//

#include "globalwarp.h"

globalwarp::globalwarp(Mat img, vector<vector<Point>> mordinate) {
    this->img = img.clone();
    this->mordinate = mordinate;
    this->seg_line = vector<vector<vector<pair<Point , Point>>>>(mordinate.size() - 1 , vector<vector<pair<Point , Point>>>(mordinate[0].size() - 1 , vector<pair<Point , Point>>()));
    this->allocate = vector<vector<vector<int>>>(this->mordinate.size() - 1,  vector<vector<int>>(this->mordinate[0].size() - 1, vector<int>()));
    this->rotate = vector<double>(50 , 0);
    this->seg_num = 0;

    get_shape_preservation();
    get_line();
    mseg_line();
    init_rotate();
    get_line_preservation();
    start_learn();
    //show();
    /*
    Point a(0 , 0);
    Point b(4 , 1);
    Point c(0 , 1);
    Point d(4 , 0);
    get_intersection(a , b , c , d);
    */
}

MatrixXd globalwarp::get_single_shape_preservation(int x , int y){
    vector<Point>p(4);
    //cout<<"start"<<'\n';
    p[0] = this->mordinate[x][y];
    p[1] = this->mordinate[x][y + 1];
    p[2] = this->mordinate[x + 1][y + 1];
    p[3] = this->mordinate[x + 1][y];
    MatrixXd Aq(8 , 4);
    for(int i = 0; i < 8 ; i++){
        if(i % 2 == 0){
            Aq(i , 0) = (double)p[i / 2].x;
            Aq(i , 1) = (double)-p[i / 2].y;
            Aq(i , 2) = 1.0;
            Aq(i , 3) = 0.0;
        }
        else {
            Aq(i , 0) = (double)p[i / 2].y;
            Aq(i , 1) = (double)p[i / 2].x;
            Aq(i , 2) = 0.0;
            Aq(i , 3) = 1.0;
        }
    }
    //cout<<Aq<<'\n';

    auto Aqt = Aq.transpose();
    auto res = Aq * (Aqt * Aq).inverse() * Aqt;
    auto I = MatrixXd::Identity(8 , 8);
    //cout<<res - I<<'\n';
    //cout<<I<<'\n';

    return res - I;

}

//SparseMatrix 用法 https://www.licc.tech/article?id=22
SparseMatrix<double> globalwarp::get_shape_preservation(){
    SparseMatrix<double> A(8 * (this->mordinate.size() - 1) * (this->mordinate[0].size() - 1) , 8 * (this->mordinate.size() - 1) * (this->mordinate[0].size() - 1));
    for(int i = 0; i < this->mordinate.size() - 1; i ++){
        for(int j = 0; j < this->mordinate[i].size() - 1 ; j++){
            //cout<<"now is "<< i <<' '<< j<<'\n';
            auto block = get_single_shape_preservation(i , j);
            for(int n = 0 ; n < 8 ; n++){
                for(int m = 0 ; m < 8 ; m++){
                    A.insert(8 * (i * ( this->mordinate[0].size() - 1 ) + j) + n , 8 * (i * ( this->mordinate[0].size() - 1 ) + j) + m) = block(n , m);
                }
            }
        }
    }
    //
    // cout<<"finish"<<'\n';
    return A;

}


void globalwarp::get_line() {

    Mat gray_img;
    cvtColor(img , gray_img , COLOR_BGR2GRAY);
    double * pixel = new double[img.rows * img.cols];
    for(int i = 0; i < img.rows ; i++){
        for(int j = 0; j < img.cols ; j++){
            pixel[i * img.cols + j] = (double)gray_img.at<uchar>(i , j);
        }
    }
    int n ;
    double *line = lsd(&n , pixel , img.cols , img.rows);
    //cout<<"ok??"<<n<<'\n';


    for(int i = 0; i < n ; i++){
        Point p1(line[i * 7 + 1] , line[i * 7 + 0]);
        Point p2(line[i * 7 + 3] , line[i * 7 + 2]);
        this->line.emplace_back(make_pair(p1 , p2));
    //        cv::line(img , Point(p1.y , p1.x) , Point(p2.y , p2.x) , Scalar(0 , 0 , 255) , 1 );


    }
    /*
    Point p1(10 , 50);
    Point p2(10 , 100);
    cv::line(img , p1 , p2 , Scalar(255 , 0 , 0) , 5);
    */

}
bool globalwarp::inmesh(Point s , int x , int y){
    Point a = this->mordinate[x][y];
    Point b = this->mordinate[x][y + 1];
    Point c = this->mordinate[x + 1][y + 1];
    Point d = this->mordinate[x + 1][y];
    //
    Point ab(b.x - a.x , b.y - a.y);
    Point as(s.x - a.x , s.y - a.y);
    int f1 = (ab.x * as.y - ab.y * as.x);

    Point bc(c.x - b.x , c.y - b.y);
    Point bs(s.x - b.x , s.y - b.y);
    int f2 = (bc.x * bs.y - bc.y * bs.x);

    Point cd(d.x - c.x , d.y - c.y);
    Point cs(s.x - c.x , s.y - c.y);
    int f3 = (cd.x * cs.y - cd.y * cs.x);

    Point da(a.x - d.x , a.y - d.y);
    Point ds(s.x - d.x , s.y - d.y);
    int f4 = (da.x * ds.y - da.y * ds.x);
    // in fact 只需要满足全部小于0即可
    if(f1 <= 0 && f2<= 0 && f3 <= 0 && f4 <= 0){
        return true;
    }
    else {
        return false;
    }
}
bool globalwarp::is_intersection(Point a , Point b , Point c , Point d){


    Point ab = (b - a);
    Point cd = (d - c);

    Point ac = (c - a);
    Point ad = (d - a);

    Point ca = (a - c);
    Point cb = (b - c);

    int res1 = ab.x * ac.y - ab.y * ac.x;
    int res2 = ab.x * ad.y - ab.y * ad.x;
    int res3 = cd.x * ca.y - cd.y * ca.x;
    int res4 = cd.x * cb.y - cd.y * cb.x;

    if(res1 * res2 <= 0 && res3 * res4 <= 0){
        return true;
    }
    else {
        return false;
    }
}

Point globalwarp::get_intersection(Point a , Point b , Point c , Point d){

    Point A = (a - c);
    Point B = (d - c);
    Point C = (b - c);
    Point D = (d - b);
    double d1 = abs(A.x * B.y - A.y * B.x);
    double d2 = abs(C.x * D.y - C.y * D.x);
//    cout<<"find::"<<' '<<d1<<' '<<d2<<'\n';

    double t = d1 / (d1 + d2);
//    cout<<"ok?"<<t<<'\n';

    Point ans = a + (b - a) * t;

//    cout<<ans.x<<' '<<ans.y<<'\n';
    return ans;

}
void globalwarp::mseg_line() {
    int num = 0;
    //cout<<"check line size "<<' '<<this->line.size()<<'\n';
    int f1[9] = {0 , 0 , 0 , 1 , 1 , 1 , -1 , -1 , -1};
    int f2[9] = {0 , -1 , 1, 0 , -1 , 1 , 0 , -1 , 1};

    for(auto it : this->line) {
        for (int i = 0; i < this->mordinate.size() - 1; i++) {
            for (int j = 0; j < this->mordinate[i].size() - 1; j++) {
                //for (auto it: this->line) {
                vector <Point> pre;
                //cout<<"start!"<<'\n';

                if (inmesh(it.first, i, j)) {
                  //  cout<<"error 1 !"<<'\n';
                    pre.emplace_back(it.first);
                }
                if (inmesh(it.second, i, j)) {
                   // cout<<"error 2 !"<<'\n';
                    pre.emplace_back(it.second);
                }
                Point a = this->mordinate[i][j];
                Point b = this->mordinate[i][j + 1];
                Point c = this->mordinate[i + 1][j + 1];
                Point d = this->mordinate[i + 1][j];


                if (is_intersection(a, b, it.first, it.second)) {
                    Point sec = get_intersection(a, b, it.first, it.second);
                   // cout<<"error 3 !"<<'\n';
                    pre.emplace_back(sec);
                }
                if (is_intersection(b, c, it.first, it.second)) {
                    Point sec = get_intersection(b, c, it.first, it.second);
                   // cout<<"error 4 !"<<'\n';
                    pre.emplace_back(sec);
                }
                if (is_intersection(c, d, it.first, it.second)) {
                    Point sec = get_intersection(c, d, it.first, it.second);
                  //  cout<<"error 5 !"<<'\n';
                    pre.emplace_back(sec);
                }
                if (is_intersection(d, a, it.first, it.second)) {
                    Point sec = get_intersection(d, a, it.first, it.second);
                  //  cout<<"error 6 !"<<'\n';
                    pre.emplace_back(sec);
                }
                if (pre.size() < 2)continue;
                if (pre.size() > 2) {
                    vector<Point>cnt;
                    for(int i = 0; i < pre.size() ; i++){
                        for(int j = i + 1; j < pre.size() ; j++ ){
                            auto now = pre[j] - pre[i];
                            if(now.x * now.x + now.y * now.y > 2){
                                cnt = vector<Point>{pre[i] , pre[j]};
                            }
                        }
                    }
                    pre = cnt;
                }
                if(pre.size() < 2)continue;

                for(int k = 0; k < 9 ; k++){
                    int x = pre[0].x + f1[k];
                    int y = pre[0].y + f2[k];
                    if(inmesh(Point(x , y) , i , j)){
                        pre[0].x = x;
                        pre[0].y = y;
                        break;
                    }
                }
                for(int k = 0; k < 9 ; k++){
                    int x = pre[1].x + f1[k];
                    int y = pre[1].y + f2[k];
                    if(inmesh(Point(x , y) , i , j)){
                        pre[1].x = x;
                        pre[1].y = y;
                        break;
                    }
                }

                Point res = pre[1] - pre[0];
                if (res.x * res.x + res.y * res.y <= 2) {
                    continue;
                }
                this->seg_line[i][j].emplace_back(make_pair(pre[0], pre[1]));
                this->seg_num++;

                /*
                Mat Pre = img.clone();
                cv::line(Pre, Point(a.y , a.x) , Point(b.y , b.x) , Scalar(255, 0 , 0) , 1);
                cv::line(Pre, Point(b.y , b.x) , Point(c.y , c.x) , Scalar(255, 0 , 0) , 1);
                cv::line(Pre, Point(c.y , c.x) , Point(d.y , d.x) , Scalar(255, 0 , 0) , 1);
                cv::line(Pre, Point(d.y , d.x) , Point(a.y , a.x) , Scalar(255, 0 , 0) , 1);
                cv::line(Pre, Point(pre[0].y , pre[0].x) , Point(pre[1].y , pre[1].x) , Scalar(0, 0 , 255) , 1);
                cout<<"check::"<<i<<' '<<j<<'\n';

                imshow("img" , Pre);
                waitKey(0);
                */
                /*
                circle(img, Point(a.y, a.x), 2, Scalar(0, 255, 255), -1);
                circle(img, Point(b.y, b.x), 2, Scalar(0, 255, 255), -1);
                circle(img, Point(c.y, c.x), 2, Scalar(0, 255, 0), -1);
                circle(img, Point(d.y, d.x), 2, Scalar(0, 255, 0), -1);

                Point now1(pre[0].y, pre[0].x);
                Point now2(pre[1].y, pre[1].x);
                cv::line(img , now1 , now2 , Scalar(255 , 0 , 0) , 1);
                //cv::line(img, Point(it.first.y, it.first.x), Point(it.second.y, it.second.x), Scalar(255, 0, 255), 2);
                cout<<"check::"<<a.x<<' '<<a.y<<' '<<b.x<<' '<<b.y<<' '<<it.first.x<<' '<<it.first.y<<' '<<it.second.x<<' '<<it.second.y<<'\n';
                imshow("img", img);
                waitKey(0);
                cout << "finish" << '\n';
                */



            }
        }
    }
    //cout<<num<<'\n';


}
void globalwarp::init_rotate(){
    vector<vector<pair<Point , Point>>>Bin(50 , vector<pair<Point , Point>>());
    for(int i = 0; i < this->mordinate.size() - 1 ; i++){
        for(int j = 0 ; j < this->mordinate[0].size() - 1 ; j++){
            vector<int>pre;
            for(int k = 0 ; k < this->seg_line[i][j].size() ; k++){
                //cout<<"find::"<<i<<' '<<j<<' '<<k<<'\n';
                Point L = seg_line[i][j][k].second - seg_line[i][j][k].first;
                //cout<<"ok??"<<'\n';
                if(L.y < 0){
                    L.x = -L.x;
                    L.y = -L.y;
                }
                //cout<<"ok2?"<<'\n';
                double res = -L.x;
                double cnt = (double)L.x * L.x + (double)L.y * L.y;
                cnt = sqrt(cnt);
                //cout<<"ok3?"<<'\n';
                double thea = asin(abs(res / cnt));
                //cout<<"ok4?"<<'\n';
                if(res < 0){
                    thea = PI / 2.0 - thea;
                }
                else thea = PI / 2.0  + thea;
                double bin = PI / 50.0;
                int belong = thea / bin;
                //cout<<"check::"<<belong<<'\n';
                Bin[belong].emplace_back(seg_line[i][j][k]);
                pre.emplace_back(belong);
            }
            this->allocate[i][j] = pre;
        }
    }
    /*
    //check
    for(int i = 0; i < 50 ; i++){
        Mat pre = img.clone();
        if(Bin[i].size() == 0)continue;
        for(auto it : Bin[i]){
            cout<<it.first.x<<' '<<it.first.y<<' '<<it.second.x<<' '<<it.second.y<<'\n';
            cv::line(pre , Point(it.first.y , it.first.x) , Point(it.second.y , it.second.x) , Scalar(255 , 0 , 0) , 2);
        }
        cout<<"bin : "<<i<<'\n';
        imshow("img" , pre);
        waitKey(0);
    }
    */
}
bool globalwarp::in_line(Vector2d &p, Vector2d &a, Vector2d &b) {
    auto ab = b - a;
    auto ap = p - a;
    if((ab(0) * ap(1) - ab(1) * ap(0)) == 0)return true;
    else return false;
}
//逆双线性差值推导 https://www.cnblogs.com/lipoicyclic/p/16338901.html
MatrixXd globalwarp::inv_biliner(Point P , int x, int y) {
    Point A = this->mordinate[x][y];
    Point B = this->mordinate[x][y + 1];
    Point C = this->mordinate[x + 1][y + 1];
    Point D = this->mordinate[x + 1][y];
    Vector2d p(P.x , P.y);
    Vector2d a(A.x , A.y);
    Vector2d b(B.x , B.y);
    Vector2d c(C.x , C.y);
    Vector2d d(D.x , D.y);
    double w1 , w2 , w3 , w4;
    if(in_line( p , a , b)){
        double u;
        if((b - a)(0) != 0) {
            u = (p - a)(0) / (b - a)(0);
        }
        else {
            u = (p - a)(1) / (b - a)(1);
        }

        w1 = 1 - u;
        w2 = u;
        w3 = 0;
        w4 = 0;
    }
    else if(in_line(p , b , c)){
        double u;
        if((c - b)(0) != 0) {
            u = (p - b)(0) / (c - b)(0);
        }
        else {
            u = (p - b)(1) / (c - b)(1);
        }
        w1 = 0;
        w2 = 1 - u;
        w3 = u;
        w4 = 0;
    }
    else if(in_line(p , c , d)){
        double u;
        if((d - c)(0) != 0) {
            u = (p - c)(0) / (d - c)(0);
        }
        else {
            u = (p - c)(1) / (d - c)(1);
        }
        w1 = 0;
        w2 = 0;
        w3 = 1 - u;
        w4 = u;
    }
    else if(in_line(p , d , a)){
        double u;

        if((a - d)(0) != 0) {
            u = (p - d)(0) / (a - d)(0);
        }
        else {
            u = (p - d)(1) / (a - d)(1);
        }
        w1 = u;
        w2 = 0;
        w3 = 0;
        w4 = 1 - u;
    }
    else {
        auto e = b - a;
        auto f = d - a;
        auto g = a - b + c - d;
        auto h = p - a;
        auto cross2d = [&](Vector2d v1, Vector2d v2) {
            return v1(0) * v2(1) - v1(1) * v2(0);
        };
        double k2 = cross2d(g, f);
        double k1 = cross2d(e, f) + cross2d(h, g);
        double k0 = cross2d(h, e);
        double u, v;
        int flag = 0;
        if (abs(k2) < 0.001) {

            if((e(1) * k1 - g(1) * k0) == 0 && (e(0) * k1 - g(0) * k0) == 0){
                cout<<"error !"<<'\n';
                exit(0);
            }

            if((e(0) * k1 - g(0) * k0) == 0){
                u = (h(1) * k1 + f(1) * k0) / (e(1) * k1 - g(1) * k0);
            }
            else {
                u = (h(0) * k1 + f(0) * k0) / (e(0) * k1 - g(0) * k0);
            }
            v = -k0 / k1;
        }
        else {
            double w = k1 * k1 - 4.0 * k0 * k2;
            if (w < 0.0) {
                cout << "no solution!" << '\n';
                exit(0);
            }
            w = sqrt(w);
            double ik2 = 0.5 / k2;
            v = (-k1 - w) * ik2;
            if((e(1) + g(1) * v) == 0 && (e(0) + g(0) * v) == 0){
                //cout<<"error ! in this!"<<'\n';
                exit(0);
            }

            if((e(0) + g(0) * v) == 0){
                //cout<<"check::"<<(e(1) + g(1) * v)<<'\n';
                u = (h(1) - f(1) * v) / (e(1) + g(1) * v);
            }
            else {

                //cout<<"check::"<<(e(0) + g(0) * v)<<'\n';
                //cout<<"ok??"<<e(0)<<' '<<g(0)<<' '<<v<<'\n';

                u = (h(0) - f(0) * v) / (e(0) + g(0) * v);
            }
            /*
            if((e(0) + g(0) * v) == 0){
                cout<<e(0) <<' '<<g(0)<<' '<<v<<'\n';
                cout<<e(1) <<' '<<g(1)<<' '<<v<<'\n';

                cout<<"wa wa wa"<<'\n';
                exit(0);

            }
             */
            if (u < 0.0 || u > 1.0 || v < 0.0 || v > 1.0) {
                v = (-k1 + w) * ik2;
                if((e(0) + g(0) * v) == 0 && (e(1) + g(1) * v) == 0){
                    cout<<"error in change !"<<'\n';
                    exit(0);

                }
                if((e(0) + g(0) * v) != 0) {
                    u = (h(0) - f(0) * v) / (e(0) + g(0) * v);
                }
                else {
                    u = (h(1) - f(1) * v) / (e(1) + g(1) * v);

                }
            }
        }
        /*
        auto xx = a + (b - a) * u + (d - a) * v + (a - b + c - d) * u * v;
        double dis = (xx(0) - p(0)) * (xx(0) - p(0)) + (xx(1) - p(1)) * (xx(1) - p(1));
        dis = sqrt(dis);
        if(dis > 0.1){
            if(flag == 1){
                cout<<"wa in pa"<<'\n';
            }
            else {
                cout<<"wa in se"<<'\n';
            }
            cout<<"dis long "<< dis<<'\n';
        }
        */
        w1 = 1 - u - v + u * v;
        w2 = u - u * v;
        w3 = u * v;
        w4 = v - u * v;
        if(isnan(w1) ||isnan(w2) || isnan(w3) || isnan(w4)  ){
            cout<<u<<' '<<v<<'\n';

            cout<<"nan in w "<<'\n';
            exit(0);

        }
    }
    MatrixXd w(2, 8);
    w << w1, 0, w2, 0, w3, 0, w4, 0,
            0, w1, 0, w2, 0, w3, 0, w4;
    return w;
}
SparseMatrix<double> globalwarp::get_line_preservation() {

    SparseMatrix<double> line_pre(2 * this->seg_num , 8 * (this->mordinate.size() - 1) * (this->mordinate[0].size() - 1));
    int idx = 0;

    for(int i = 0; i < this->mordinate.size() - 1 ; i ++){
        for(int j = 0 ; j < this->mordinate[0].size() - 1 ; j++){
            for(int k = 0; k < this->seg_line[i][j].size() ; k++) {
                auto it = this->seg_line[i][j][k];
                auto w1 = inv_biliner(it.first , i , j);
                auto w2 = inv_biliner(it.second, i , j);

                int B = this->allocate[i][j][k];
                double theta = this->rotate[B];
                MatrixXd R(2 , 2);
                R << cos(theta) , -sin(theta) , sin(theta) , cos(theta);
                MatrixXd E(2 , 1);
                E << (it.second - it.first).x , (it.second - it.first).y;
                auto C = R * E * (E.transpose() * E).inverse() * E.transpose() * R.transpose() - MatrixXd::Identity(2 , 2);
                auto w = w2 - w1;
                auto Ce = C * w;
                for(int g = 0 ; g < 8 ; g++){
                    if(isnan(Ce(0 , g))){

                        cout<<"Ce(0 , g) is nan" <<'\n';
                        for(int i = 0; i < 2 ; i++){
                            for(int j = 0 ;j < 8 ; j++){
                                cout<<w1(i , j)<<' ';
                            }
                            cout<<'\n';

                        }
                        for(int i = 0; i < 2 ; i++){
                            for(int j = 0 ;j < 8 ; j++){
                                cout<<w2(i , j)<<' ';
                            }
                            cout<<'\n';

                        }
                        Point a = this->mordinate[i][j];
                        Point b = this->mordinate[i][j + 1];
                        Point c = this->mordinate[i + 1][j + 1];
                        Point d = this->mordinate[i + 1][j];

                        cout<<a.x<<' '<<a.y<<' '<<b.x<<' '<<b.y<<' '<<c.x<<" "<<c.y<<' '<<d.x<<' '<<d.y<<'\n';
                        cout<<"point1:" << ' '<<it.first.x <<' '<<it.first.y<<'\n';
                        cout<<"point2:" << ' '<<it.second.x <<' '<<it.second.y<<'\n';

                        exit(0);
                    }
                    line_pre.insert(idx , 8 *  ( i * (this->mordinate[0].size() - 1) + j ) + g) = Ce(0 , g);
                }
                idx++;
                for(int g = 0 ; g < 8 ; g++){
                    if(isnan(Ce(0 , g))){
                        cout<<"Ce(1 , g) is nan" <<'\n';
                        exit(0);
                    }
                    line_pre.insert(idx , 8 *  ( i * (this->mordinate[0].size() - 1) + j ) + g) = Ce(1 , g);
                }
                idx++;

            }

                /*
                // check
                MatrixXd p(8 , 1);
                Point a = this->mordinate[i][j];
                Point b = this->mordinate[i][j + 1];
                Point c = this->mordinate[i + 1][j + 1];
                Point d = this->mordinate[i + 1][j];
                p << a.x , a.y , b.x , b.y , c.x , c.y , d.x , d.y;
                auto ans1 = w1 * p;
                auto ans2 = w2 * p;
                double dis1 = (ans1(0) - it.first.x) * (ans1(0) - it.first.x) + (ans1(1) - it.first.y) * (ans1(1) - it.first.y);
                double dis2 = (ans2(0) - it.second.x) * (ans2(0) - it.second.x) + (ans2(1) - it.second.y) * (ans2(1) - it.second.y);
                dis1 = sqrt(dis1);
                dis2 = sqrt(dis2);
                */

        }
    }
    return line_pre;

}
pair<SparseMatrix<double> , VectorXd> globalwarp::get_boundary_constraints() {
    int vernum = (this->mordinate.size() ) * (this->mordinate[0].size() );
    SparseMatrix<double> bound(2 * vernum , 2 * vernum);
    VectorXd b(2 * vernum , 1);
    // top and bottom
    //map<int , int>ma;

    for(int i = 0; i < this->mordinate[0].size() ; i++){
        int idx1 = 2 * i;
        int idx2 = 2 * ((this->mordinate[0].size()) * (this->mordinate.size() - 1) + i);
        // ma[idx1] = 1;
        // ma[idx2] = 1;
        bound.insert(idx1 , idx1) = 1;
        bound.insert(idx2 , idx2) = 1;
        b(idx1) = 0;
        b(idx2) = this->img.rows - 1;
    }
    //left and right
    //cout<<"check::left ans right"<<'\n';

    for(int i = 0; i < this->mordinate.size(); i++){

        int idx1 = 2 * i * (this->mordinate[0].size());
        int idx2 = 2 * ( i * (this->mordinate[0].size()) + this->mordinate[0].size() - 1 );

        bound.insert(idx1 + 1, idx1 + 1) = 1;
        bound.insert(idx2 + 1, idx2 + 1) = 1;
        b(idx1 + 1) = 0;
        b(idx2 + 1) = this->img.cols - 1;
    }

    return make_pair(bound , b);

}
SparseMatrix<double> globalwarp::get_position_information() {
    int quad_num = (this->mordinate.size() - 1) * (this->mordinate[0].size() - 1);
    int ver_num = this->mordinate.size() * this->mordinate[0].size();

    SparseMatrix<double> pos(8 * quad_num , ver_num * 2);

    for(int i = 0; i < this->mordinate.size() - 1 ; i++){
        for(int j = 0; j < this->mordinate[i].size() - 1 ; j++){

            int idx1 = i * this->mordinate[0].size() + j;
            int idx2 = idx1 + 1;
            int idx3 = idx2 + this->mordinate[0].size();
            int idx4 = idx1 + this->mordinate[0].size();
            int q_num = 8 * ( i * (this->mordinate[0].size() - 1) + j );
            pos.insert(q_num + 0 , 2 * idx1) = 1;
            pos.insert(q_num + 1 , 2 * idx1 + 1) = 1;
            pos.insert(q_num + 2 , 2 * idx2) = 1;
            pos.insert(q_num + 3 , 2 * idx2 + 1) = 1;
            pos.insert(q_num + 4 , 2 * idx3) = 1;
            pos.insert(q_num + 5 , 2 * idx3 + 1) = 1;
            pos.insert(q_num + 6 , 2 * idx4) = 1;
            pos.insert(q_num + 7 , 2 * idx4 + 1) = 1;
        }
    }
    return pos;

}
void globalwarp::update_rotate(VectorXd &V){
    vector<double>r(50);
    vector<int>num(50 , 0);

    for(int i = 0; i < this->mordinate.size() - 1; i ++){
        for(int j = 0; j < this->mordinate[0].size() - 1 ; j++){
            int idx1 = (i * this->mordinate[0].size()) + j;
            int idx2 = idx1 + 1;
            int idx3 = idx2 + this->mordinate[0].size();
            int idx4 = idx1 + this->mordinate[0].size();

            Vector2d a(V(idx1 * 2) , V(idx1 * 2 + 1));
            Vector2d b(V(idx2 * 2) , V(idx2 * 2 + 1));
            Vector2d c(V(idx3 * 2) , V(idx3 * 2 + 1));
            Vector2d d(V(idx4 * 2) , V(idx4 * 2 + 1));

            VectorXd p(8 , 1);
            p << a(0) , a(1) , b(0) , b(1) ,c(0) , c(1) , d(0) , d(1) ;


           int k = 0;
            for(auto it : this->seg_line[i][j]){
                Vector2d orp = Vector2d((it.second - it.first).x , (it.second - it.first).y);
                auto w1 = inv_biliner(it.first , i , j);
                auto w2 = inv_biliner(it.second , i , j);
                Vector2d p1 = w1 * p;
                Vector2d p2 = w2 * p;
                auto P = p2 - p1;
                double len1 = orp.transpose() * orp;
                double len2 = P.transpose() * P;
                //double cnt = (double)dot(p , orp) / (len1 * len2);
                //cout<<cnt<<'\n';

                double cnt = (P(0) * orp(0) + P(1) * orp(1)) / (len1 + len2);
                double theta = acos(cnt);
                cout<<"check::rotate " << ' ' << P(0) <<' '<< P(1)<<' '<<orp(0)<<' '<<orp(1)<< ' ' << theta / PI * 180 << '\n';
                int B = this->allocate[i][j][k];
                if((orp(0) * P(1) - orp(1) * P(0)) < 0){
                    r[B] += theta;
                }
                else {
                    r[B] -= theta;
                }
                num[B]++;

                k++;
            }
        }
    }
    for(int i = 0; i < 50 ; i++){
        rotate[i] /= (double)num[i];
    }
    this->rotate = rotate;
}

SparseMatrix<double> globalwarp::Connect_mat(SparseMatrix<double> m1 , SparseMatrix<double> m2){
    //cout<<m1.rows() <<' '<<m1.cols()<<' '<<m2.rows() <<' '<<m2.cols()<<'\n';

    SparseMatrix<double> m(m1.rows() + m2.rows()  , m1.cols());
    for (int k = 0; k < m1.outerSize(); ++k){
        for (SparseMatrix<double>::InnerIterator it(m1, k); it; ++it){
            m.insert(it.row(), it.col()) = it.value();
        }
    }
    for (int k = 0; k < m2.outerSize(); ++k){
        for (SparseMatrix<double>::InnerIterator it(m2, k); it; ++it){
            m.insert(it.row() + m1.rows(), it.col()) = it.value();
        }
    }
    m.makeCompressed();
    return m;

}

void globalwarp::start_learn() {
    const double lamada1 = 100;
    const double lamada2 = 1e8;
    double quad_num = (this->mordinate.size() - 1) * (this->mordinate[0].size() - 1);

    SparseMatrix<double> pos = get_position_information();

    SparseMatrix<double> shape_p = get_shape_preservation() * pos / (quad_num);

    auto Bb = get_boundary_constraints();
    SparseMatrix<double> bound_p = Bb.first * sqrt(lamada2);
    SparseMatrix<double> line_p = lamada1 * get_line_preservation() * pos / ((double)this->seg_num) ;


    cout<<"check::bound" << ' '<<bound_p.rows()<<' '<<bound_p.cols()<<'\n';
    cout<<"check::shape "<<' '<<shape_p.rows() <<' '<<shape_p.cols()<<'\n';
    cout<<"check::line" << ' '<<line_p.rows()<<' '<<line_p.cols()<<'\n';


    VectorXd b = VectorXd::Zero(bound_p.rows() + shape_p.rows() + line_p.rows() , 1);
    b.block(0 , 0 , Bb.second.rows() , 1) = (Bb.second * sqrt(lamada2)) ;

    cout<<"check?"<<'\n';

    SparseMatrix<double>k1 = Connect_mat(bound_p , shape_p);
    SparseMatrix<double>k2 = Connect_mat(k1 , line_p);
    cout<<"ok?"<<'\n';

    auto k = k2;
    for (int g = 0; g < k.outerSize(); ++g){
        for (SparseMatrix<double>::InnerIterator it(k, g); it; ++it){
            if(isnan(it.value())){
                cout<<"value is nan"<<'\n';
                cout<<it.row()<<'\n';

                exit(0);

            }
          //  cout<<"check:val "<< it.value()<<'\n';

        }
    }
    auto K = k.transpose() * k;

    /*
    cout<<"check::k transpose "<<k.transpose().rows()<<' '<<k.transpose().cols()<<'\n';
    cout<<"check::b "<<' '<<b.rows()<<' '<<b.cols()<<'\n';
    */

    auto B = k.transpose() * b;

    auto *solver = new SimplicialCholesky<SparseMatrix<double>>(K);
    VectorXd V = solver->solve(B);

    auto ans = k * V;
    cout<<"ans:: "<<ans.rows() << ' ' << ans.cols() << ' ' << ans(0)<<'\n';

    cout<<"finish?"<<'\n';

    Mat pre = img.clone();

    for(int i = 0 ; i < V.rows(); i += 2){
        Point now((int)V(i + 1) , (int)V(i));
        circle(pre , now , 2 , Scalar(0 , 0 ,255 ) , -1);

    }
    imshow("img" , pre);
    waitKey(0);
    //update_rotate(V);
    //cout<<b.rows() <<' '<<b.cols()<<'\n';


}

void globalwarp::show(){

    for(int i = 0; i < this->mordinate.size() - 1 ; i++){
        for(int j = 0; j < this->mordinate[0].size() - 1; j ++){
            //cv::circle(img , this->mordinate[i][j] , 2 , Scalar(0 , 255 , 0 ) , 1);

            for(auto it : this->seg_line[i][j]){
                Point a(it.first.y , it.first.x);
                Point b(it.second.y , it.second.x);

                cv::line(img , a , b , Scalar(255 , 0 , 0 ) , 1);
            }
        }
    }
    cout<<"show img"<<'\n';

    imshow("img" , img);
    waitKey(0);

}

