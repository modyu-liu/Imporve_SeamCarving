//
// Created by 刘思语 on 2023/6/21.
//

#include "globalwarp.h"
#include <sys/time.h>

globalwarp::globalwarp(Mat img, vector<vector<Point>> mordinate) {
    this->img = img.clone();
    this->mordinate = mordinate;
    this->seg_line = vector<vector<vector<pair<Pointd , Pointd>>>>(mordinate.size() - 1 , vector<vector<pair<Pointd , Pointd>>>(mordinate[0].size() - 1 , vector<pair<Pointd , Pointd>>()));
    this->allocate = vector<vector<vector<int>>>(this->mordinate.size() - 1,  vector<vector<int>>(this->mordinate[0].size() - 1, vector<int>()));
    this->rotate = vector<double>(50 , 0);
    this->seg_line_w = vector<vector<vector<pair<MatrixXd , MatrixXd>>>>(mordinate.size() - 1 , vector<vector<pair<MatrixXd , MatrixXd>>>(mordinate[0].size() - 1 , vector<pair<MatrixXd , MatrixXd>>()));
    this->ver = vector<vector<Point>>(this->mordinate.size(), vector<Point>(this->mordinate[0].size()));

    this->seg_num = 0;
    get_line();
    mseg_line();
    init_rotate();
    start_learn();
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

    MatrixXd Aqt = Aq.transpose();
    MatrixXd res = Aq * (Aqt * Aq).inverse() * Aqt;
    MatrixXd I = MatrixXd::Identity(8 , 8);

    return res - I;

}

//SparseMatrix 用法 https://www.licc.tech/article?id=22
SparseMatrix<double> globalwarp::get_shape_preservation(){
    SparseMatrix<double> A(8 * (this->mordinate.size() - 1) * (this->mordinate[0].size() - 1) , 8 * (this->mordinate.size() - 1) * (this->mordinate[0].size() - 1));
    for(int i = 0; i < this->mordinate.size() - 1; i ++){
        for(int j = 0; j < this->mordinate[i].size() - 1 ; j++){
            auto block = get_single_shape_preservation(i , j);
            for(int n = 0 ; n < 8 ; n++){
                for(int m = 0 ; m < 8 ; m++){
                    A.insert(8 * (i * ( this->mordinate[0].size() - 1 ) + j) + n , 8 * (i * ( this->mordinate[0].size() - 1 ) + j) + m) = block(n , m);
                }
            }
        }
    }
    A.makeCompressed();
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

    //Mat pre = img.clone();

    for(int i = 0; i < n ; i++){
        Pointd p1(line[i * 7 + 1] , line[i * 7 + 0]);
        Pointd p2(line[i * 7 + 3] , line[i * 7 + 2]);
        this->line.emplace_back(make_pair(p1 , p2));
        //cv::line(pre , Point(p1.y , p1.x) , Point(p2.y , p2.x) , Scalar(255 , 0 , 0) , 2);

    }
    //imshow("img" , pre);
    //waitKey(0);


}
bool globalwarp::inmesh(Pointd& s , int x , int y){
    Pointd a = Pointd(this->mordinate[x][y]);
    Pointd b = Pointd(this->mordinate[x][y + 1]);
    Pointd c = Pointd(this->mordinate[x + 1][y + 1]);
    Pointd d = Pointd(this->mordinate[x + 1][y]);

    Pointd ab(b.x - a.x , b.y - a.y);
    Pointd as(s.x - a.x , s.y - a.y);

    double f1 = (ab.x * as.y - ab.y * as.x);

    Pointd bc(c.x - b.x , c.y - b.y);
    Pointd bs(s.x - b.x , s.y - b.y);
    double f2 = (bc.x * bs.y - bc.y * bs.x);

    Pointd cd(d.x - c.x , d.y - c.y);
    Pointd cs(s.x - c.x , s.y - c.y);
    double f3 = (cd.x * cs.y - cd.y * cs.x);

    Pointd da(a.x - d.x , a.y - d.y);
    Pointd ds(s.x - d.x , s.y - d.y);
    double f4 = (da.x * ds.y - da.y * ds.x);
    // in fact 只需要满足全部小于0即可
    if(f1 <= 0 && f2 <= 0 && f3 <= 0 && f4 <= 0){
        return true;
    }
    else {
        return false;
    }
}
bool globalwarp::is_intersection(Pointd& a , Pointd& b , Pointd& c , Pointd& d){


    Pointd ab = (b - a);
    Pointd cd = (d - c);

    Pointd ac = (c - a);
    Pointd ad = (d - a);

    Pointd ca = (a - c);
    Pointd cb = (b - c);

    double res1 = ab.x * ac.y - ab.y * ac.x;
    double res2 = ab.x * ad.y - ab.y * ad.x;
    double res3 = cd.x * ca.y - cd.y * ca.x;
    double res4 = cd.x * cb.y - cd.y * cb.x;

    if(res1 * res2 <= 0 && res3 * res4 <= 0){
        return true;
    }
    else {
        return false;
    }
}

Pointd globalwarp::get_intersection(Pointd& a , Pointd& b , Pointd& c , Pointd& d){

    Pointd A(a.x - c.x , a.y - c.y) ;
    Pointd B(d.x - c.x , d.y - c.y) ;
    Pointd C(b.x - c.x , b.y - c.y) ;
    Pointd D(d.x - b.x , d.y - b.y) ;
    double d1 = abs(A.x * B.y - A.y * B.x);
    double d2 = abs(C.x * D.y - C.y * D.x);
    double t = d1 / (d1 + d2);
    double x = a.x + (double)(b.x - a.x) * t;
    double y = a.y + (double)(b.y - a.y) * t;

    Pointd ans(x , y);
    return ans;

}

Pointd globalwarp::modify_line(Pointd &s, int i, int j) {
    int x = int(s.x);
    int y = int(s.y);
    Pointd p1 = Pointd(x + 1 , y);
    Pointd p2 = Pointd(x , y + 1);
    Pointd p3 = Pointd(x - 1 , y);
    Pointd p4 = Pointd(x , y - 1);
    Pointd p5 = Pointd(x - 1 , y - 1);
    Pointd p6 = Pointd(x + 1 , y + 1);
    Pointd p7 = Pointd(x + 1 , y - 1);
    Pointd p8 = Pointd(x - 1 , y + 1);

    if(inmesh(p1 , i , j)){
        return p1;
    }
    else if(inmesh(p2 , i , j)){
        return p2;
    }
    else if(inmesh(p3 , i , j)){
        return p3;
    }
    else if(inmesh(p4 , i , j)){
        return p4;
    }
    else if(inmesh(p5 , i , j)){
        return p5;
    }
    else if(inmesh(p6 , i , j)){
        return p6;
    }
    else if(inmesh(p7 , i , j)){
        return p7;
    }
    else if(inmesh(p8 , i , j)){
        return p8;
    }
    else {
        cout<<"have some bug!"<<'\n';
        cout<<s.x<<' '<<s.y<<'\n';
        Pointd a = this->mordinate[i][j];
        Pointd b = this->mordinate[i][j + 1];
        Pointd c = this->mordinate[i + 1][j + 1];
        Pointd d = this->mordinate[i + 1][j];
        cout<<a.x<<' '<<a.y<<' '<<b.x<<' '<<b.y<<' '<<c.x<<' '<<c.y<<" "<<d.x<<' '<<d.y<<'\n';
        exit(0);
    }

}

void globalwarp::mseg_line() {
    for(auto it : this->line) {
        for (int i = 0; i < this->mordinate.size() - 1; i++) {
            for (int j = 0; j < this->mordinate[i].size() - 1; j++) {
                vector <Pointd> pre;
                bool ok = 0;
                if (inmesh(it.first, i, j)) {
                    pre.emplace_back(it.first);
                }
                if (inmesh(it.second, i, j)) {
                    pre.emplace_back(it.second);
                }
                if(pre.size() < 2) {
                    Pointd a = Pointd(this->mordinate[i][j]);
                    Pointd b = Pointd(this->mordinate[i][j + 1]);
                    Pointd c = Pointd(this->mordinate[i + 1][j + 1]);
                    Pointd d = Pointd(this->mordinate[i + 1][j]);
                    Pointd sec;
                    if (is_intersection(a, b, it.first, it.second ) ) {
                        sec = get_intersection(a, b, it.first, it.second);
                        pre.emplace_back(sec);
                    }
                    if (is_intersection(b, c, it.first, it.second)) {
                        sec = get_intersection(b, c, it.first, it.second);
                        pre.emplace_back(sec);
                    }
                    if (is_intersection(c, d, it.first, it.second)) {
                        sec = get_intersection(c, d, it.first, it.second);
                        pre.emplace_back(sec);
                    }
                    if (is_intersection(d, a, it.first, it.second)) {
                        sec = get_intersection(d, a, it.first, it.second);
                        pre.emplace_back(sec);
                    }
                    if (pre.size() < 2)continue;
                    if (pre.size() > 2) {
                        vector <Pointd> cnt;
                        for (int i = 0; i < pre.size(); i++) {
                            for (int j = i + 1; j < pre.size(); j++) {
                                auto now = pre[j] - pre[i];
                                if (now.x * now.x + now.y * now.y >= 4) {
                                    cnt = vector<Pointd>{pre[i], pre[j]};
                                }
                            }
                        }
                        pre = cnt;
                    }
                    if (pre.size() < 2)continue;

                }
                Pointd res = pre[1] - pre[0];
                if (res.x * res.x + res.y * res.y < 4) {
                    continue;
                }
                if(!inmesh(pre[0] , i , j)){
                    pre[0] = modify_line(pre[0] , i , j);
                }
                if(!inmesh(pre[1] , i , j)){
                    pre[1] = modify_line(pre[1] , i , j);
                }
                auto w1 = inv_biliner(pre[0] , i , j);
                auto w2 = inv_biliner(pre[1] , i , j);
                this->seg_line_w[i][j].emplace_back(make_pair(w1 , w2));
                this->seg_line[i][j].emplace_back(make_pair(pre[0], pre[1]));
                this->seg_num++;
            }
        }
    }

}
void globalwarp::init_rotate(){
    vector<vector<pair<Pointd , Pointd>>>Bin(50 , vector<pair<Pointd , Pointd>>());
    for(int i = 0; i < this->mordinate.size() - 1 ; i++){
        for(int j = 0 ; j < this->mordinate[0].size() - 1 ; j++){
            vector<int>pre;
            for(int k = 0 ; k < this->seg_line[i][j].size() ; k++){
                Pointd L = seg_line[i][j][k].second - seg_line[i][j][k].first;
                if(L.y < 0){
                    L.x = -L.x;
                    L.y = -L.y;
                }
                double res = -L.x;
                double cnt = (double)L.x * L.x + (double)L.y * L.y;
                cnt = sqrt(cnt);
                double thea = asin(abs(res / cnt));
                if(res < 0){
                    thea = PI / 2.0 - thea;
                }
                else thea = PI / 2.0  + thea;
                double bin = PI / 50.0;
                int belong = thea / bin;
                Bin[belong].emplace_back(seg_line[i][j][k]);
                pre.emplace_back(belong);
            }
            this->allocate[i][j] = pre;
        }
    }
}
bool globalwarp::in_line(Vector2d &p, Vector2d &a, Vector2d &b) {
    Vector2d ab = b - a;
    Vector2d ap = p - a;
    if(abs((ab(0) * ap(1) - ab(1) * ap(0))) <= 1e-8 && ab.norm() > ap.norm())return true;
    else return false;
}

//逆双线性差值推导 https://www.cnblogs.com/lipoicyclic/p/16338901.html
MatrixXd globalwarp::inv_biliner(Pointd& P , int x, int y) {
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
    int way = 0;
    if(in_line( p , a , b)){
        way = 1;
        double u;
        if((b - a)(0) != 0) {
            u = (p - a)(0) / (b - a)(0);
        }
        else {
            u = (p - a)(1) / (b - a)(1);
        }
        if(u > 1){
            cout<<"in this wa"<<'\n';
            cout<<(p - a)(0)<<' '<<(b - a)(0)<<'\n';
            cout<<p(0)<<' '<<p(1)<<' '<<a(0)<<' '<<a(1)<<' '<<b(0)<<' '<<b(1)<<'\n';
            cout<<"check::"<<A.x<<' '<<A.y<<' '<<B.x<<' '<<B.y<<' '<<C.x<<' '<<C.y<<' '<<D.x<<' '<<D.y<<'\n';
            cout<<"po::"<<P.x<<' '<<P.y<<'\n';
            cout<<"error in inv_bi!"<<'\n';
            exit(0);
        }
        w1 = 1 - u;
        w2 = u;
        w3 = 0;
        w4 = 0;
    }
    else if(in_line(p , b , c)){
        way = 2;
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
        way = 3;
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
        way = 4;
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
        way = 5;
        Vector2d e = b - a;
        Vector2d f = d - a;
        Vector2d g = a - b + c - d;
        Vector2d h = p - a;
        auto cross2d = [&](Vector2d v1, Vector2d v2) {
            return v1(0) * v2(1) - v1(1) * v2(0);
        };
        double k2 = cross2d(g, f);
        double k1 = cross2d(e, f) + cross2d(h, g);
        double k0 = cross2d(h, e);
        double u, v;
        int flag = 0 ;

        bool ppp = 0;
        if (abs(k2) < 0.001) {
            flag = 1;
            if((e(0) * k1 - g(0) * k0) == 0){
                u = (h(1) * k1 + f(1) * k0) / (e(1) * k1 - g(1) * k0);
            }
            else {
                u = (h(0) * k1 + f(0) * k0) / (e(0) * k1 - g(0) * k0);
            }
            v = -k0 / k1;
        }
        else {
            flag = 2;
            double w = k1 * k1 - 4.0 * k0 * k2;
            if (w < 0.0) {
                cout << "no solution!" << '\n';
                exit(0);
            }
            w = sqrt(w);
            double ik2 = 0.5 / k2;
            v = (-k1 - w) * ik2;
            if((e(1) + g(1) * v) == 0 && (e(0) + g(0) * v) == 0){
                exit(0);
            }
            double u1 = -1, u2 = -1;
            if((e(1) + g(1) * v) != 0){
                u1 = (h(1) - f(1) * v) / (e(1) + g(1) * v);
            }
            if((e(0) + g(0) * v) != 0) {
                u2 = (h(0) - f(0) * v) / (e(0) + g(0) * v);
            }
            if(u1 >= 0 && u1 <= 1)u = u1;
            else u = u2;

            if (u < 0.0 || u > 1.0 || v < 0.0 || v > 1.0) {
                ppp = 1;
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
        if(u > 1 || v > 1){
            cout<<u<<' '<<v<< ' ' << flag << ' ' << ppp << '\n';
            cout<<"check::"<<A.x<<' '<<A.y<<' '<<B.x<<' '<<B.y<<' '<<C.x<<' '<<C.y<<' '<<D.x<<' '<<D.y<<'\n';
            cout<<"po::"<<P.x<<' '<<P.y<<'\n';
            cout<<"error in inv_bi!"<<'\n';
            exit(0);
        }
        w1 = 1 - u - v + u * v;
        w2 = u - u * v;
        w3 = u * v;
        w4 = v - u * v;

    }
    if(isnan(w1) ||isnan(w2) || isnan(w3) || isnan(w4)  ){
        cout<<"nan in w "<<'\n';
        exit(0);
    }
    if(abs(w1) < 1e-9)w1 = 0;
    if(abs(w2) < 1e-9)w2 = 0;
    if(abs(w3) < 1e-9)w3 = 0;
    if(abs(w4) < 1e-9)w4 = 0;

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
                pair<Pointd , Pointd> it = this->seg_line[i][j][k];
                MatrixXd w1 = this->seg_line_w[i][j][k].first;
                MatrixXd w2 = this->seg_line_w[i][j][k].second;
                int B = this->allocate[i][j][k];
                double theta = this->rotate[B];
                MatrixXd R(2 , 2);
                R << cos(theta) , -sin(theta) , sin(theta) , cos(theta);
                MatrixXd E(2 , 1);
                E << (it.second - it.first).x , (it.second - it.first).y;
                MatrixXd C = R * E * (E.transpose() * E).inverse() * E.transpose() * R.transpose() - MatrixXd::Identity(2 , 2);

                MatrixXd w = w2 - w1;
                MatrixXd Ce = C * w;
                for(int g = 0 ; g < 8 ; g++){
                    line_pre.insert(idx , 8 *  ( i * (this->mordinate[0].size() - 1) + j ) + g) = Ce(0 , g);
                }
                idx++;
                for(int g = 0 ; g < 8 ; g++){
                    line_pre.insert(idx , 8 *  ( i * (this->mordinate[0].size() - 1) + j ) + g) = Ce(1 , g);
                }
                idx++;
            }
        }
    }
    line_pre.makeCompressed();
    return line_pre;

}
pair<SparseMatrix<double> , VectorXd> globalwarp::get_boundary_constraints() {
    int vernum = (this->mordinate.size() ) * (this->mordinate[0].size() );
    SparseMatrix<double> bound(2 * vernum , 2 * vernum);
    VectorXd b(2 * vernum , 1);

    for(int i = 0; i < this->mordinate[0].size() ; i++){
        int idx1 = 2 * i;
        int idx2 = 2 * ((this->mordinate[0].size()) * (this->mordinate.size() - 1) + i);
        bound.insert(idx1 , idx1) = 1;
        bound.insert(idx2 , idx2) = 1;
        b(idx1) = 0;
        b(idx2) = this->img.rows - 1;
    }

    for(int i = 0; i < this->mordinate.size(); i++){

        int idx1 = 2 * i * (this->mordinate[0].size());
        int idx2 = 2 * ( i * (this->mordinate[0].size()) + this->mordinate[0].size() - 1 );

        bound.insert(idx1 + 1, idx1 + 1) = 1;
        bound.insert(idx2 + 1, idx2 + 1) = 1;
        b(idx1 + 1) = 0;
        b(idx2 + 1) = this->img.cols - 1;
    }
    bound.makeCompressed();
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
    pos.makeCompressed();
    return pos;

}
void globalwarp::update_rotate(){

    vector<int>num(50 , 0);
    for(int i = 0; i < 50 ; i++){
        this->rotate[i] = 0 ;
    }
    for(int i = 0; i < this->mordinate.size() - 1; i ++){
        for(int j = 0; j < this->mordinate[0].size() - 1 ; j++){
           Mat pre = img.clone();
           Vector2d a(this->ver[i][j].x , this->ver[i][j].y);
           Vector2d b(this->ver[i][j + 1].x , this->ver[i][j + 1].y);
           Vector2d c(this->ver[i + 1][j + 1].x , this->ver[i + 1][j + 1].y);
           Vector2d d(this->ver[i + 1][j].x , this->ver[i + 1][j].y);
           VectorXd p(8 , 1);
           p << a(0) , a(1) , b(0) , b(1) ,c(0) , c(1) , d(0) , d(1);
           int k = 0;
           for(auto it : this->seg_line[i][j]){
               Vector2d orp = Vector2d((it.second - it.first).x , (it.second - it.first).y);
               MatrixXd w1 = this->seg_line_w[i][j][k].first;
               MatrixXd w2 = this->seg_line_w[i][j][k].second;
               Vector2d p1 = w1 * p;
               Vector2d p2 = w2 * p;
               Vector2d P = p2 - p1;
               double theta1 = atan2(orp(1) , orp(0));
               double theta2 = atan2(P(1) , P(0));
               double theta = theta2 - theta1 ;
               if(theta > PI || theta < -PI){
                   if(theta > PI){
                       theta -= 2 * PI;
                   }
                   else {
                       theta += 2 * PI;
                   }
               }
               int B = this->allocate[i][j][k];
               this->rotate[B] += theta;
               num[B]++;
               k++;
           }

        }
    }
    for(int i = 0; i < 50 ; i++){
        if(num[i] == 0)continue;
        this->rotate[i] /= (double)num[i];
        if(isnan(this->rotate[i]) ) {
            cout<<"rotate is nan!"<<'\n';
        }
    }

}

SparseMatrix<double> globalwarp::Connect_mat(SparseMatrix<double>& m1 , SparseMatrix<double>& m2){

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
    pair<SparseMatrix<double> , VectorXd> Bb = get_boundary_constraints();
    SparseMatrix<double> bound_p = Bb.first * sqrt(lamada2);
    SparseMatrix<double> k1 = Connect_mat(bound_p , shape_p);
    for(int i = 0; i < 10 ; i++) {
        SparseMatrix<double> line_p = lamada1 * get_line_preservation() * pos / ((double)this->seg_num) ;
        VectorXd b = VectorXd::Zero(bound_p.rows() + shape_p.rows() + line_p.rows(), 1);
        b.block(0, 0, Bb.second.rows(), 1) = (Bb.second * sqrt(lamada2));
        SparseMatrix<double> k = Connect_mat(k1, line_p);
        SparseMatrix<double> k_trans = k.transpose();
        SparseMatrix<double> K = k_trans * k;
        auto B = k_trans * b;
        auto *solver = new SimplicialCholesky <SparseMatrix<double>>(K);
        VectorXd V = solver->solve(B);
        for(int i = 0 ; i < this->mordinate.size() ; i++){
            for(int j = 0 ; j < this->mordinate[0].size() ; j++){
                int idx = i * this->mordinate[0].size() + j;
                this->ver[i][j] = Point(V(2 * idx) , V(2 * idx + 1));
            }
        }
        update_rotate();
    }

}

vector<vector<Point>> globalwarp::get_ordinate() {
    return this->ver;
}
void globalwarp::show_line() {
    Mat pre = img.clone();
    for(auto it : this->line){
        cv::line(pre , Point(it.first.y , it.first.x) , Point(it.second.y , it.second.x) , Scalar(255 , 0 , 0) , 1);

    }
    imshow("img" , pre);
    waitKey(0);
}
void globalwarp::show_seg_line() {
    Mat pre = img.clone();
    for(int i = 0; i < this->mordinate.size() - 1 ; i++){
        for(int j = 0; j < this->mordinate[0].size() - 1; j ++){
            for(auto it : this->seg_line[i][j]){
                cv::line(pre , Point(it.first.y , it.first.x) , Point(it.second.y , it.second.x) , Scalar(255 , 0 , 0) , 1);
            }
        }
    }
    imshow("img" , pre);
    waitKey(0);
}