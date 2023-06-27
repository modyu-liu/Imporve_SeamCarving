#include "globalwarp.h"
#include "SeamCarving.h"
#include <sys/time.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "GL/glut.h"

void fillHole(Mat src , Mat &dst){
    cv::Size m_Size = src.size();
    Mat Temp = Mat::zeros(m_Size.height + 2, m_Size.width + 2, src.type());
    src.copyTo(Temp(cv::Range(1, m_Size.height + 1), cv::Range(1, m_Size.width + 1)));

    Rect ccomp;
    cv::floodFill(Temp, cv::Point(0, 0), cv::Scalar(255) , &ccomp ,cv::Scalar(5));
    Mat cutImg;
    Temp(cv::Range(1, m_Size.height + 1), cv::Range(1, m_Size.width + 1)).copyTo(cutImg);
    dst = src | (~cutImg);
}

Mat clean(Mat src ){
    Mat bw;
    cv::cvtColor(src, src, cv::COLOR_BGR2GRAY);
    uchar thr = 252;
    Mat mask = Mat::zeros(src.size(), CV_8UC1);
    for (int row = 0; row < src.rows; row++) {
        for (int col = 0; col < src.cols; col++) {
            if (src.at<uchar>(row, col)<thr) {
                mask.at<uchar>(row, col) = 255;
            }
        }
    }
    fillHole(mask, bw);
    bw = ~bw;

    Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(5, 5));
    Mat dilate_out;
    cv::dilate(bw, dilate_out, element);
    cv::dilate(dilate_out, dilate_out, element);
    cv::dilate(dilate_out, dilate_out, element);
    Mat erode_out;
    erode(dilate_out, erode_out, element);
    return erode_out;
}

#define clamp(x,a,b) (((a)<(b))? ((x)<(a))?(a):(((x)>(b))?(b):(x)):((x)<(b))?(b):(((x)>(a))?(a):(x)))
GLuint texGround;

// 纹理基础学习 https://www.bilibili.com/video/BV1UW411A7Vh/?spm_id_from=333.337.search-card.all.click&vd_source=f7a4925df10e6949b85d47acd866d063
GLuint matToTexture(cv::Mat mat, GLenum minFilter = GL_LINEAR,
                    GLenum magFilter = GL_LINEAR, GLenum wrapFilter = GL_REPEAT) {
    // Generate a number for our textureID's unique handle
    GLuint textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapFilter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapFilter);
    GLenum inputColourFormat = GL_BGR_EXT;
    // Create the texture
    glTexImage2D(GL_TEXTURE_2D,     // Type of texture
                 0,                 // Pyramid level (for mip-mapping) - 0 is the top level
                 GL_RGB,            // Internal colour format to convert to
                 mat.cols,          // Image width  i.e. 640 for Kinect in standard mode
                 mat.rows,          // Image height i.e. 480 for Kinect in standard mode
                 0,                 // Border width in pixels (can either be 1 or 0)
                 inputColourFormat, // Input image format (i.e. GL_RGB, GL_RGBA, GL_BGR etc.)
                 GL_UNSIGNED_BYTE,  // Image data type
                 mat.ptr());        // The actual image data itself
    return textureID;
}

vector<vector<Point>>co1 , co2;
Mat img ;

//纹理映射 https://zhuanlan.zhihu.com/p/369977849
//thanks https://github.com/guyuchao/rectangle-panoramic-image
void display() {
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBindTexture(GL_TEXTURE_2D, texGround);
    vector<pair<double , double>> d_mesh;
    vector<pair<double , double>> d_final_mesh;
    int vertexRow =  co1.size();
    int vertexCol = co1[0].size();
    for (int row = 0; row < vertexRow; row ++) {
        for (int col = 0; col < vertexCol; col ++) {
            int index = row * vertexCol + col;
            pair<double , double> d_coord = {co2[row][col].x, co2[row][col].y};
            pair<double , double> d_localcoord = {co1[row][col].x, co1[row][col].y};
            d_coord.first /= img.rows;
            d_coord.second /= img.cols;
            d_coord.first -= 0.5;
            d_coord.second -= 0.5;
            d_coord.first *= 2;
            d_coord.second *= 2;
            d_coord =  {clamp(d_coord.first, -1, 1), clamp(d_coord.second, -1, 1)};
            d_localcoord.first /= img.rows;
            d_localcoord.second /= img.cols;
            d_localcoord = {clamp(d_localcoord.first, 0, 1), clamp(d_localcoord.second, 0, 1)};
            d_final_mesh.push_back(d_coord);
            d_mesh.push_back(d_localcoord);
        }
    }
    for (int row = 0; row < vertexRow - 1; row ++) {
        for (int col = 0; col < vertexCol - 1; col ++) {
            int index = row * vertexCol + col;
            pair<double , double> local_left_top = d_mesh[index];
            pair<double , double> local_right_top = d_mesh[index + 1];
            pair<double , double> local_left_bottom = d_mesh[index + vertexCol];
            pair<double , double> local_right_bottom = d_mesh[index + vertexCol + 1];
            pair<double , double> global_left_top = d_final_mesh[index];
            pair<double , double> global_right_top = d_final_mesh[index + 1];
            pair<double , double> global_left_bottom = d_final_mesh[index + vertexCol];
            pair<double , double> global_right_bottom = d_final_mesh[index + vertexCol + 1];
            glBegin(GL_QUADS);
            glTexCoord2d(local_right_top.second, local_right_top.first);
            glVertex3d(global_right_top.second,  -1 * global_right_top.first, 0.0f);
            glTexCoord2d(local_right_bottom.second, local_right_bottom.first);
            glVertex3d(global_right_bottom.second,  -1 * global_right_bottom.first, 0.0f);
            glTexCoord2d(local_left_bottom.second, local_left_bottom.first);
            glVertex3d(global_left_bottom.second,  -1 * global_left_bottom.first, 0.0f);
            glTexCoord2d(local_left_top.second, local_left_top.first);
            glVertex3d(global_left_top.second,  -1 * global_left_top.first, 0.0f);
            glEnd();

        }
    }
    //glutSwapBuffers();

}

void resizeimg(Mat src, Mat &dst) {
    int num = src.rows * src.cols;
    double scale = sqrt(num / 1000000.0);
    resize(src, dst, Size(), 1 / scale, 1 / scale);
}

int main(int argc, char* argv[]){

    cout<<"please chose picture:"<<'\n';
    int num ;
    cin>>num;
    string filename1 = "./data/" + to_string(num) + "_input.jpg";
    string filename2 = "./data/" + to_string(num) + "_result.jpg";
    img = imread(filename1);
    Mat out_img;
    resizeimg(img , out_img);
    img = out_img;
    Mat mask = clean(img);

    for(int i = 0; i < mask.rows ; i++){
        for(int j = 0 ; j < mask.cols ; j++){
            if(mask.at<uchar>(i , j) == 255){
                mask.at<uchar>(i , j) = BG;
            }
            else {
                mask.at<uchar>(i , j) = FG;
            }
        }
    }

    auto prepare = [&](){
        int up = 0 , bottom = 0 , left = 0 , right = 0;
        bool ok = 1;
        while(ok) {
            for (int i = 0; i < img.cols; i++) {
                if (mask.at<uchar>(up , i) != BG) {
                    ok = 0;
                    break;
                }
            }
            if(ok)up++;
        }
        ok = 1;
        while(ok) {
            for (int i = 0; i < img.cols; i++) {
                if (mask.at<uchar>(img.rows - 1 - bottom , i) != BG) {
                    ok = 0;
                    break;
                }
            }
            if(ok)bottom++;
        }
        ok = 1;
        while(ok) {
            for (int i = 0; i < img.rows; i++) {
                if (mask.at<uchar>(i, left) != BG) {
                    ok = 0;
                    break;
                }
            }
            if(ok)left++;
        }
        ok = 1;
        while(ok) {
            for (int i = 0; i < img.rows; i++) {
                if (mask.at<uchar>(i, img.cols - 1 - right) != BG) {
                    ok = 0;
                    break;
                }
            }
            if(ok)right++;
        }
        Rect roi(left , up , img.cols - left - right , img.rows - up - bottom);
        img = img(roi).clone();
        mask = mask(roi).clone();
    };

    prepare();
    clock_t time1 , time2 , t1 , t2;
    Mat input_img = img.clone();
    time1 = clock();
    t1 = clock();

    SeamCarving seam(input_img , mask);
    time2 = clock();
    cout<<"localwarp cost time:: "<<(double)(time2 - time1) / CLOCKS_PER_SEC<<" s"<<'\n';
    co1 = seam.get_ordinate();
    time1 = clock();
    globalwarp globalwarp(img , seam.get_ordinate());
    //globalwarp.show_seg_line();
    time2 = clock();
    cout<<"globalwarp cost time:: "<<double(time2 - time1) / CLOCKS_PER_SEC <<" s"<<'\n';
    co2 = globalwarp.get_ordinate();
    t2 = clock();
    cout<<"total cost time:: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<" s"<<'\n';
    cout<<"check::"<<img.rows<<' '<<img.cols<<'\n';

    Mat author_img = imread(filename2);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(img.cols, img.rows);
    glutCreateWindow("result");
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // 防止图片倾斜
    glPixelStorei(GL_PACK_SWAP_BYTES, 1);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_2D);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    texGround = matToTexture(img);

    display();
    //glutSwapBuffers();


    Mat result_img(img.rows, img.cols, CV_8UC3);
    glReadPixels(0, 0, img.cols, img.rows, GL_BGR , GL_UNSIGNED_BYTE, result_img.data);
    flip(result_img, result_img, 0);

    imwrite("result.jpg" , result_img);
    cv::resize(result_img , result_img , cv::Size(author_img.cols , author_img.rows));
    Mat result ;
    cv::vconcat(author_img , result_img , result);

    imwrite("contrast.jpg" , result);


    return 0 ;
}