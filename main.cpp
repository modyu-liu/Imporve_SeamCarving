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

#define clamp(x,a,b)    (  ((a)<(b)) \
? ((x)<(a))?(a):(((x)>(b))?(b):(x))	\
: ((x)<(b))?(b):(((x)>(a))?(a):(x))	\
)
GLuint texGround;

GLuint matToTexture(cv::Mat mat, GLenum minFilter = GL_LINEAR,
                    GLenum magFilter = GL_LINEAR, GLenum wrapFilter = GL_REPEAT) {
    //cv::flip(mat, mat, 0);
    // Generate a number for our textureID's unique handle
    GLuint textureID;
    glGenTextures(1, &textureID);

    // Bind to our texture handle
    glBindTexture(GL_TEXTURE_2D, textureID);

    // Catch silly-mistake texture interpolation method for magnification
    if (magFilter == GL_LINEAR_MIPMAP_LINEAR ||
        magFilter == GL_LINEAR_MIPMAP_NEAREST ||
        magFilter == GL_NEAREST_MIPMAP_LINEAR ||
        magFilter == GL_NEAREST_MIPMAP_NEAREST)
    {
        //cout << "You can't use MIPMAPs for magnification - setting filter to GL_LINEAR" << endl;
        magFilter = GL_LINEAR;
    }

    // Set texture interpolation methods for minification and magnification
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);

    // Set texture clamping method
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapFilter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapFilter);

    // Set incoming texture format to:
    // GL_BGR       for CV_CAP_OPENNI_BGR_IMAGE,
    // GL_LUMINANCE for CV_CAP_OPENNI_DISPARITY_MAP,
    // Work out other mappings as required ( there's a list in comments in main() )
    GLenum inputColourFormat = GL_BGR_EXT;
    if (mat.channels() == 1)
    {
        inputColourFormat = GL_LUMINANCE;
    }

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

#define PDD pair<double , double>

vector<vector<Point>>co1 , co2;

Mat img ;

void display() {
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBindTexture(GL_TEXTURE_2D, texGround);
    vector<PDD> d_mesh;
    vector<PDD> d_final_mesh;
    int vertexRow =  co1.size();
    int vertexCol = co1[0].size();
    for (int row = 0; row < vertexRow; row ++) {
        for (int col = 0; col < vertexCol; col ++) {
            int index = row * vertexCol + col;
            PDD d_coord = {co2[row][col].x, co2[row][col].y};
            PDD d_localcoord = {co1[row][col].x, co1[row][col].y};
            //cout<<"check::"<<d_coord.first<<' '<<d_coord.second<<' '<<d_localcoord.first<<' '<<d_localcoord.second<<'\n';

            d_coord.first /= img.rows;
            d_coord.second /= img.cols;
            d_coord.first -= 0.5;
            d_coord.second -= 0.5;
            d_coord.first *= 2;
            d_coord.second *= 2;
            d_coord =  {clamp(d_coord.first, -1, 1), clamp(d_coord.second, -1, 1)};
            //cout << coord << " ";

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
            PDD local_left_top = d_mesh[index];
            PDD local_right_top = d_mesh[index + 1];
            PDD local_left_bottom = d_mesh[index + vertexCol];
            PDD local_right_bottom = d_mesh[index + vertexCol + 1];


            PDD global_left_top = d_final_mesh[index];
            PDD global_right_top = d_final_mesh[index + 1];
            PDD global_left_bottom = d_final_mesh[index + vertexCol];
            PDD global_right_bottom = d_final_mesh[index + vertexCol + 1];


            glBegin(GL_QUADS);
            glTexCoord2d(local_right_top.second, local_right_top.first); glVertex3d(global_right_top.second,  -1 * global_right_top.first, 0.0f);
            glTexCoord2d(local_right_bottom.second, local_right_bottom.first); glVertex3d(global_right_bottom.second,  -1 * global_right_bottom.first, 0.0f);
            glTexCoord2d(local_left_bottom.second, local_left_bottom.first);	glVertex3d(global_left_bottom.second,  -1 * global_left_bottom.first, 0.0f);
            glTexCoord2d(local_left_top.second, local_left_top.first); glVertex3d(global_left_top.second,  -1 * global_left_top.first, 0.0f);
            glEnd();

        }
    }

    //GLFWwindow* window;
    glutSwapBuffers();
}
double shrinkImage(Mat src, Mat &dst) {
    int pixel_numbers = src.rows * src.cols;
    double scale = sqrt(pixel_numbers / 1000000.0);
    resize(src, dst, Size(), 1 / scale, 1 / scale);
    return scale;
}

int main(int argc, char* argv[]){


    img = imread("./data/10a_input.jpg");
    Mat out_img;
    double scale = shrinkImage(img , out_img);
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

    /*
    //绘制网格图
    Mat pre = img.clone();
    for(int i = 0; i < co1.size() - 1 ; i ++){
        for(int j = 0 ; j < co1[0].size() - 1 ; j++){
            Point a = co1[i][j];
            Point b = co1[i][j + 1];
            Point c = co1[i + 1][j + 1];
            Point d = co1[i + 1][j];
            cv::line(pre , Point(a.y , a.x) , Point(b.y , b.x) , Scalar(255 , 0 , 0) , 1);
            cv::line(pre , Point(b.y , b.x) , Point(c.y , c.x) , Scalar(255 , 0 , 0) , 1);
            cv::line(pre , Point(c.y , c.x) , Point(d.y , d.x) , Scalar(255 , 0 , 0) , 1);
            cv::line(pre , Point(d.y , d.x) , Point(a.y , a.x) , Scalar(255 , 0 , 0) , 1);

        }
    }
    imshow("img" , pre);
    waitKey(0);
    */
    time1 = clock();
    globalwarp globalwarp(img , seam.get_ordinate());
    time2 = clock();

    cout<<"globalwarp cost time:: "<<double(time2 - time1) / CLOCKS_PER_SEC <<" s"<<'\n';

    co2 = globalwarp.get_ordinate();

    t2 = clock();
    cout<<"total cost time:: "<<(double)(t2 - t1) / CLOCKS_PER_SEC<<" s"<<'\n';


    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(img.cols, img.rows);
    glutCreateWindow("img");
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // 防止图片倾斜
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_2D);
    texGround = matToTexture(img);
    glutDisplayFunc(&display);
    glutMainLoop();

    return 0 ;
}