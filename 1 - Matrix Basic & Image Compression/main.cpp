#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/SVD>
#include<iostream>
#include<opencv2/opencv.hpp>
#include<opencv2/core.hpp>
#include<opencv2/core/eigen.hpp>

#define PI 3.14159265

using namespace Eigen;
using namespace std;
using namespace cv;

MatrixXd image_compression(int n, MatrixXd U, VectorXd S, MatrixXd V){
    int rows = U.rows();
    int cols = V.rows();
    MatrixXd img = MatrixXd::Zero(rows, cols);
    for(int i = 0; i < n; ++i){
        img += S[i] * U.col(i) * V.col(i).transpose();
    }
    img *= 255; // to show the color: we convert the range [0,1] back to [0,255]
    return img;
}

Matrix3d rotate_3d(double x, double y, double z){
    Matrix3d rx, ry, rz;
    double cosx = cos(x/180 * PI);
    double sinx = sin(x/180 * PI);
    double cosy = cos(y/180 * PI);
    double siny = sin(y/180 * PI);
    double cosz = cos(z/180 * PI);
    double sinz = sin(z/180 * PI);
    
    rx << 1, 0,    0,
          0, cosx, -sinx,
          0, sinx, cosx;
    
    ry << cosy,  0, siny,
          0,     1, 0,
          -siny, 0, cosy;
    
    rz << cosz, -sinz, 0,
          sinz, cosz,  0,
          0,    0,     1;

    return rx * ry * rz;
}

int main(){

    // QUESTION 1 //
    cout << "#######################" << endl;
    cout << "QUESTION 1:" << endl;
    cout << "Basic vector operations" << endl;
    cout << "#######################" << endl << endl;

    Vector4d v(1, 1.5, 2, 3);
    Vector4d w(0, 1, 2, 4);

    cout << "vector add: " << endl << "v + w =" << endl << v + w << endl;
    cout << "vector inner product: " << endl << "v · w =" << endl << v.dot(w) << endl;

    // for cross product: only apply for 3D vector
    // convert the 4D vector back to 3D
    Vector3d _v(v[0]/v[3], v[1]/v[3], v[2]/v[3]);
    Vector3d _w(w[0]/w[3], w[1]/w[3], w[2]/w[3]);
    cout << "after convertion to 3D:" << endl;
    cout << "v_3 =" << endl << _v << endl;
    cout << "w_3 =" << endl << _w << endl;

    cout << "vector cross product (in 3D): " << endl << "v x w =" << endl << _v.cross(_w) << endl;


    // QUESTION 2 //
    cout << endl << "#######################" << endl;
    cout << "QUESTION 2:" << endl;
    cout << "Basic matrix operations" << endl;
    cout << "#######################" << endl << endl;

    Matrix4d i(4, 4);
    i << 1,  2,  3,  4,
            5,  6,  7,  8,
            9,  10, 11, 12,
            13, 14, 15, 16;

    Matrix4d j(4, 4);
    j << 4,  3,  2,  1,
            8,  7,  6,  5,
            12, 11, 10, 9,
            16, 15, 14, 13;
    
    cout << "matrix add: " << endl << "i + j =" << endl << i + j << endl;
    cout << "matrix multiply: " << endl << "i * j =" << endl << i * j << endl;
    cout << "matrix multiply vector: " << endl << "i * v =" << endl << i * v << endl;


    // QUESTION 3 //
    cout << endl << "#######################" << endl;
    cout << "QUESTION 3:" << endl;
    cout << "SVD decomposition of \"lenna\"" << endl;
    cout << "#######################" << endl << endl;

    Mat img = imread("lenna.png", IMREAD_GRAYSCALE);
    MatrixXd lenna(512, 512);
    cv2eigen(img, lenna);
    lenna /= 255;

    JacobiSVD<MatrixXd> svd(lenna, ComputeThinU | ComputeThinV );  
    MatrixXd V = svd.matrixV(), U = svd.matrixU();  
    VectorXd S = svd.singularValues();

    // cout << "matrix U =" << endl << U << endl;
    // cout << "matrix V =" << endl << V << endl;

    MatrixXd feature_1 = image_compression(1, U, S, V);
    MatrixXd feature_10 = image_compression(10, U, S, V);
    MatrixXd feature_50 = image_compression(50, U, S, V);

    Mat save_1, save_10, save_50;
    eigen2cv(feature_1, save_1);
    eigen2cv(feature_10, save_10);
    eigen2cv(feature_50, save_50);
    // imshow("feature_1.png", save_1);
    // imshow("feature_10.png", save_10);
    // imshow("feature_50.png", save_50);
    imwrite("feature_1.png", save_1);
    imwrite("feature_10.png", save_10);
    imwrite("feature_50.png", save_50);
    cout << "the feature map images are saved in the build folder" << endl;
    cout << "if necessary, please check the images in the build folder" << endl;


    // QUESTION 4 //
    cout << endl << "#######################" << endl;
    cout << "QUESTION 4:" << endl;
    cout << "Basic transformation operations" << endl;
    cout << "#######################" << endl << endl;

    Vector3d p1(1, 2, 3);
    Vector3d p2(4, 5, 6);
    Vector3d l = p1 - p2;

    Matrix3d R = rotate_3d(45, 30, 60);
    p1 = R * l + p2;

    cout << "after rotation, the new point is: " << endl << p1 << endl;
    cout << "Note: PI takes 3.14159265" << endl;

    return 0;
}

