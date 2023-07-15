#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <cmath>
// add some other header files you need

constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 
                 0, 1, 0, -eye_pos[1], 
                 0, 0, 1, -eye_pos[2],
                 0, 0, 0, 1;

    view = translate * view;
    // std::clog << "view" << std::endl << view << std::endl;  // check data

    return view;
}


Eigen::Matrix4f get_model_matrix(float rotation_angle, Eigen::Vector3f T, Eigen::Vector3f S, Eigen::Vector3f P0, Eigen::Vector3f P1)
{

    //Step 1: Build the Translation Matrix M_trans:
    Eigen::Matrix4f Translation;
    Translation <<  1, 0, 0, -T[0], 
                    0, 1, 0, -T[1], 
                    0, 0, 1, -T[2],
                    0, 0, 0, 1;

    //Step 2: Build the Scale Matrix S_trans:
    Eigen::Matrix4f Scale;
    Scale << S[0], 0, 0, 0, 
             0, S[1], 0, 0, 
             0, 0, S[2], 0,
             0, 0,    0, 1;


    //Step 3: Implement Rodrigues' Rotation Formular, rotation by angle theta around axis u, then get the model matrix
	// The axis u is determined by two points, u = P1-P0: Eigen::Vector3f P0 ,Eigen::Vector3f P1  
    // Create the model matrix for rotating the triangle around a given axis. // Hint: normalize axis first
    Eigen::Vector3f u = P1 - P0;
    // normalize to length 1
    double r = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
    u /= r;
    Eigen::Matrix3f I = Eigen::Matrix3f::Identity();
    Eigen::Matrix3f U;
    U << 0, -u[2], u[1],
         u[2], 0, -u[0],
         -u[1], u[0], 0;
    // Implement Rodrigues' Rotation Formular
    double c = cos(rotation_angle * MY_PI/180);
    double s = sin(rotation_angle * MY_PI/180);
    Eigen::Matrix3f R = c * I + s * U + (1 - c) * u * u.transpose();
    // homogeneous coord
    Eigen::Matrix4f R_h = Eigen::Matrix4f::Identity();
    R_h.block<3,3>(0,0) = R;

	//Step 4: Use Eigen's "AngleAxisf" to verify your Rotation
	//Eigen::AngleAxisf rotation_vector(radian, Vector3f(axis[0], axis[1], axis[2]));  
	//Eigen::Matrix3f rotation_matrix;
	//rotation_m = rotation_vector.toRotationMatrix();

    Eigen::Matrix4f model = Translation * R_h * Scale;
	return model;
}



Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    // Implement this function

    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.

    Eigen::Matrix4f Persp, Ortho;

    // frustum -> cubic
    Persp << zNear, 0, 0, 0,
             0, zNear, 0, 0,
             0, 0, zNear+zFar, -zNear*zFar,
             0, 0, 1, 0;

    // orthographic projection
    float top = tan(eye_fov/2 * MY_PI/180) * abs(zNear);
    float right = aspect_ratio * top;
    Ortho << 1/right, 0, 0, 0,
             0, 1/top, 0, 0,
             0, 0, 2/(zFar-zNear), -(zNear+zFar)/2,
             0, 0, 0, 1;

    // squash all transformations
    Eigen::Matrix4f projection = Ortho * Persp;

    // std::clog << "projection" << std::endl << projection << std::endl; //check

    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "result.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
        else
            return 0;
    }

    rst::rasterizer r(1024, 1024);
    // define your eye position "eye_pos" to a proper position
    Eigen::Vector3f eye_pos(0.0, 0.0, 10.0);

    // define a triangle named by "pos" and "ind"
    std::vector<Eigen::Vector3f> pos{{3, 0, -3}, {0, 3, -3}, {-3, 0, -3}};
    // std::vector<Eigen::Vector3f> pos{{2, 0, 0}, {0, 2, 0}, {0, 0, 0}};
    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    // added parameters for get_model_matrix(float rotation_angle, Eigen::Vector3f T, Eigen::Vector3f S, Eigen::Vector3f P0, Eigen::Vector3f P1)
    Eigen::Vector3f T(0.0, 0.0, 0.0);
    Eigen::Vector3f S(1.0, 1.0, 1.0);
    Eigen::Vector3f P0(0, 0, 0);
    Eigen::Vector3f P1(0, 1, 1);
    // Eigen::Vector3f P1(0, 0, 1);

    // added parameters for get_projection_matrix(float eye_fov, float aspect_ratio,float zNear, float zFar)
    float eye_fov = 45;
    float aspect_ratio = 1;
    float zNear = 1;
    float zFar = 10;

    Eigen::Vector3f axis(0, 0, 1);
    // Eigen::Vector3f axis(1, 0, 0);

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle, T, S, P0, P1));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(eye_fov, aspect_ratio, zNear, zFar));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(1024, 1024, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle, T, S, P0, P1));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(eye_fov, aspect_ratio, zNear, zFar));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(1024, 1024, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';
        std::clog << "angle: " << angle << std::endl;
    

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
    }

    return 0;
}
