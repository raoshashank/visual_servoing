    /* Visual Servoing Demo 1
    *
    * Date: May 20, 2016
    * ----------------------------------- */

    #include <vs_demo_swagat/ur5model.h>
    #include<cmath>
    #include</home/shashank/catkin_ws/src/vs_demo_swagat/gnuplot_ci/include/gnuplot_ci.h>
    using namespace cv;
    using namespace MODUR5;
    using namespace std;
    using namespace gnuplot_ci;
    #include "ros/ros.h"
    #include "object_tracking/coordinates.h"
    #define ball_size = 190 //target ball size in pixels
    #define NF 2

    const double lambda = 1.0;


    //#define CAMERA_TEST

    #define VS_TEST1 0
    object_tracking::coordinates coo;

    //-----------------------------------------------------------------

    void computer_image_jacobian(double u[], double z, cv::Mat &L)
    {

        if(z == 0)
            z = z + 0.001;

        L.at<double>(0,0) = lambda / z;
        L.at<double>(0,1) = 0.0;
        L.at<double>(0,2) = - u[0] / z;
        L.at<double>(0,3) = -u[0]*u[1]/lambda;
        L.at<double>(0,4) = (lambda*lambda+u[0]*u[0])/lambda;
        L.at<double>(0,5) = -u[1];

        L.at<double>(1,0) = 0;
        L.at<double>(1,1) = lambda/z;
        L.at<double>(1,2) = -u[1]/z;
        L.at<double>(1,3) = -(lambda*lambda+u[1]*u[1])/lambda;
        L.at<double>(1,4) = u[0]*u[1]/lambda;
        L.at<double>(1,5) = u[0];

    }

    //-----------------------------------------------------

    void camera_model(double x[], double u[])
    {
        double z = x[2];
        if(z == 0) z = z + 0.001;



        u[0] = (lambda/z) * x[0];
        u[1] = (lambda/z) * x[1];
    }

    //----------------------------------------------

    void callback(const object_tracking::coordinates& c)
    {
        cout<<"x : " <<coo.x<<endl;
        cout<<"y : " <<coo.y<<endl;
        cout<<"radius" << coo.radius<<endl;
        
        coo = c;
    
    }

    int main(int argc,char** argv)
    {
        ros::init(argc,argv,"ibvs");
        ros::NodeHandle nh;
        ros::Subscriber sub = nh.subscribe("this_topic",1000,callback);
    //   ros::Subscriber sub = nh.Subscriber("/this_topic",10,callback)
        // Define an UR5 object
        UR5 robot;

       ros::AsyncSpinner spinner(1);
       spinner.start();
       ros::waitForShutdown();


    #ifdef CAMERA_TEST


        double u1[NF], u2[NF];  // Pixel Centroid of the object being tracked
        double udot[NF];

        //Image Jacobian Matrix
        cv::Mat L = cv::Mat(NF, NW, CV_64F, 0.0);


        double xc1[NC], xc2[NC];

        double theta[NL] = {0, 0, 0, 0, 0, 0};
        double pose[NW];

        robot.set_joint_angles(theta);
        robot.fwd_kin_pose(pose);

        for(int i = 0; i < NC; i++)
        {
            xc1[i] = pose[i];
            xc2[i] = pose[i];
        }

        camera_model(xc1, u1);
        camera_model(xc2, u2);

        cv::Mat Jg = cv::Mat(NW, NL, CV_64F, 0.0);
        cv::Mat ThDot = cv::Mat(NL,1, CV_64F, 0.0);
        cv::Mat UDot = cv::Mat(NF, 1, CV_64F,0.0);

        double theta_dot[NL] = {DEG2RAD(0), DEG2RAD(10), 0, 0, 0, 0};

        for(int i = 0; i < NL; i++)
            ThDot.at<double>(i) = theta_dot[i];

        double Tmax = 2.0;


        ofstream f1("camera1.txt");
        ofstream f2("pose2.txt");
        ofstream f3("angles.txt");
        ofstream f4("camera2.txt");




        for(int i = 0; i < NF; i++)
        {
            f1 << u1[i] << "\t";
            f4 << u2[i] << "\t";
        }
        f1 << endl;
        f4 << endl;


        double t = 0.0;
        do
        {
        while(ros::ok)
        {
            cout<<::coo.x<<endl<<::coo.y<<endl<<::coo.radius<<endl;
            robot.compute_geometric_jacobian(Jg);

            computer_image_jacobian(u2, xc1[2], L);

            UDot = L * Jg * ThDot;

            for(int i = 0; i < NF; i++)
            {
                udot[i] = UDot.at<double>(i);
                u2[i] = u2[i] + 0.5*udot[i]*dt;

                f4 << u2[i] << "\t";
            }
            f4 << endl;

            for(int i = 0; i < NL; i++)
            {
                theta[i] = theta[i] + ThDot.at<double>(i)*dt;

                f3 << theta[i] << "\t";
            }
            f3 << endl;

            robot.set_joint_angles(theta);DEG2RAD(0), DEG2RAD(10), 0, 0, 0, 0};

            robot.fwd_kin_pose(pose);

            for(int i = 0; i < NW; i++)
            {
                if(i < NC)
                    xc1[i] = pose[i];

                f2 << pose[i] << "\t";
            }
            f2 << endl;


            camera_model(xc1, u1);


            for(int i = 0; i < NF; i++)
                f1 << u1[i] << "\t";
            f1 << endl;

            t = t + dt;

        }}while(t <Tmax);
        f1.close();
        f2.close();
        f3.close();DEG2RAD(0), DEG2RAD(10), 0, 0, 0, 0};

        f4.close();
        //------------------------

        GP_handle G1("/usr/bin/", "X (m)", "Y (m)", "Z (m)");
        G1.gnuplot_cmd("set terminal wxt");
        G1.gnuplot_cmd("splot 'pose2.txt' u 1:2:3 w p t 'with theta'");

        GP_handle G2("/usr/bin/", "X (Pixels)", "Y (Pixels)");
        G2.gnuplot_cmd("plot 'camera1.txt' u 1:2 w p ps 3 t 'with FK', 'camera2.txt' u 1:2 w p t 'with L'");
        G2.gnuplot_cmd("set title 'Camera View'");

        getchar();

    #endif

        //---------------------------------------------------------------------

    #ifdef VS_TEST1

        cv::Mat Jg = cv::Mat(NW, NL, CV_64F, 0.0);
        cv::Mat Jgt = cv::Mat(NL, NW, CV_64F, 0.0);
        cv::Mat SC_Matinv = cv::Mat(NL, NF, CV_64F, 0.0);
        cv::Mat C = cv::Mat(NW, NL, CV_64F, 0.0);


        cv::Mat Theta_dot = cv::Mat(NL,1, CV_64F, 0.0);
        cv::Mat Error = cv::Mat(NF, 1, CV_64F, 0.0);
        cv::Mat L = cv::Mat(NF, NW, CV_64F, 0.0);
        cv::Mat H = cv::Mat(NF, NL, CV_64F, 0.0);
        cv::Mat Hinv = cv::Mat(NL, NF, CV_64F, 0.0);
        cv::Mat q0_dot = cv::Mat(NL, 1, CV_64F, 0.0);
        cv::Mat Cj = cv::Mat(NR,NC,CV_64F,0.0);

        cv::Mat JJt = cv::Mat(NW,NL,CV_64F,0.0);

        double theta_max[NL], theta_min[NL];

        robot.get_angle_range(theta_max, theta_min);

        ofstream f1("ctarget.txt");
        ofstream f2("cactual.txt");
        ofstream f3("pose.txt");
        ofstream f4("config.txt");
        ofstream f5("angles.txt");
        ofstream f6("error.txt");
        
        //Target location
        double theta_t[NL] = {DEG2RAD(-20), DEG2RAD(0), DEG2RAD(0), 0, DEG2RAD(0), DEG2RAD(0)};
        double pose_t[NW], x_t[NC], u_t[NF];
        robot.set_joint_angles(theta_t);
        robot.fwd_kin_pose(pose_t);

        for(int i = 0; i < NC; i++)
            x_t[i] = pose_t[i];

        camera_model(x_t, u_t);
        cv::Mat UT = cv::Mat(NF, 1, CV_64F, u_t);

        for(int i = 0; i < NC; i++)
            f1 << x_t[i] << "\t";

        for(int i = 0; i < NF; i++)
            f1 << u_t[i] << "\t";
        f1 << endl;
        f1.close();


        //current location

        double theta[NL] = {0, 0, 0, 0, 0, 0};
        double pose[NW], x[NC], u[NF];
        robot.set_joint_angles(theta);
        robot.fwd_kin_pose(pose);
        robot.joint_position(Cj);

        for(int i = 0; i < NC; i++)
            x[i] = pose[i];

        camera_model(x, u);
        cv::Mat U = cv::Mat(NF, 1, CV_64F, u);

        //Initial robot configuration
        for(int i = 0; i < NR; i++)
        {
            for(int j = 0; j < NC; j++)
                f4 << Cj.at<double>(i,j) << "\t";
            f4 << endl;
        }
        f4 << endl << endl;

        for(double t = 0; t < TMAX; t = t+dt)
        {
            Error = UT - U;

            robot.compute_geometric_jacobian(Jg);

            computer_image_jacobian(u, x[2], L);

            H = L * Jg;

            ///testing//////////////
    //        cv::transpose(Jg,Jgt);

    //        JJt= Jg*Jgt;

    //        //    std::cout<<"Jg"<<std::endl<<Jg<<std::endl<<"Jgt"<<std::endl<<Jgt<<std::endl<<"JJT"<<std::endl<<JJt<<std::endl;

    //        cv::Mat S, V, Vt;
    //        cv::SVD::compute(JJt, S, V, Vt, cv::SVD::FULL_UV);
    //        //        std::cout << "S" << std::endl << S << std::endl << "V" << std::endl << V << "Vt" << std::endl << Vt << std::endl << std::endl;
    //        //        std::cout << "V * Vt" << std::endl << V * Vt << std::endl;
    //        //        std::cout << "S[0]" << std::endl << S.at<double>(0)<< std::endl;
    //        //        double data[6][6] = {{ S.at<double>(0),0,0,0,0,0 },
    //        //                             { 0,S.at<double>(1),0,0,0,0},
    //        //                             { 0,0,S.at<double>(2),0,0,0},
    //        //                             { 0,0,0,S.at<double>(3),0,0},
    //        //                             { 0,0,0,0,S.at<double>(4),0},
    //        //                             { 0,0,0,0,0,S.at<double>(5)}};
    //        double data[6][6] = {{ std::sqrt(S.at<double>(0)+0.0001),0,0,0,0,0},
    //                             { 0,std::sqrt(S.at<double>(1)+0.0001),0,0,0,0},
    //                             { 0,0,std::sqrt(S.at<double>(2)+0.0001),0,0,0},
    //                             { 0,0,0,std::sqrt(S.at<double>(3)+0.0001),0,0},
    //                             { 0,0,0,0,std::sqrt(S.at<double>(4)+0.0001),0},
    //                             { 0,0,0,0,0,std::sqrt(S.at<double>(5)+0.0001)}};
    //        cv::Mat SC_Mat( 6, 6, CV_64F, data);
    //        cv::invert(SC_Mat, SC_Matinv, cv::DECOMP_SVD);

    ////        std::cout<<"S C Matrix inv \n"<<SC_Matinv<<"\n";

    //        C = V*SC_Matinv*Vt;
    //        std::cout<<" C Matrix  \n"<<C<<"\n";


            /////////////////////////////////////////////////////////////////////////////////////

            //Compute Jacobian inverse
            cv::invert(H, Hinv, cv::DECOMP_SVD);



            for(int i = 0; i < NL; ++i)
            {
                //q0_dot.at<double>(i,0) = 2*(theta[i+1]-pref_config[i]) / (NL * pow((theta_max[i] - theta_min[i]),2.0));
                q0_dot.at<double>(i,0) = 2*(theta[i]) / (NL * pow((theta_max[i] - theta_min[i]),2.0));
            }

            Theta_dot = Hinv * Error - (cv::Mat::eye(NL, NL, CV_64F) - Hinv*H) * q0_dot;

            for(int i = 0; i < NL; i++)
            {
                theta[i] = theta[i] + dt * Theta_dot.at<double>(i);

                f5 << theta[i] << "\t";
            }
            f5 << endl;

            robot.set_joint_angles(theta);
            robot.fwd_kin_pose(pose);

            for(int i = 0; i < NW; i++)
            {
                if(i < NC) x[i] = pose[i];

                f3 << pose[i] << "\t";
            }
            f3 << endl;


            camera_model(x, u);

            for(int i = 0; i < NF; i++)
            {
                U.at<double>(i) = u[i];

                f2 << U.at<double>(i) << "\t";
            }
            f2 << endl;

            f6 << sqrt(cv::Mat(Error.t()*Error).at<double>(0,0)/NF) << endl;
        }



        /////////////////////////////////////////////////////checking/////////////////////////

        //float dataH[3][3] = {{ 2,-1,-1 },
        //                     { -1,2,-1},
        //                     { -1,-1,2}};
        //cv::Mat Hi( 3, 3, CV_32F, dataH );

        //cv::Mat S, V, Vt;
        //cv::SVD::compute(Hi, S, V, Vt, cv::SVD::FULL_UV);
        //std::cout << "S" << std::endl << S << std::endl << "V" << std::endl << V << "Vt" << std::endl << Vt << std::endl << std::endl;

        ////cv::Mat Ut;
        ////cv::transpose( Ui, Ut );
        //std::cout << "V * Vt" << std::endl << V * Vt << std::endl;



        //////////////////////////////////////////////////////////////////////////////////////

        //Final robot configuration
        robot.joint_position(Cj);

        for(int i = 0; i < NR; i++)
        {
            for(int j = 0; j < NC; j++)
                f4 << Cj.at<double>(i,j) << "\t";
            f4 << endl;
        }
        f4 << endl << endl;


        f2.close();
        f3.close();
        f4.close();

        //----------------------------------------------------------------

        // Plotting

        // GP_handle G1("/usr/bin/", "Time (second)", "Error (in Pixels)");
        // G1.gnuplot_cmd("set terminal wxt");
        // G1.gnuplot_cmd("plot 'error.txt' w l lw 2");


        // GP_handle G2("/usr/bin/", "X (Pixels)", " Y (Pixels)");
        // G2.gnuplot_cmd("set terminal wxt");
        // G2.gnuplot_cmd("plot 'ctarget.txt' u 4:5 w p ps 3 t 'Target' , 'cactual.txt' u 1:2 w p lt 3 t 'Actual'");

        // GP_handle G3("/usr/bin/", "X (m)", "Y (m)", "Z (m)");
        // G3.gnuplot_cmd("set terminal wxt");
        // G3.gnuplot_cmd("splot 'pose.txt' u 1:2:3 w p t 'End-effector position', 'ctarget.txt' u 1:2:3 w p ps 3 t 'Target'");
        // G3.gnuplot_cmd("replot 'config.txt' index 0 u 1:2:3 w lp lw 2 t 'Initial Config', '' index 1 u 1:2:3 w lp lw 2 t 'Final Config'");

        // cout << "\n Press Return Key to Exit .... " << endl;

        getchar();




    #endif

        return 0;

    }
