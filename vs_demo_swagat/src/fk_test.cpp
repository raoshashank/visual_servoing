/* TESTING of GWAM CLASS
 * Date: March 22, 2016
 *
 * Both types of Jacobian are tested
 *
 * Orientation Jacobian is tested through Simulation
 * Angular Velocity Jacobian is tested by Matching it with GWAM manipulator jacobian obtained from the system.
 *
 * ----------------------------------- */

#include <iostream>
#include <ur5model.h>
#include<gnuplot_ci.h>

using namespace std;
using namespace gnuplot_ci;
using namespace MODUR5;

#define JG_TEST // Testing the geometric Jacobian

//#define IK_TEST_J1

//#define CONFIG_TEST


int main()
{
    UR5 robot;

#ifdef CONFIG_TEST

    //double theta[NL] = {DEG2RAD(10), DEG2RAD(30), DEG2RAD(30), 0, DEG2RAD(30), DEG2RAD(20)};
    double theta[NL] = {0,0,0,0,0,0};
    cv::Mat X(NR,NC, CV_64F, 0.0);
    double xef[NC], eepose[NW];

    robot.set_joint_angles(theta);
    robot.joint_position(X);
    robot.fwd_kin_pose(eepose);

    ofstream f1("rconfig.txt");
    ofstream f2("eepose.txt");

    for(int i = 0; i < NR; i++)
    {
        for(int j = 0; j < 3; j++)
            f1 << X.at<double>(i,j) << "\t";
        f1 << endl;
    }

    for(int i = 0; i < NW; i++)
    {
        if(i < NC)
            xef[i] = eepose[i];

        f2 << eepose[i] << "\t";
        cout << eepose[i] << "\t";
    }
    f2 << endl;
    cout << endl;

    f1.close();
    f2.close();

    // Coordinate axes
    double xb[3] = {0,0,0};
    double dx[3] = {0.2, 0.2, 0.2};

    cv::Mat xeef = cv::Mat(3,3, CV_64F, 0.0);
    robot.rotated_EE_axis_with_rot_matrix(dx, xeef);
    //robot.rotated_EE_axis_with_rpy(dx, xeef);


    cout << "xeef = " << xeef << endl;

    double xe[3][3];
    for(int i = 0; i < NC; i++)
        for(int j = 0; j < NC; j++)
            xe[i][j] = xeef.at<double>(j,i);




    GP_handle G4("/usr/bin/", "X (m)", "Y (m)", "Z(m)");
    G4.gnuplot_cmd("set terminal wxt");
    G4.gnuplot_cmd("set border");
    G4.gnuplot_cmd("set ticslevel 0");
    G4.gnuplot_cmd("set yrange [-0.5:0.5]");
    G4.gnuplot_cmd("set zrange [-0.1:0.5]");
  //  G4.gnuplot_cmd("set xrange [-0.5:0.5]");
    G4.gnuplot_cmd("splot 'rconfig.txt' u 1:2:3 w lp lw 2");
    G4.gnuplot_cmd("replot 'eepose.txt' u 1:2:3 w p ps 2");
    G4.draw3dcoordAxis(xb, dx,true);
    G4.draw3dcoordAxis(xef, xe[0], xe[1], xe[2], true,2);


    cout << "Press Enter to exit .." << endl;
    getchar();

#endif


#ifdef JACOB_TEST
    //double theta[NL] = {0, 0, DEG2RAD(30), DEG2RAD(60), 0, DEG2RAD(30), 0};
    double theta[NL] = {DEG2RAD(10), DEG2RAD(30), DEG2RAD(30), 0, DEG2RAD(30), DEG2RAD(20), DEG2RAD(-10)};
    //double theta[NL] = {0, DEG2RAD(0), DEG2RAD(0), 0, DEG2RAD(0), 0, 0};

    robot.set_joint_angles(theta);
    robot.position_jacobian();
    robot.orientation_jacobian();
    robot.linear_velocity_jacobian();
    robot.angular_velocity_jacobian();
    robot.display_variables();


#endif

#ifdef IK_TEST_J1

    int num = 2000;
    double angle[num][NL];

    double init_config[NL] = {0, 0, 0, 0, 0, 0};
    //double pref_config[NL] = {DEG2RAD(-30), DEG2RAD(20), DEG2RAD(-80), 0, 0, 0}; // minimize angle norm
    double pref_config[NL] = {0,0,0, 0, 0, 0}; // minimize angle norm

    double theta[NL];


    //Target pose
    double theta_t[NL] = { DEG2RAD(30), DEG2RAD(-20), DEG2RAD(40), DEG2RAD(10), DEG2RAD(-20), DEG2RAD(30)};
    double pose_t[NW];
    robot.set_joint_angles(theta_t);
    robot.fwd_kin_pose(pose_t);

    ofstream f1("target_pose.txt");
    ofstream f2("rconfig.txt");
    ofstream f3("actual_pose.txt");

    double err = robot.ik_traj_withpose(init_config, pose_t, pref_config, angle, num);


    cout << "error = " <<  err << endl;

    for(int i = 0; i < NW; ++i)
        f1 << pose_t[i] << "\t";
    f1 << endl;

    double pose[NW];
    cv::Mat Jpos(NR,3,CV_64F,0.0);
    for(int cnt = 0; cnt < num; ++cnt)
    {
        for(int i = 0; i < NL; ++i)
        {
            theta[i] = angle[cnt][i];
        }
        robot.set_joint_angles(theta);
        robot.joint_position(Jpos);
        robot.fwd_kin_pose(pose);

        for(int i = 0; i < NR; ++i)
        {
            for(int j = 0; j < 3; ++j)
                f2 << Jpos.at<double>(i,j) << "\t";
            f2 << endl;
        }

        f2 << endl << endl;

        for(int i = 0; i < NW; ++i)
            f3 << pose[i] << "\t";
        f3 << endl;
    }



    f1.close();
    f2.close();
    f3.close();


    // Base coordinate frame
    double xb[3] = {0,0,0};
    double dx[3] = {0.2, 0.2, 0.2};

    // End-effector coordinate frame
    cv::Mat xef = cv::Mat(3,3, CV_64F, 0.0);  //Actual orientation
    robot.rotated_EE_axis_with_rpy(dx, xef);

    cv::Mat xed = cv::Mat(3,3, CV_64F, 0.0); // Desired Orientation

    double rpyd[3];
    for(int i = 0; i < NC; i++)
        rpyd[i] = pose_t[i+NC];

    robot.rotated_EE_axis_with_rpy(dx, rpyd, xed);


    // End-effector position
    double xp[NC];
    for(int i = 0; i < NC; i++)
        xp[i] = pose[i];

    double xe1[3][3], xe2[3][3];
    for(int i = 0; i < NC; i++)
        for(int j = 0; j < NC; j++)
        {
            xe1[i][j] = xef.at<double>(j,i);
            xe2[i][j] = xed.at<double>(j,i);
        }



    cout << "\nJoint Angle Computed (in Degrees) = " ;
    for(int i = 0; i < NL; i++)
        cout << RAD2DEG(theta[i]) << "\t";
    cout << endl;



    GP_handle G4("/usr/bin/", "X (m)", "Y (m)", "Z(m)");
    G4.gnuplot_cmd("set terminal wxt");
    G4.gnuplot_cmd("set border");
    G4.gnuplot_cmd("set ticslevel 0");
    G4.gnuplot_cmd("splot 'target_pose.txt' u 1:2:3 w p ps 2 pt 4 t 'target','actual_pose.txt' u 1:2:3 w p t 'actual'");
    G4.gnuplot_cmd("replot 'rconfig.txt' index 0 u 1:2:3 w lp lw 2 t 'First Config'");
    G4.gnuplot_cmd("replot 'rconfig.txt' index 1999 u 1:2:3 w lp lw 2 t 'Last config'");
    G4.draw3dcoordAxis(xb, dx, true);
    G4.draw3dcoordAxis(xp, xe1[0], xe1[1], xe1[2], true,2);
    G4.draw3dcoordAxis(xp, xe2[0], xe2[1], xe2[2], true,3);

    GP_handle G1("/usr/bin/", "Time (seconds)", "Error");
    G1.gnuplot_cmd("set terminal wxt");
    G1.gnuplot_cmd("plot 'error.txt' u 1 w l t 'Position Error', '' u 2 w l t 'Orientation Error'");

    cout << "Press Enter to exit .." << endl;
    getchar();

#endif


#ifdef JG_TEST




    double xdot[NW];

    double theta[NL] = {0, 0, 0, 0, 0, 0};
    double pose[NW] = {0,0,0,0,0,0};

    robot.set_joint_angles(theta);
    robot.fwd_kin_pose(pose);

    cv::Mat Jg = cv::Mat(NW, NL, CV_64F, 0.0);
    cv::Mat ThDot = cv::Mat(NL,1, CV_64F, 0.0);
    cv::Mat XDot = cv::Mat(NW,1,CV_64F,0.0);

    double theta_dot[NL] = {0, DEG2RAD(10), 0, 0, 0, 0};

    for(int i = 0; i < NL; i++)
         ThDot.at<double>(i) = theta_dot[i];

    double Tmax = 2.0;


    ofstream f1("pose.txt");
    ofstream f2("pose2.txt");
    ofstream f3("angles.txt");

    double x[NW];

    for(int i = 0; i < NW; i++)
        x[i] = pose[i];

    double t = 0.0;
    do
    {

        robot.compute_geometric_jacobian(Jg);

        XDot = Jg * ThDot;



        for(int i = 0; i < NL; i++)
        {
            theta[i] = theta[i] + ThDot.at<double>(i)*dt;

            f3 << theta[i] << "\t";
        }
        f3 << endl;

        robot.set_joint_angles(theta);
        robot.fwd_kin_pose(pose);

        for(int i = 0; i < NW; i++)
        {
            xdot[i] = XDot.at<double>(i);
            x[i] = x[i] + xdot[i]*dt;


            f1 << x[i] << "\t";
            f2 << pose[i] << "\t";
        }
        f1 << endl;
        f2 << endl;

        t = t + dt;

    }while(t <Tmax);
    f1.close();
    f2.close();
    f3.close();

    //------------------------

    GP_handle G1("/usr/bin/", "X (m)", "Y (m)", "Z (m)");
    G1.gnuplot_cmd("set terminal wxt");
    G1.gnuplot_cmd("splot 'pose.txt'  u 1:2:3 w p ps 3 t 'with Jg' , 'pose2.txt' u 1:2:3 w p t 'with theta'");

    getchar();


#endif






    return 0;

}
