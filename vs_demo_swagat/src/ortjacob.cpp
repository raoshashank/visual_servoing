
#include<ur5model.h>

using namespace std;
using namespace MODUR5;


void UR5::orientation_jacobian(cv::Mat &Jort)
{

    // Rotation Matrix
    rotation_matrix();

    double r21 = R.at<double>(1,0);
    double r11 = R.at<double>(0,0);
    double r31 = R.at<double>(2,0);
    double r32 = R.at<double>(2,1);
    double r33 = R.at<double>(2,2);
    double r22 = R.at<double>(1,1);
    double r12 = R.at<double>(0,1);
    double r13 = R.at<double>(0,2);
    double r23 = R.at<double>(1,2);

    double beta = atan2(-r31, abs(sqrt(r11*r11 + r21*r21)));



    if(beta == M_PI/2.0)
    {
        Jo.at<double>(0,0) = -(((cos(theta[1])*sin(theta[5])+((sin(theta[1])*cos(theta[2])*sin(theta[3])+sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (sin(theta[1])*sin(theta[2])*sin(theta[3])-sin(theta[1])*cos(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[5]))*sin(theta[6])+(
                (sin(theta[1])*sin(theta[2])*sin(theta[3])-sin(theta[1])*cos(theta[2])*cos(theta[3]))*sin(theta[4])+
                (-sin(theta[1])*cos(theta[2])*sin(theta[3])-sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[6]))*r22+(((
                (cos(theta[1])*cos(theta[2])*sin(theta[3])+cos(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (cos(theta[1])*sin(theta[2])*sin(theta[3])-cos(theta[1])*cos(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[5])-sin(theta[1])*sin(theta[5]))*
                sin(theta[6])+((cos(theta[1])*sin(theta[2])*sin(theta[3])-cos(theta[1])*cos(theta[2])*cos(theta[3]))*sin(theta[4])+
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[6]))*r12)*pow(pow(r22,2)+(
                pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+((pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(
                (-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5])*
                sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                );





        Jo.at<double>(0,1) = -((((cos(theta[1])*sin(theta[2])*sin(theta[3])-cos(theta[1])*cos(theta[2])*cos(theta[3]))*sin(theta[4])+
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[5])*sin(theta[6])+(
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (cos(theta[1])*cos(theta[2])*cos(theta[3])-cos(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[6]))*r22+(((-cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*
                pow(sin(theta[4]),2)+(cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[5]),2)+cos(theta[5])*sin(theta[5])*(sin(theta[4])*
                (sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))+cos(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))))*pow(sin(theta[6]),2)+((cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*
                pow(sin(theta[4]),2)+(-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[6]),2)+cos(theta[6])*sin(theta[6])*(cos(theta[5])*((
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(sin(theta[4]),2)+(
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-4*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(4*cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-4*cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))))+
                sin(theta[5])*(cos(theta[4])*(cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2))+sin(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)))))*pow(pow(r22,2)+(
                pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+((pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(
                (-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5])*
                sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                );


        Jo.at<double>(0,2) = -((((cos(theta[1])*sin(theta[2])*sin(theta[3])-cos(theta[1])*cos(theta[2])*cos(theta[3]))*sin(theta[4])+
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[5])*sin(theta[6])+(
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (cos(theta[1])*cos(theta[2])*cos(theta[3])-cos(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[6]))*r22+(((-cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*
                pow(sin(theta[4]),2)+(cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[5]),2)+cos(theta[5])*sin(theta[5])*(sin(theta[4])*
                (sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))+cos(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))))*pow(sin(theta[6]),2)+((cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*
                pow(sin(theta[4]),2)+(-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[6]),2)+cos(theta[6])*sin(theta[6])*(cos(theta[5])*((
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(sin(theta[4]),2)+(
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-4*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(4*cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-4*cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))))+
                sin(theta[5])*(cos(theta[4])*(cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2))+sin(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)))))*pow(pow(r22,2)+(
                pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+((pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(
                (-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5])*
                sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                );


        Jo.at<double>(0,3) = -((((cos(theta[1])*sin(theta[2])*sin(theta[3])-cos(theta[1])*cos(theta[2])*cos(theta[3]))*sin(theta[4])+
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[5])*sin(theta[6])+(
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (cos(theta[1])*cos(theta[2])*cos(theta[3])-cos(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[6]))*r22+(((-cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*
                pow(sin(theta[4]),2)+(cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[5]),2)+cos(theta[5])*sin(theta[5])*(sin(theta[4])*
                (sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))+cos(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))))*pow(sin(theta[6]),2)+((cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*
                pow(sin(theta[4]),2)+(-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[6]),2)+cos(theta[6])*sin(theta[6])*(cos(theta[5])*((
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(sin(theta[4]),2)+(
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-4*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(4*cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-4*cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))))+
                sin(theta[5])*(cos(theta[4])*(cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2))+sin(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)))))*pow(pow(r22,2)+(
                pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+((pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(
                (-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5])*
                sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                );


        Jo.at<double>(0,4) = -pow(pow(r22,2)+(pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+((pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+
                pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*
                pow(sin(theta[4]),2)+(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*
                pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*(2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*
                pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))
                )*pow(cos(theta[5]),2)+((-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*
                sin(theta[4])+(2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*
                cos(theta[5])*sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                )*(((sin(theta[1])*sin(theta[5])+((-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (cos(theta[1])*cos(theta[2])*cos(theta[3])-cos(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5]))*pow(sin(theta[6]),2)+(
                (cos(theta[1])*cos(theta[2])*cos(theta[3])-cos(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (cos(theta[1])*cos(theta[2])*sin(theta[3])+cos(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[6])*sin(theta[6]))*r23+sin(theta[6])*
                r13*r22);



        Jo.at<double>(0,5) = -(r11*r22+(cos(theta[5])*((-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*
                cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*pow(sin(theta[4]),2)+(cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*
                pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))+sin(theta[5])*(sin(theta[4])*
                (sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))+cos(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))))*pow(sin(theta[6]),2)+(cos(theta[5])*((
                cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*
                pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))+sin(theta[5])*(sin(theta[4])*
                (sin(theta[2])*sin(theta[3])*pow(cos(theta[1]),2)-cos(theta[2])*cos(theta[3])*pow(cos(theta[1]),2))+cos(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(cos(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(cos(theta[1]),2))))*pow(cos(theta[6]),2)+cos(theta[6])*
                sin(theta[6])*(-cos(theta[1])*sin(theta[1])*pow(sin(theta[5]),2)+((cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+
                cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*
                sin(theta[3]))*pow(sin(theta[4]),2)+(cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*
                pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-2*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-2*cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(-cos(theta[1])*
                sin(theta[1])*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(sin(theta[4]),2)+(-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)*
                pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(sin(theta[3]),2)-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-2*cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))+cos(theta[5])*sin(theta[5])*(
                cos(theta[4])*(cos(theta[2])*cos(theta[3])*(pow(sin(theta[1]),2)-pow(cos(theta[1]),2))+sin(theta[2])*sin(theta[3])*
                (pow(cos(theta[1]),2)-pow(sin(theta[1]),2)))+sin(theta[4])*(cos(theta[2])*sin(theta[3])*(pow(cos(theta[1]),2)-pow(sin(theta[1]),2))+
                sin(theta[2])*cos(theta[3])*(pow(cos(theta[1]),2)-pow(sin(theta[1]),2))))))*pow(pow(r22,2)+(pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+(
                (pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)+(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*
                pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*
                pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*
                pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(
                (-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5])*
                sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                );


        //----------------------------
        Jo.at<double>(1,0) = 0.0;
        Jo.at<double>(1,1) = 0.0;
        Jo.at<double>(1,2) = 0.0;
        Jo.at<double>(1,3) = 0.0;
        Jo.at<double>(1,4) = 0.0;
        Jo.at<double>(1,5) = 0.0;


        Jo.at<double>(2,0) = 0.0;
        Jo.at<double>(2,1) = 0.0;
        Jo.at<double>(2,2) = 0.0;
        Jo.at<double>(2,3) = 0.0;
        Jo.at<double>(2,4) = 0.0;
        Jo.at<double>(2,5) = 0.0;

    }
    else if (Beta == -M_PI/2.0)
    {

        //--------------------------------

        Jo.at<double>(0,0) = (((cos(theta[1])*sin(theta[5])+((sin(theta[1])*cos(theta[2])*sin(theta[3])+sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (sin(theta[1])*sin(theta[2])*sin(theta[3])-sin(theta[1])*cos(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[5]))*sin(theta[6])+(
                (sin(theta[1])*sin(theta[2])*sin(theta[3])-sin(theta[1])*cos(theta[2])*cos(theta[3]))*sin(theta[4])+
                (-sin(theta[1])*cos(theta[2])*sin(theta[3])-sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[6]))*r22+(((
                (cos(theta[1])*cos(theta[2])*sin(theta[3])+cos(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (cos(theta[1])*sin(theta[2])*sin(theta[3])-cos(theta[1])*cos(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[5])-sin(theta[1])*sin(theta[5]))*
                sin(theta[6])+((cos(theta[1])*sin(theta[2])*sin(theta[3])-cos(theta[1])*cos(theta[2])*cos(theta[3]))*sin(theta[4])+
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[6]))*r12)*pow(pow(r22,2)+(
                pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+((pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(
                (-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5])*
                sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                );


        Jo.at<double>(0,1) = ((((cos(theta[1])*sin(theta[2])*sin(theta[3])-cos(theta[1])*cos(theta[2])*cos(theta[3]))*sin(theta[4])+
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[5])*sin(theta[6])+(
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (cos(theta[1])*cos(theta[2])*cos(theta[3])-cos(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[6]))*r22+(((-cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*
                pow(sin(theta[4]),2)+(cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[5]),2)+cos(theta[5])*sin(theta[5])*(sin(theta[4])*
                (sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))+cos(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))))*pow(sin(theta[6]),2)+((cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*
                pow(sin(theta[4]),2)+(-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[6]),2)+cos(theta[6])*sin(theta[6])*(cos(theta[5])*((
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(sin(theta[4]),2)+(
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-4*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(4*cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-4*cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))))+
                sin(theta[5])*(cos(theta[4])*(cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2))+sin(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)))))*pow(pow(r22,2)+(
                pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+((pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(
                (-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5])*
                sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                );






        Jo.at<double>(0,2) = ((((cos(theta[1])*sin(theta[2])*sin(theta[3])-cos(theta[1])*cos(theta[2])*cos(theta[3]))*sin(theta[4])+
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[5])*sin(theta[6])+(
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (cos(theta[1])*cos(theta[2])*cos(theta[3])-cos(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[6]))*r22+(((-cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*
                pow(sin(theta[4]),2)+(cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[5]),2)+cos(theta[5])*sin(theta[5])*(sin(theta[4])*
                (sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))+cos(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))))*pow(sin(theta[6]),2)+((cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*
                pow(sin(theta[4]),2)+(-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[6]),2)+cos(theta[6])*sin(theta[6])*(cos(theta[5])*((
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(sin(theta[4]),2)+(
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-4*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(4*cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-4*cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))))+
                sin(theta[5])*(cos(theta[4])*(cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2))+sin(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)))))*pow(pow(r22,2)+(
                pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+((pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(
                (-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5])*
                sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                );





        Jo.at<double>(0,3) = ((((cos(theta[1])*sin(theta[2])*sin(theta[3])-cos(theta[1])*cos(theta[2])*cos(theta[3]))*sin(theta[4])+
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[5])*sin(theta[6])+(
                (-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (cos(theta[1])*cos(theta[2])*cos(theta[3])-cos(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[6]))*r22+(((-cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*
                pow(sin(theta[4]),2)+(cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[5]),2)+cos(theta[5])*sin(theta[5])*(sin(theta[4])*
                (sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))+cos(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))))*pow(sin(theta[6]),2)+((cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*
                pow(sin(theta[4]),2)+(-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))*pow(cos(theta[6]),2)+cos(theta[6])*sin(theta[6])*(cos(theta[5])*((
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(sin(theta[4]),2)+(
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-4*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(4*cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-4*cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))))+
                sin(theta[5])*(cos(theta[4])*(cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2))+sin(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2)))))*pow(pow(r22,2)+(
                pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+((pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(
                (-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5])*
                sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                );


        Jo.at<double>(0,4) = pow(pow(r22,2)+(pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+((pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+
                pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*
                pow(sin(theta[4]),2)+(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*
                pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*(2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*
                pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))
                )*pow(cos(theta[5]),2)+((-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*
                sin(theta[4])+(2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*
                cos(theta[5])*sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                )*(((sin(theta[1])*sin(theta[5])+((-cos(theta[1])*cos(theta[2])*sin(theta[3])-cos(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (cos(theta[1])*cos(theta[2])*cos(theta[3])-cos(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5]))*pow(sin(theta[6]),2)+(
                (cos(theta[1])*cos(theta[2])*cos(theta[3])-cos(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (cos(theta[1])*cos(theta[2])*sin(theta[3])+cos(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*cos(theta[6])*sin(theta[6]))*r23+sin(theta[6])*
                r13*r22);

        Jo.at<double>(0,5) = (r11*r22+(cos(theta[5])*((-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*
                cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*pow(sin(theta[4]),2)+(cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*
                pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(cos(theta[3]),2)+4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))+sin(theta[5])*(sin(theta[4])*
                (sin(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-cos(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))+cos(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(sin(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(sin(theta[1]),2))))*pow(sin(theta[6]),2)+(cos(theta[5])*((
                cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+
                cos(theta[3])*sin(theta[3])*(cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)))*
                pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(
                (cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-4*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])))+sin(theta[5])*(sin(theta[4])*
                (sin(theta[2])*sin(theta[3])*pow(cos(theta[1]),2)-cos(theta[2])*cos(theta[3])*pow(cos(theta[1]),2))+cos(theta[4])*
                (-cos(theta[2])*sin(theta[3])*pow(cos(theta[1]),2)-sin(theta[2])*cos(theta[3])*pow(cos(theta[1]),2))))*pow(cos(theta[6]),2)+cos(theta[6])*
                sin(theta[6])*(-cos(theta[1])*sin(theta[1])*pow(sin(theta[5]),2)+((cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+
                cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*
                sin(theta[3]))*pow(sin(theta[4]),2)+(cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+cos(theta[1])*sin(theta[1])*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*
                pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(sin(theta[3]),2)-2*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-2*cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(-cos(theta[1])*
                sin(theta[1])*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)+2*cos(theta[1])*
                sin(theta[1])*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(sin(theta[4]),2)+(-cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)*
                pow(sin(theta[3]),2)-cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*cos(theta[3])*sin(theta[3]))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[1])*sin(theta[1])*cos(theta[2])*
                sin(theta[2])*pow(sin(theta[3]),2)-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[2])*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*cos(theta[1])*sin(theta[1])*pow(sin(theta[2]),2)-2*cos(theta[1])*sin(theta[1])*pow(cos(theta[2]),2)))+cos(theta[5])*sin(theta[5])*(
                cos(theta[4])*(cos(theta[2])*cos(theta[3])*(pow(sin(theta[1]),2)-pow(cos(theta[1]),2))+sin(theta[2])*sin(theta[3])*
                (pow(cos(theta[1]),2)-pow(sin(theta[1]),2)))+sin(theta[4])*(cos(theta[2])*sin(theta[3])*(pow(cos(theta[1]),2)-pow(sin(theta[1]),2))+
                sin(theta[2])*cos(theta[3])*(pow(cos(theta[1]),2)-pow(sin(theta[1]),2))))))*pow(pow(r22,2)+(pow(sin(theta[1]),2)*pow(sin(theta[5]),2)+(
                (pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)+(pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*
                pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*
                pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*
                pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))))*pow(cos(theta[5]),2)+(
                (-2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*cos(theta[4]))*cos(theta[5])*
                sin(theta[5]))*pow(sin(theta[6]),2)+((pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*
                pow(cos(theta[2]),2)*pow(cos(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(sin(theta[4]),2)
                +(pow(cos(theta[1]),2)*pow(cos(theta[2]),2)*pow(sin(theta[3]),2)+pow(cos(theta[1]),2)*pow(sin(theta[2]),2)*pow(cos(theta[3]),2)+2*
                cos(theta[2])*sin(theta[2])*cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2))*pow(cos(theta[4]),2)+cos(theta[4])*sin(theta[4])*(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))))*pow(cos(theta[6]),2)+
                cos(theta[6])*sin(theta[6])*(cos(theta[5])*((2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)-2*cos(theta[2])*sin(theta[2])*
                pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*sin(theta[3])*
                (2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)))*pow(sin(theta[4]),2)+(-2*cos(theta[2])*
                sin(theta[2])*pow(cos(theta[1]),2)*pow(sin(theta[3]),2)+2*cos(theta[2])*sin(theta[2])*pow(cos(theta[1]),2)*pow(cos(theta[3]),2)+cos(theta[3])*
                sin(theta[3])*(2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)))*pow(cos(theta[4]),2)+cos(theta[4])*
                sin(theta[4])*((2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2))*pow(sin(theta[3]),2)+
                (2*pow(cos(theta[1]),2)*pow(cos(theta[2]),2)-2*pow(cos(theta[1]),2)*pow(sin(theta[2]),2))*pow(cos(theta[3]),2)-8*cos(theta[2])*sin(theta[2])*
                cos(theta[3])*sin(theta[3])*pow(cos(theta[1]),2)))+(
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*cos(theta[3])-2*cos(theta[1])*sin(theta[1])*sin(theta[2])*sin(theta[3]))*sin(theta[4])+
                (2*cos(theta[1])*sin(theta[1])*cos(theta[2])*sin(theta[3])+2*cos(theta[1])*sin(theta[1])*sin(theta[2])*cos(theta[3]))*cos(theta[4]))*sin(theta[5])),-1
                );




        //------------------------------------

        Jo.at<double>(1,0) = 0.0;
        Jo.at<double>(1,1) = 0.0;
        Jo.at<double>(1,2) = 0.0;
        Jo.at<double>(1,3) = 0.0;
        Jo.at<double>(1,4) = 0.0;
        Jo.at<double>(1,5) = 0.0;


        Jo.at<double>(2,0) = 0.0;
        Jo.at<double>(2,1) = 0.0;
        Jo.at<double>(2,2) = 0.0;
        Jo.at<double>(2,3) = 0.0;
        Jo.at<double>(2,4) = 0.0;
        Jo.at<double>(2,5) = 0.0;


    }
    else // if beta > -M_PI/2.0 && beta < M_PI/2.0
    {

        Jo.at<double>(0,0) = 0.0;

        Jo.at<double>(0,1) = ((8*sin(theta[5]+theta[4]+theta[3]+theta[2])+8*sin(theta[5]-theta[4]-theta[3]-theta[2]))*r32-
                    sin(theta[6]+2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+sin(theta[6]+2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-2*
                    sin(theta[6]+theta[5]+2*theta[4]+2*theta[3]+2*theta[2])-2*sin(theta[6]+theta[5]-2*theta[4]-2*theta[3]-2*theta[2])+4*sin(theta[6]+theta[5])+2*
                    sin(theta[6]-theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+2*sin(theta[6]-theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-4*sin(theta[6]-theta[5])+
                    sin(theta[6]-2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])-sin(theta[6]-2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2]))/(16*r32^2+2*
                    cos(2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+2*cos(2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-4*cos(2*theta[5])-4*
                    cos(2*theta[4]+2*theta[3]+2*theta[2])+4);

        Jo.at<double>(0,2) = ((8*sin(theta[5]+theta[4]+theta[3]+theta[2])+8*sin(theta[5]-theta[4]-theta[3]-theta[2]))*r32-
                    sin(theta[6]+2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+sin(theta[6]+2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-2*
                    sin(theta[6]+theta[5]+2*theta[4]+2*theta[3]+2*theta[2])-2*sin(theta[6]+theta[5]-2*theta[4]-2*theta[3]-2*theta[2])+4*sin(theta[6]+theta[5])+2*
                    sin(theta[6]-theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+2*sin(theta[6]-theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-4*sin(theta[6]-theta[5])+
                    sin(theta[6]-2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])-sin(theta[6]-2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2]))/(16*r32^2+2*
                    cos(2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+2*cos(2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-4*cos(2*theta[5])-4*
                    cos(2*theta[4]+2*theta[3]+2*theta[2])+4);


        Jo.at<double>(0,3) = ((8*sin(theta[5]+theta[4]+theta[3]+theta[2])+8*sin(theta[5]-theta[4]-theta[3]-theta[2]))*r32-
                    sin(theta[6]+2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+sin(theta[6]+2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-2*
                    sin(theta[6]+theta[5]+2*theta[4]+2*theta[3]+2*theta[2])-2*sin(theta[6]+theta[5]-2*theta[4]-2*theta[3]-2*theta[2])+4*sin(theta[6]+theta[5])+2*
                    sin(theta[6]-theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+2*sin(theta[6]-theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-4*sin(theta[6]-theta[5])+
                    sin(theta[6]-2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])-sin(theta[6]-2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2]))/(16*r32^2+2*
                    cos(2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+2*cos(2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-4*cos(2*theta[5])-4*
                    cos(2*theta[4]+2*theta[3]+2*theta[2])+4);


        Jo.at<double>(0,4) = ((8*sin(theta[5]+theta[4]+theta[3]+theta[2])-8*sin(theta[5]-theta[4]-theta[3]-theta[2]))*r32-
                    sin(theta[6]+2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])-sin(theta[6]+2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])+2*sin(theta[6]+2*theta[5])-
                    sin(theta[6]-2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])-sin(theta[6]-2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])+2*sin(theta[6]-2*theta[5])+2*
                    sin(theta[6]+2*theta[4]+2*theta[3]+2*theta[2])+2*sin(theta[6]-2*theta[4]-2*theta[3]-2*theta[2])-4*sin(theta[6]))/(16*r32^2+2*
                    cos(2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+2*cos(2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-4*cos(2*theta[5])-4*
                    cos(2*theta[4]+2*theta[3]+2*theta[2])+4);


        Jo.at<double>(0,5) = -(sin(theta[6]+2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+sin(theta[6]+2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-2*
                                                                                                                                                   sin(theta[6]+2*theta[5])+2*sin(theta[6]+theta[5]+2*theta[4]+2*theta[3]+2*theta[2])-2*sin(theta[6]+theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-2*
                                                                                                                                                   sin(theta[6]-theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+2*sin(theta[6]-theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-
                                                                                                                                                   sin(theta[6]-2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])-sin(theta[6]-2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])+2*sin(theta[6]-2*theta[5]))/(16*
                                                                                                                                                   r32^2+2*cos(2*theta[5]+2*theta[4]+2*theta[3]+2*theta[2])+2*cos(2*theta[5]-2*theta[4]-2*theta[3]-2*theta[2])-4*cos(2*theta[5])-4*
                                                                                                                                                   cos(2*theta[4]+2*theta[3]+2*theta[2])+4);

        //----------------

        Jo.at<double>(1,0) = (8*r11*r21+(-sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])+sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+2*sin(theta[6]+theta[5]+theta[1])
                    +2*sin(theta[6]+theta[5]-theta[1])-sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])+
                    sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])-2*sin(theta[6]-theta[5]+theta[1])-2*sin(theta[6]-theta[5]-theta[1])-2*
                    sin(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])+2*sin(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])+2*
                    sin(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])-2*sin(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r11)*pow(pow(r21,2)+pow(r11,2),1/2)*r31*
                    pow((8*pow(r21,2)+8*pow(r11,2))*pow(r31,2)+8*pow(r21,4)+16*pow(r11,2)*pow(r21,2)+8*pow(r11,4),-1);

        Jo.at<double>(1,1) = pow(pow(r21,2)+pow(r11,2),1/2)*(((cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+2*cos(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])-2*
                    cos(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])+2*cos(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])-2*
                    cos(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r21+(-sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])+sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])-sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])+sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])-2*sin(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])-2*
                    sin(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])-2*sin(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])-2*
                    sin(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r11)*r31+(-2*cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2])-2*
                    cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2])-2*cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2])-2*
                    cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2])-4*cos(theta[6]+theta[4]+theta[3]+theta[2])+4*cos(theta[6]-theta[4]-theta[3]-theta[2]))*pow(r21,2)+
                    (-2*cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2])-2*cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2])-2*
                    cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2])-2*cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2])-4*cos(theta[6]+theta[4]+theta[3]+theta[2])+4*
                    cos(theta[6]-theta[4]-theta[3]-theta[2]))*pow(r11,2))*
                    pow((8*pow(r21,2)+8*pow(r11,2))*pow(r31,2)+8*pow(r21,4)+16*pow(r11,2)*pow(r21,2)+8*pow(r11,4),-1);


        Jo.at<double>(1,2) =pow(pow(r21,2)+pow(r11,2),1/2)*(((cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+2*cos(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])-2*
                    cos(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])+2*cos(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])-2*
                    cos(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r21+(-sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])+sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])-sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])+sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])-2*sin(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])-2*
                    sin(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])-2*sin(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])-2*
                    sin(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r11)*r31+(-2*cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2])-2*
                    cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2])-2*cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2])-2*
                    cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2])-4*cos(theta[6]+theta[4]+theta[3]+theta[2])+4*cos(theta[6]-theta[4]-theta[3]-theta[2]))*pow(r21,2)+
                    (-2*cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2])-2*cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2])-2*
                    cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2])-2*cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2])-4*cos(theta[6]+theta[4]+theta[3]+theta[2])+4*
                    cos(theta[6]-theta[4]-theta[3]-theta[2]))*pow(r11,2))*
                    pow((8*pow(r21,2)+8*pow(r11,2))*pow(r31,2)+8*pow(r21,4)+16*pow(r11,2)*pow(r21,2)+8*pow(r11,4),-1);

        Jo.at<double>(1,3) = pow(pow(r21,2)+pow(r11,2),1/2)*(((cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+2*cos(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])-2*
                    cos(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])+2*cos(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])-2*
                    cos(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r21+(-sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])+sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])-sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-
                    sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])+sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+
                    sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])-2*sin(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])-2*
                    sin(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])-2*sin(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])-2*
                    sin(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r11)*r31+(-2*cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2])-2*
                    cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2])-2*cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2])-2*
                    cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2])-4*cos(theta[6]+theta[4]+theta[3]+theta[2])+4*cos(theta[6]-theta[4]-theta[3]-theta[2]))*pow(r21,2)+
                    (-2*cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2])-2*cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2])-2*
                    cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2])-2*cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2])-4*cos(theta[6]+theta[4]+theta[3]+theta[2])+4*
                    cos(theta[6]-theta[4]-theta[3]-theta[2]))*pow(r11,2))*
                    pow((8*pow(r21,2)+8*pow(r11,2))*pow(r31,2)+8*pow(r21,4)+16*pow(r11,2)*pow(r21,2)+8*pow(r11,4),-1);


        Jo.at<double>(1,4) = pow(pow(r21,2)+pow(r11,2),1/2)*((4*cos(theta[6])*r21*r23+4*cos(theta[6])*r11*r13)*r31+(-
                                                                                                                        cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2])+cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2])+cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2])-
                                                                                                                        cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]))*pow(r21,2)+(-cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2])+
                                                                                                                        cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2])+cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2])-cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]))*
                                                                                                                        pow(r11,2))*pow((4*pow(r21,2)+4*pow(r11,2))*pow(r31,2)+4*pow(r21,4)+8*pow(r11,2)*pow(r21,2)+4*pow(r11,4),-1);



        Jo.at<double>(1,5) = -pow(pow(r21,2)+pow(r11,2),1/2)*
                              pow((pow(r21,2)+pow(r11,2))*pow(r31,2)+pow(r21,4)+2*pow(r11,2)*pow(r21,2)+pow(r11,4),-1)*
                              ((pow(r21,2)+pow(r11,2))*r32+(-r21*r22-r11*r12)*r31);

        //----------------------

        Jo.at<double>(2,0) = ((sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])+
                sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])-sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])-2*sin(theta[6]+theta[5]+theta[1])-
                2*sin(theta[6]+theta[5]-theta[1])+sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])
                +sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])-sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+2*sin(theta[6]-theta[5]+theta[1])
                +2*sin(theta[6]-theta[5]-theta[1])+2*sin(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])-2*sin(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])-2*
                sin(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])+2*sin(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r21+8*pow(r11,2))*
                pow(8*pow(r21,2)+8*pow(r11,2),-1);


        Jo.at<double>(2,1) = ((sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])+sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])-sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+
                    sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])+sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])-sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+2*
                    sin(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])+2*sin(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])+2*
                    sin(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])+2*sin(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r21+(
                    cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+
                    cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+2*
                    cos(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])-2*cos(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])+2*
                    cos(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])-2*cos(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r11)*pow(8*pow(r21,2)+8*pow(r11,2),-1);

        Jo.at<double>(2,2) = ((sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])+sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])-sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+
                    sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])+sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])-sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+2*
                    sin(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])+2*sin(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])+2*
                    sin(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])+2*sin(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r21+(
                    cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+
                    cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+2*
                    cos(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])-2*cos(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])+2*
                    cos(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])-2*cos(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r11)*pow(8*pow(r21,2)+8*pow(r11,2),-1);


        Jo.at<double>(2,3) = ((sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])+sin(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])-sin(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+
                    sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])+sin(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])-sin(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+2*
                    sin(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])+2*sin(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])+2*
                    sin(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])+2*sin(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r21+(
                    cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-cos(theta[6]+theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+cos(theta[6]+theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+
                    cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]+theta[1])-cos(theta[6]-theta[5]+theta[4]+theta[3]+theta[2]-theta[1])-
                    cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]+theta[1])+cos(theta[6]-theta[5]-theta[4]-theta[3]-theta[2]-theta[1])+2*
                    cos(theta[6]+theta[4]+theta[3]+theta[2]+theta[1])-2*cos(theta[6]+theta[4]+theta[3]+theta[2]-theta[1])+2*
                    cos(theta[6]-theta[4]-theta[3]-theta[2]+theta[1])-2*cos(theta[6]-theta[4]-theta[3]-theta[2]-theta[1]))*r11)*pow(8*pow(r21,2)+8*pow(r11,2),-1);

        Jo.at<double>(2,4) = pow(pow(r21,2)+pow(r11,2),-1)*(cos(theta[6])*r11*r23-cos(theta[6])*r13*r21);

        Jo.at<double>(2,5) =pow(pow(r21,2)+pow(r11,2),-1)*(r11*r22-r12*r21);

    }


    //-------------------------------

    if(&Jort != NULL)
    {
        for(int i = 0; i < NC; i++)
            for(int j = 0; j < NL; j++)
                Jort.at<double>(i,j) = Jo.at<double>(i,j);
    }

}
