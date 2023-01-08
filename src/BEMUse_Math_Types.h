#ifndef VPML_TYPES_H
#define VPML_TYPES_H

#include <Eigen/Eigen>          // Eigen data types
#include <tr1/cmath>            // Special functions
#include <memory>               // Shared ptr.

#define OpenMPfor _Pragma("omp parallel for")

#ifdef SinglePrec

//--------Single Precision----------

typedef float                       Real;
typedef std::complex<float>         CReal;
typedef Eigen::MatrixXf             Matrix;
typedef Eigen::Matrix<float,\
                Eigen::Dynamic,\
                Eigen::Dynamic,\
                Eigen::RowMajor>    MatrixRM;
typedef Eigen::SparseMatrix<float>  SparseMatrix;
typedef Eigen::MatrixXcf            CMatrix;
typedef Eigen::Matrix3f             Matrix3;
typedef Eigen::VectorXf             Vector;
typedef Eigen::RowVectorXf          RowVector;
typedef Eigen::Vector3f             Vector3;
typedef Eigen::Quaternionf          Quat;
typedef Eigen::MatrixXf::Index      Matrix_ID;
typedef std::vector<float>          StdVector;

#endif

#ifdef DoublePrec

//--------Double Precision----------

typedef double                      Real;
typedef std::complex<double>        CReal;
typedef Eigen::MatrixXd             Matrix;
typedef Eigen::Matrix<double,\
                Eigen::Dynamic,\
                Eigen::Dynamic,\
                Eigen::RowMajor>    MatrixRM;
typedef Eigen::SparseMatrix<double> SparseMatrix;
typedef Eigen::MatrixXcd            CMatrix;
typedef Eigen::Matrix3d             Matrix3;
typedef Eigen::VectorXd             Vector;
typedef Eigen::RowVectorXd          RowVector;
typedef Eigen::Vector3d             Vector3;
typedef Eigen::Quaterniond          Quat;
typedef Eigen::MatrixXd::Index      Matrix_ID;
typedef std::vector<double>         StdVector;

#endif

//-------- Cartesian ID types

typedef Eigen::Vector2i             Cart_ID2;

//-------- Additional std::vector variabels and capabilities

typedef std::vector<int>        IntVector;
typedef std::vector<Vector>     StateVector;
static  StateVector             Empty_SV;
static  Vector                  Empty_V;

//-------- Additional deque variables and capabilities

template <class T>          // Hack to avoid using insert everytime
static void StdAppend(std::vector<T>& lhs, const std::vector<T>& rhs)       {lhs.insert(lhs.end(),rhs.begin(), rhs.end());}

template <class T>          // Hack to avoid using insert everytime
static void StdAppend(std::vector<T>& lhs, const std::vector<T>& rhs, const int &P1, const int &P2)     {lhs.insert(lhs.end(),rhs.begin()+P1, rhs.begin()+P2);}

//-------- Time stamp

#include <sys/time.h>
typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

//---- Mathematical Constants

static Real const  EUL          =  0.5772156649;
static Real const  PI           =  M_PI;
static Real const  TwoPI        =  2.0*M_PI;
static Real const  FourPI       =  4.0*M_PI;
static Real const  TwoPIinv     =  0.5/M_PI;
static Real const  FourPIinv    =  0.25/M_PI;

static Real const  D2R          =  M_PI/180;
static Real const  R2D          =  180/M_PI;
static CReal const Im(0.0,1.0);

//---- Kinematic parameters

static Real const  P_atm        = 101325.0;
static Real const  Rho_wat      =  1025;
static Real const  Rho_air      =  1.225;
static Real const  Gravity      =  9.81;
static Real const  Kin_Visc_wat =  1.004e-6;
static Real const  Kin_Visc_air =  1.5571e-5;

#endif // VPML_TYPES_H
