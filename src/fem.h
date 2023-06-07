//
// Created by Roman Wallner- Silberhuber on 22.04.23.
//



#ifndef SIMPLEFEMCPP_FEM_H
#define SIMPLEFEMCPP_FEM_H

// Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>

// C++
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <tuple>



// Fast Matrix Market
#include "fast_matrix_market/fast_matrix_market.hpp"
#include "fast_matrix_market/app/Eigen.hpp"
#include "fem.h"
//!  A test class.
/*!
  A more elaborate class description.
*/

class FemParams
{
public:
    /** An enum type.
     *  The documentation block cannot be put after the enum!
     */
    int mesh_ex;  /*!< Detailed description after the member */
    int mesh_ey;  /*!< Detailed description after the member */
    double mesh_lx;
    double mesh_ly;
    int mesh_nx;
    int mesh_ny;
    int num_nodes;
    int num_elements;
    double mesh_hx;
    double mesh_hy;

    //! A constructor.
    /*!
      A more elaborate description of the constructor.
    */
    FemParams(int ex = 9, int ey = 49, double lx = 10.0, double ly = 50.0);
};

class MaterialParams
{
public:
    double E;
    double v;

    MaterialParams(double c_E = 100.0, double c_v = 0.48);
};


template<typename T, int Rows, int Cols>
std::array<std::array<T, Cols>, Rows> eigenArrayToCppArray( const Eigen::Array<T, Rows, Cols>& eigenArray );


Eigen::Vector4d shape(double xi, double eta);


Eigen::Matrix<double, 2, 4> grad_shape(double xi, double eta);


std::tuple <Eigen::MatrixXd, Eigen::MatrixXi> create_mesh(const FemParams &fem_params);


Eigen::Matrix3d material_matrix(const MaterialParams &mat_params);

#endif //SIMPLEFEMCPP_FEM_H
