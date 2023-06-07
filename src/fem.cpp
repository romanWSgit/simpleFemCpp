// header
#include "fem.h"

// Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>

// C++
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <tuple>

// VTK :
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkQuad.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridWriter.h>

// Fast Matrix Market
#include "fast_matrix_market/fast_matrix_market.hpp"
#include "fast_matrix_market/app/Eigen.hpp"
#include "fem.h"
FemParams::FemParams(int ex, int ey, double lx, double ly)
{
    mesh_ex = ex;
    mesh_ey = ey;
    mesh_lx = lx;
    mesh_ly = ly;
    mesh_nx = ex + 1;
    mesh_ny = ey + 1;
    num_nodes = (ex + 1) * (ey + 1);
    num_elements = ex * ey;
    mesh_hx = lx / ex;
    mesh_hy = ly / ey;
}

MaterialParams::MaterialParams(double c_E, double c_v)
{
    E = c_E;
    v = c_v;
}


template<typename T, int Rows, int Cols>
std::array<std::array<T, Cols>, Rows> eigenArrayToCppArray( const Eigen::Array<T, Rows, Cols>& eigenArray )
{
    std::array<std::array<T, Cols>, Rows> cppArray;
    for (int i = 0; i < Rows; i++) {
        for (int j = 0; j < Cols; j++) {
            cppArray[i][j] = eigenArray(i, j);
        }
    }
    return cppArray;
}


Eigen::Vector4d shape(double xi, double eta) {
    Eigen::Vector4d N;
    N <<    (1.0-xi)*(1.0-eta), (1.0+xi)*(1.0-eta),(1.0+xi)*(1.0+eta), (1.0-xi)*(1.0+eta);
    N *= 0.25;
    return N;
}

Eigen::Matrix<double, 2, 4> grad_shape(double xi, double eta) {
    Eigen::Matrix<double, 2, 4> dN;
    dN << -(1.0 - eta), (1.0 - eta), (1.0 + eta), -(1.0 + eta),
            -(1.0 - xi), -(1.0 + xi), (1.0 + xi), (1.0 - xi);
    dN *= 0.25;
    return dN;
}

std::tuple <Eigen::MatrixXd, Eigen::MatrixXi> create_mesh(const FemParams &fem_params)
{

    std::cout << "number of nodes: " << fem_params.num_nodes << std::endl;

    // Nodes:

    Eigen::MatrixXd nodes(fem_params.num_nodes , 2);
    int cnt_rows = 0;
    for (auto y : Eigen::ArrayXd::LinSpaced(50, 0.0, 50.0)) {
        for (auto x : Eigen::ArrayXd::LinSpaced(10, 0.0, 10.0)) {
            nodes(cnt_rows, 0) = x;
            nodes(cnt_rows, 1) = y;
            cnt_rows++;
        }
    }


    // Connectivity

    Eigen::MatrixXi conn(fem_params.num_elements , 4);
    int n0;
    cnt_rows = 0;
    for (int j = 0; j < fem_params.mesh_ey; j++) {
        for (int i = 0; i < fem_params.mesh_ex; i++) {
            n0 = i + j * fem_params.mesh_nx;
            conn(cnt_rows, 0) = n0;
            conn(cnt_rows, 1) = n0 + 1;
            conn(cnt_rows, 2) = n0 + 1 + fem_params.mesh_nx;
            conn(cnt_rows, 3) = n0 + fem_params.mesh_nx;
            cnt_rows++;
        }
    }

    return {nodes, conn};
}

Eigen::Matrix3d material_matrix(const MaterialParams &mat_params) {
    Eigen::Matrix3d D;
    double poisson_ratio = mat_params.v;
    double young_modulus = mat_params.E;

    D <<
        1.0 - poisson_ratio,     poisson_ratio,     0.0,
        poisson_ratio, 1.0-poisson_ratio,     0.0,
        0.0,   0.0,   0.5-poisson_ratio;
    D *= young_modulus/((1.0+poisson_ratio)*(1.0-2.0*poisson_ratio)) ;
    return D;
}

int main() {
    std::system("(cd .. && mkdir plot)");
    std::system("(cd .. && mkdir data)");
    std::system("(cd .. && mkdir ParaView)");



    FemParams fem_params;
    MaterialParams mat_params;

    std::cout << fem_params.mesh_nx << std::endl;

    Eigen::Matrix3d D = material_matrix(mat_params);



    Eigen::Array<double, 500 ,2> nodes;
    Eigen::Array<int, 441, 4> connectivity;
    std::tie(nodes, connectivity)  = create_mesh(fem_params);



    std::cout << "create global stiffness matrix... " << std::endl;
    Eigen::MatrixXd K(2 * fem_params.num_nodes,  2 * fem_params.num_nodes);
    K.setZero();

    Eigen::SparseMatrix<double> K_global(2 * fem_params.num_nodes, 2 * fem_params.num_nodes);
    std::vector<Eigen::Triplet<double> > triplets;



    Eigen::Matrix<double, 4 , 2> q4;
    int cnt_rows = 0;
    for (double y : {-1.0, 1.0}) {
        for (double x : {-1.0, 1.0}) {
            Eigen::Vector2d coord_V(x/sqrt(3.0), y/sqrt(3.0));
            q4.row(cnt_rows) = coord_V;
            cnt_rows++;
        }
    }


    Eigen::Matrix<double, 3 , 8> B = Eigen::Matrix<double, 3 , 8>::Zero();

    for (int r = 0; r < connectivity.rows(); r++) {
        Eigen::Matrix<double, 4 , 2> x_Ie = nodes(connectivity.row(r),Eigen::all);
        Eigen::Matrix<double, 8 , 8> K_e = Eigen::Matrix<double, 8 , 8>::Zero();

        for (int q = 0; q < q4.rows(); q++) {

            Eigen::Matrix<double, 2, 4> d_N = grad_shape(q4.row(q)(0), q4.row(q)(1));
            Eigen::Matrix2d J = (d_N*x_Ie).transpose();

            Eigen::Matrix<double, 2, 4> d_N_new = J.inverse() * d_N;

            B(0,Eigen::seq(0,7,2)) = d_N_new.row(0);
            B(1,Eigen::seq(1,7,2))= d_N_new.row(1);
            B(2,Eigen::seq(0,7,2)) = d_N_new.row(1);
            B(2,Eigen::seq(1,7,2)) = d_N_new.row(0);


            K_e +=  ((B.transpose()) * D * B)  * J.determinant();    // np.dot(np.dot(B.T,C),B) * np.linalg.det(J)

        }

        Eigen::Vector4i c_index = connectivity.row(r);
        for (int i = 0; i < c_index.size(); i++) {
            int I = c_index[i];
            for (int j = 0; j < c_index.size(); j++) {
                int J = c_index[j];
                K(2*I,2*J)     += K_e(2*i,2*j);
                K(2*I+1,2*J)   += K_e(2*i+1,2*j);
                K(2*I+1,2*J+1) += K_e(2*i+1,2*j+1);
                K(2*I,2*J+1)   += K_e(2*i,2*j+1);
                Eigen::Triplet<double> trplt_11(2 * I + 0, 2 * J + 0, K_e(2*i,2*j));
                Eigen::Triplet<double> trplt_12(2 * I + 0, 2 * J + 1, K_e(2*i,2*j+1));
                Eigen::Triplet<double> trplt_21(2 * I + 1, 2 * J + 0, K_e(2*i+1,2*j));
                Eigen::Triplet<double> trplt_22(2 * I + 1, 2 * J + 1, K_e(2*i+1,2*j+1));

                triplets.push_back(trplt_11);
                triplets.push_back(trplt_12);
                triplets.push_back(trplt_21);
                triplets.push_back(trplt_22);
            }
        }
    }







    std::cout << "assign nodal forces and boundary conditions... " << std::endl;
    K_global.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::MatrixXd F(2 * fem_params.num_nodes, 1);
    F.setZero();
    std::vector<int> indicesToConstraint;
    for (int i=0; i < fem_params.num_nodes; i++) {
        if ((nodes.row(i)(1) < 0.001) && (nodes.row(i)(1) > -0.001)) {
            K.col(2*i).setZero();
            K.col(2*i + 1).setZero();
            K_global.col(2*i) *= 0.0;
            K_global.col(2*i + 1) *= 0.0;

            K(2*i,2*i)   = 1.0;
            K(2*i+1,2*i+1) = 1.0;

        }
        if (nodes(i,1) == fem_params.mesh_ly) {
            double x =  nodes(i,0);
            F(2*i+1) = 20.0;
            if (x == 0.0 or x == fem_params.mesh_lx) {
                F(2*i+1) *= 0.5;
            }
        }
    }
    // Sparce Matrix Boundary Conditions
    for (int k = 0; k < K_global.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K_global, k); it; ++it) {
            for (std::vector<int>::iterator idit = indicesToConstraint.begin(); idit != indicesToConstraint.end(); ++idit) {
                if (it.row() == *idit || it.col() == *idit) {
                    it.valueRef() = it.row() == it.col() ? 1.0 : 0.0;
                }
            }
        }
    }



    Eigen::MatrixXd dense_from_sparse_KMat;
    dense_from_sparse_KMat = Eigen::MatrixXd(K_global);

    // write stiffness Matrices to mtx
    {
        std::ofstream K_out("K.mtx");
        fast_matrix_market::write_matrix_market_eigen_dense(K_out,K);

        std::ofstream K_out2("dense_from_spare_K.mtx");
        fast_matrix_market::write_matrix_market_eigen_dense(K_out2,dense_from_sparse_KMat);

    }


    std::cout << "solve linear system... " << std::endl;
    Eigen::VectorXd u(2 * fem_params.num_nodes);
    u = K.ldlt().solve(F);

    std::cout << "solve sparce linear system... " << std::endl;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver(K_global);
    Eigen::VectorXd displacements = solver.solve(F);

    double u_max = u.maxCoeff();
    double u_min = u.minCoeff();
    std::cout << "umax: "<< u_max <<std::endl;
    std::cout << "umin: "<< u_min <<std::endl;
    double displ_max = displacements.maxCoeff();
    double displ_min = displacements.minCoeff();
    std::cout << "displmax: "<< displ_max <<std::endl;
    std::cout << "displmin: "<< displ_min <<std::endl;

    Eigen::MatrixXd u_matrix2d(fem_params.num_nodes, 2);
    std::vector<std::vector<double>> u_array2d(fem_params.num_nodes, std::vector<double>(2));
    std::vector<std::vector<double>> nodes_puls_displ_2d(fem_params.num_nodes, std::vector<double>(2));
//    std::array<std::array<double, 2>, fem_params.num_nodes> u_array2d;
//    std::array<std::array<double, 2>, fem_params.num_nodes> nodes_puls_displ_2d;


    for (int i = 0; i < 2*fem_params.num_nodes; i += 2) {
        int row = i / 2;
        u_array2d[row][0] = u(i);
        u_array2d[row][1] = u(i + 1);
        u_matrix2d.row(row)(0) = u(i);
        u_matrix2d.row(row)(1) = u(i + 1);
        nodes_puls_displ_2d[row][0] = nodes.row(row)(0) + u(i);
        nodes_puls_displ_2d[row][1] = nodes.row(row)(1) + u(i + 1);
    }



    Eigen::MatrixXd u_x(fem_params.num_nodes, 1 );
    u_x = u(Eigen::seq(0, 2*fem_params.num_nodes-2, 2));
    Eigen::MatrixXd u_y(fem_params.num_nodes, 1 );
    u_y = u(Eigen::seq(1,  2*fem_params.num_nodes, 2));
    Eigen::MatrixXd displacement(fem_params.num_nodes,  fem_params.num_nodes);




    auto cpp_Nodes = eigenArrayToCppArray(nodes);
    auto cpp_Connectivity = eigenArrayToCppArray(connectivity);


    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (const auto &point : cpp_Nodes) {
        points->InsertNextPoint(point[0], point[1], 0.0);
    }

    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    grid->Allocate();



    // Define points and cells for deformed mesh
    vtkSmartPointer<vtkPoints> deformedPoints = vtkSmartPointer<vtkPoints>::New();
    // Add points to deformedPoints
    for (const auto &point : nodes_puls_displ_2d) {
        deformedPoints->InsertNextPoint(point[0], point[1], 0.0);
    }


    // Create deformed mesh
    vtkSmartPointer<vtkUnstructuredGrid> deformed_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    deformed_grid->SetPoints(deformedPoints);
    deformed_grid->Allocate();

    for (const auto &el : cpp_Connectivity) {
        vtkSmartPointer<vtkQuad> cell = vtkSmartPointer<vtkQuad>::New();
        for (size_t i = 0; i < el.size(); ++i) {
            cell->GetPointIds()->InsertId(i, el[i]);
        }
        grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
    }

    for (const auto &el : cpp_Connectivity) {
        vtkSmartPointer<vtkQuad> cell_def = vtkSmartPointer<vtkQuad>::New();
        for (size_t i = 0; i < el.size(); ++i) {
            cell_def->GetPointIds()->InsertId(i, el[i]);
        }
        deformed_grid->InsertNextCell(cell_def->GetCellType(), cell_def->GetPointIds());
    }


    vtkSmartPointer<vtkDoubleArray> normalData = vtkSmartPointer<vtkDoubleArray>::New();
    normalData->SetNumberOfComponents(2);
    normalData->SetName("displacement");
    for (const auto &normal : u_array2d) {
        normalData->InsertNextTuple2(normal[0], normal[1]);
    }
    grid->GetPointData()->AddArray(normalData);
    deformed_grid->GetPointData()->AddArray(normalData);


    // Create writer and set file name
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName("../ParaView/grid_1.vtu");

    // Write undeformed mesh to file
    writer->SetInputData(grid);
    writer->Write();

    writer->SetFileName("../ParaView/grid_2.vtu");
    // Write deformed mesh to file
    writer->SetInputData(deformed_grid);
    writer->Write();



    fstream my_file;
    my_file.open("../ParaView/tension_rod.pvd", ios::out);
    if (!my_file) {
        cout << ".pvd File not created!";
    }
    else {
        cout << ".pvd File created successfully!";
        my_file << "<?xml version=\"1.0\"?>\n";
        my_file << "<VTKFile type=\"Collection\" version=\"0.1\" "
                   "compressor=\"vtkZLibDataCompressor\" byte_order=\"BigEndian\">\n";
        my_file << "\t <Collection>\n";
        my_file << "\t \t<DataSet timestep=\"0\" group=\"\" part=\"0\" file=\"grid_1.vtu\" />\n";
        my_file << "\t \t<DataSet timestep=\"1\" group=\"\" part=\"0\" file=\"grid_2.vtu\" />\n";
        my_file << "\t </Collection>\n";
        my_file << "</VTKFile>";

        my_file.close();
    }

//    plt::plot({1,3,2,4});
//    plt::show();

    Eigen::MatrixXd ux = Eigen::Map<Eigen::MatrixXd>(u.data(), fem_params.mesh_ny, fem_params.mesh_nx);
    Eigen::MatrixXd uy = Eigen::Map<Eigen::MatrixXd>(u.data() + 1, fem_params.mesh_ny, fem_params.mesh_nx);





    std::ofstream data_file("../data/data.dat");
    if (data_file.is_open()) {
        for (int r = 0; r < connectivity.rows(); r++) {
            Eigen::MatrixXd x_Ie = nodes(connectivity.row(r), Eigen::all);
            x_Ie.conservativeResize(5, 2);
            x_Ie.row(x_Ie.rows()-1) = x_Ie.row(0);


            Eigen::MatrixXd u_Ie = u_matrix2d(connectivity.row(r), Eigen::all);
            u_Ie.conservativeResize(5, 2);
            u_Ie.row(x_Ie.rows()-1) = u_Ie.row(0);
            for (int i = 0; i < 5; i++) {
                data_file << x_Ie.row(i)(0) << " " << x_Ie.row(i)(1) << " "
                          << u_Ie.row(i)(0) << " "<<  u_Ie.row(i)(1) << " " << std::endl;
            }
            data_file << std::endl;
        }
        data_file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }




    std::ofstream data_file2("../data/data2.dat");
    if (data_file2.is_open()) {
        for (int r = 0; r < nodes.rows(); r++) {

            data_file2 << nodes.row(r)(0) << " " <<  nodes.row(r)(1)<< " "<< u_matrix2d.row(r)(1)  << std::endl;
        }
        data_file2.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }

    std::system("(cd ../gnuplot/ && gnuplot plot.gnuplot)");
    std::system("(cd ../gnuplot/ && gnuplot plot2.gnuplot)");
    return 0;

}


