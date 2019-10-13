#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <memory>

namespace ProjDyn {

	class Constraint;
	typedef std::shared_ptr<Constraint> ConstraintPtr;


	typedef double Scalar;										// A scalar type, double or float, as defined
	typedef size_t Index;										// Type for array and vector sizes and indices.

	//Dense
	template < int Rows, int Cols, int Options = (Eigen::ColMajor) >
	using MatrixT = Eigen::Matrix<Scalar, Rows, Cols, Options>; // A typedef of the dense matrix of Eigen.

	typedef MatrixT<Eigen::Dynamic, 3> Positions;				// A n by 3 matrix of scalars, often used for positions or m_velocities
	typedef MatrixT<Eigen::Dynamic, 1> Vector;					// Vector of scalars
	typedef MatrixT<3, 1> Vector3;								// Vector of 3 scalars
	typedef MatrixT<2, 1> Vector2;								// Vector of 2 scalars
	typedef MatrixT<Eigen::Dynamic, Eigen::Dynamic> Matrix;		// Arbitrary size Matrix

	typedef Eigen::Matrix<Index, Eigen::Dynamic, 3> Triangles;	// List of three indices per row, used for triangles
	typedef Eigen::Matrix<Index, 1, 3> Triangle;				// A single triangle (i.e. four indices).
	typedef Eigen::Matrix<Index, Eigen::Dynamic, 4>				// List of four indices per row, used for tetrahedrons
		Tetrahedrons;
	typedef Eigen::Matrix<Index, 1, 4> Tetrahedron;				// A single tetrahedron (i.e. four indices).

	//Sparse
	template<int Options = Eigen::ColMajor>
	using SparseMatrixT = Eigen::SparseMatrix<Scalar, Options>;	// A typedef of the sparse matrix of Eigen.
	typedef SparseMatrixT<> SparseMatrix;						// A column-major sparse matrix.
	typedef SparseMatrixT<Eigen::RowMajor> SparseMatrixRM;		// A row-major sparse matrix.
	typedef Eigen::Triplet<Scalar> Triplet;						// A triplet, used in the sparse triplet representation for matrices.

	//Quaternions
	using Quaternion = Eigen::Quaternionf;
    template < int Rows, int Cols>
	using ArrayQ = Eigen::Array<Quaternion, Rows, Cols>;
    typedef ArrayQ<Eigen::Dynamic, 1> Orientations;

    template < int Rows, int Cols>
	using FlatM = Eigen::Array<Scalar, Rows, Cols>;
    typedef FlatM<Eigen::Dynamic, 1> FlatPos;

}
