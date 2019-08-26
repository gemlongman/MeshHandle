#pragma once
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "InheritOMesh.h"

class MatrixHandle {

public:
	//MatrixHandle();
	MatrixHandle(OMTriMesh & iMesh);//todo
	~MatrixHandle();

	int FillMat();
	int FillMatWithFix(int fixRing = 1);
	int SolveLSCG();
	int SolveLU();
	int UpdateMesh();
	int SmoothMesh(int fixPointType);
private:
	int Dim, N;//dim=3*n
	//int DimB, Nb;
	OMTriMesh * mesh;
	std::vector<Eigen::Triplet<float>> TripletList;
	Eigen::VectorXf b_sparse;
	Eigen::VectorXf x_sparse;
	Eigen::SparseMatrix<float> A_sparse;

	int calulateWeight(const std::vector<OpenMesh::VertexHandle> & neighborVs, std::vector<float> & weights);
	void addAPointWeights(int row, int col, float w);
};