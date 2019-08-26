#include "MatrixHandle.h"

#include <cstdio>
#include <iostream>



using namespace Eigen;
using namespace std;

#define Tolerance 0.000001

///////////////////////////////////////////////////////////////////////////////// test
#include <ctime>
void test() {
	MatrixXd m = MatrixXd::Random(3, 3);              //随机生成3*3的double型矩阵
	m = (m + MatrixXd::Constant(3, 3, 1.2)) * 50;      //MatrixXd::Constant(3,3,1.2)表示生成3*3的double型矩阵，该矩阵所有元素均为1.2
	cout << "m =" << endl << m << endl;
	VectorXd v(3);        // 定义v为3*1的double型向量
	v << 1, 2, 3;         // 向量赋值
	cout << "m * v =" << endl << m * v << endl;

}

void axb() {
	// 解方程
	// 我们求解 A * x = b 这个方程
	// 直接求逆自然是最直接的，但是求逆运算量大

	Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > A1;
	A1 = Eigen::MatrixXd::Random(128, 128);

	Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > b1;
	b1 = Eigen::MatrixXd::Random(128, 1);

	clock_t time_stt = clock(); // 计时
	// 直接求逆
	Eigen::Matrix<double, 128, 1> x = A1.inverse()*b1;
	cout << "time use in normal inverse is " << 1000 * (clock() - time_stt) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << x << endl;

	// QR分解colPivHouseholderQr()
	time_stt = clock();
	x = A1.colPivHouseholderQr().solve(b1);
	cout << "time use in Qr decomposition is " << 1000 * (clock() - time_stt) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << x << endl;

	//QR分解fullPivHouseholderQr()
	time_stt = clock();
	x = A1.fullPivHouseholderQr().solve(b1);
	cout << "time use in Qr decomposition is " << 1000 * (clock() - time_stt) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << x << endl;

	//lu分解 partialPivLu()
	time_stt = clock();
	x = A1.partialPivLu().solve(b1);
	cout << "time use in lu decomposition is " << 1000 * (clock() - time_stt) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << x << endl;

	//lu分解（fullPivLu()
	time_stt = clock();
	x = A1.fullPivLu().solve(b1);
	cout << "time use in lu decomposition is " << 1000 * (clock() - time_stt) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << x << endl;

}

//solve A*x=b return x
int solverLSCG(const Eigen::SparseMatrix<float> & A1_sparse ,const Eigen::VectorXf & b1_sparse,Eigen::VectorXf & x1_sparse) {
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float> > Solver_sparse;
	Solver_sparse.setTolerance(Tolerance);// 设置迭代精度
	
	Solver_sparse.compute(A1_sparse);//x1_sparse 即为解
	
	if (Solver_sparse.info() !=Eigen::ComputationInfo::Success) {
		cerr << "solver error!" << endl;
		return Solver_sparse.info();
	}
	x1_sparse = Solver_sparse.solve(b1_sparse);
	//debug cout << x1_sparse << endl;
	return 0;
}

//error? solve A*x=b return x
int solverLU(const Eigen::SparseMatrix<float> & A1_sparse, const Eigen::VectorXf & b1_sparse, Eigen::VectorXf & x1_sparse) {
	Eigen::SparseLU<Eigen::SparseMatrix<float> , COLAMDOrdering<int>> Solver_sparse;
	Solver_sparse.analyzePattern(A1_sparse);
	Solver_sparse.factorize(A1_sparse);

	
	if (Solver_sparse.info() != Eigen::ComputationInfo::Success) {
		cerr << "solver error!" << endl;
		return Solver_sparse.info();
	}
	x1_sparse = Solver_sparse.solve(b1_sparse);
	//debug 
	cout << x1_sparse << endl;
	return 0;
}

void testAxb() {
	int n = 10;

	std::vector<Eigen::Triplet<float>> tripletlist;
	Eigen::VectorXf b1_sparse(n);
	Eigen::VectorXf x1_sparse(n);
	
	//test value:
	for (int i = 0; i < n; i++)
	{
		int j = i;
		tripletlist.push_back(Eigen::Triplet<float>(i, j, j));//对角线为1
		b1_sparse[i] = i * j;//b
	}

	Eigen::SparseMatrix<float> A1_sparse(n, n);
	A1_sparse.setFromTriplets(tripletlist.begin(), tripletlist.end());
	A1_sparse.makeCompressed();// 压缩优化矩阵

	int ret = solverLSCG(A1_sparse, b1_sparse, x1_sparse);
	//debug 
	cout << x1_sparse << endl;
}

/////////////////////////////////////////////////////////////////////////////////end test

//MatrixHandle::MatrixHandle()
//{
//	mesh = NULL;
//}

MatrixHandle::~MatrixHandle()
{
	mesh = NULL;//!!! don't delete
}

int MatrixHandle::calulateWeight(const vector<OpenMesh::VertexHandle> & neighborVs, vector<float> & weights) {
	//unify
	for ( int i = 0; i < neighborVs.size(); i++) {
		weights[i] = 1.0 / neighborVs.size();
	}
	//todo: face center? according distance or area or 双端滤波? 
	return 0;
}

void MatrixHandle::addAPointWeights(int row,int col,float w) {
	for (int i = 0; i < 3;i++) {
		TripletList.push_back(Eigen::Triplet<float>(row+i, col+i, w));
	}
}

//this error
int MatrixHandle::FillMat()
{
	int fixRing = 1;
	int row, col;
	float value;
	OpenMesh::VertexHandle curVH;

	//sure that index since 0 to dim-1
	//unknow x: w1*x1 +...+ wn*pn =xi so:
	//know p:1 * pi=pi w_cur=-1,b=0 or w_cur=0,b=-xi
	for (int v_i = 0; v_i < N; v_i++) {
		curVH = mesh->vertex_handle(v_i);	//v_i==curVH.idx()		
		row = v_i * 3;
		if (mesh->data(curVH).type < fixRing) {// -1 is new point// 0 is boundary // >0 are out ring of boundary

			//calucate a line weights
			vector<OpenMesh::VertexHandle> neighborVs;
			for (auto vv = mesh->vv_cwbegin(curVH); vv != mesh->vv_cwend(curVH); vv++) {// 1 ring neighbor 
				neighborVs.push_back(*vv);
			}
			vector<float> weights(neighborVs.size());
			calulateWeight(neighborVs, weights);
			//debug
			//cout << "calulateWeight" << neighborVs.size() << endl;
			//for (int j = 0; j< neighborVs.size(); j++) {
			//	cout << weights[j] << endl;
			//}
			//fill a line weight
			addAPointWeights(row, row, 1.0);
			int i = 0;
			for (vector<OpenMesh::VertexHandle>::iterator v = neighborVs.begin(); v != neighborVs.end(); v++, i++) {
				col = v->idx() * 3;
				addAPointWeights(row, col, -weights[i]);
			}
			//b=0
		}
		else {	// 1*p=p
			addAPointWeights(row, row, 1.0);
			//fill a ponit into b[vi]
			for (int i = 0; i < 3; i++) {
				b_sparse[row + i] = mesh->point(curVH)[i];
			}
		}

	}
	A_sparse.setFromTriplets(TripletList.begin(), TripletList.end());
	A_sparse.makeCompressed();// 压缩优化矩阵
	//debug 
	// cout << "A_sparse" << A_sparse << endl;
	//cout << "b_sparse" << b_sparse << endl;
	return 0;
}

// -1 is new point// 0 is boundary // >0 are out ring of boundary
int MatrixHandle::FillMatWithFix(int fixRing)
{
	int row, col;
	float value;
	OpenMesh::VertexHandle curVH;

	//sure that index since 0 to dim-1
	//unknow x: w1*x1 +...+ wn*pn =xi so: w1*x1 +...+wj*xj - xi = wk*pk +...+ wn*pn so: (w1 w2 ..)x=b
	//know p:1 * pi=pi w_cur=-1,b=0 or w_cur=0,b=-xi
	for (int v_i = 0; v_i < N; v_i++) {
		curVH = mesh->vertex_handle(v_i);	//v_i==curVH.idx()		
		row = v_i * 3;
		if (mesh->data(curVH).type < fixRing) {// -1 is new point// 0 is boundary // >0 are out ring of boundary
			//todo is it three weight of a point?
			//calucate a line weights
			vector<OpenMesh::VertexHandle> neighborVs;
			for (auto vv = mesh->vv_cwbegin(curVH); vv != mesh->vv_cwend(curVH); vv++) {// 1 ring neighbor 
				neighborVs.push_back(*vv);
			}
			vector<float> weights(neighborVs.size());
			calulateWeight(neighborVs, weights);
			//debug
			//cout << "calulateWeight" << neighborVs.size() << endl;
			//for (int j = 0; j< neighborVs.size(); j++) {
			//	cout << weights[j] << endl;
			//}
			//fill a line weight
			addAPointWeights(row, row, 1);
			int i = 0;
			for (vector<OpenMesh::VertexHandle>::iterator v = neighborVs.begin(); v != neighborVs.end(); v++, i++) {
				if (mesh->data(*v).type < fixRing) {//unknow p
					col = v->idx() * 3;
					addAPointWeights(row, col, -weights[i]);
				}
				else {//know p
					//fill w * ponit into b[vi]
					for (int i = 0; i < 3; i++) {
						b_sparse[row + i] += weights[i] * mesh->point(curVH)[i];
					}
				}
			}
			
		}
		else {	// 1*p=p
			col = row;
			addAPointWeights(row, col, 1.0);
			//fill a ponit into b[vi]
			for (int i = 0; i < 3; i++) {
				b_sparse[row + i] = mesh->point(curVH)[i];
			}
		}

	}
	A_sparse.setFromTriplets(TripletList.begin(), TripletList.end());
	A_sparse.makeCompressed();// 压缩优化矩阵
	//debug 
	// cout << "A_sparse" << A_sparse << endl;
	//cout << "b_sparse" << b_sparse << endl;
	return 0;
}

int MatrixHandle::SolveLSCG()
{
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float> > Solver_sparse;
	Solver_sparse.setTolerance(Tolerance);// 设置迭代精度

	Solver_sparse.compute(A_sparse);

	if (Solver_sparse.info() != Eigen::ComputationInfo::Success) {
		cerr << "solver error!" << endl;
		return Solver_sparse.info();
	}
	x_sparse = Solver_sparse.solve(b_sparse);
	//debug 
	//cout<<"x_sparse" << x_sparse << endl;
	return 0;
}

int MatrixHandle::SolveLU()
{
	Eigen::SparseLU<Eigen::SparseMatrix<float>, COLAMDOrdering<int>> Solver_sparse;
	Solver_sparse.analyzePattern(A_sparse);
	Solver_sparse.factorize(A_sparse);

	if (Solver_sparse.info() != Eigen::ComputationInfo::Success) {
		cerr << "solver error!" << endl;
		return Solver_sparse.info();
	}
	x_sparse = Solver_sparse.solve(b_sparse);
	//debug 
	//cout<<"x_sparse" << x_sparse << endl;
	return 0;
}

int MatrixHandle::UpdateMesh()
{
	OpenMesh::VertexHandle curVH;
	OpenMesh::Vec3f point;
	for (int v_i = 0; v_i < N; v_i++) {
		curVH = mesh->vertex_handle(v_i);
		if (mesh->data(curVH).type>0) {//skip boundary
			continue;
		}
		Eigen::VectorXf temp = x_sparse.block<3, 1>(v_i * 3, 0);
		mesh->set_point(curVH, { temp[0], temp[1], temp[2] });

		//float * arraya;
		//arraya = x_sparse.block<3, 1>(v_i * 3, 0).data();
		//point = static_cast<OpenMesh::Vec3f>(arraya);
		//mesh->set_point(curVH, point);
	}
	return 0;
}

MatrixHandle::MatrixHandle(OMTriMesh & iMesh)
{
	mesh = & iMesh;
	N = mesh->n_vertices();
	Dim = 3 * mesh->n_vertices();
	A_sparse.resize(Dim, Dim);
	x_sparse.resize(Dim);
	b_sparse.resize(Dim);
}

int MatrixHandle::SmoothMesh(int fixPointType)
{
	int ret = 0;
	ret |= FillMat();
	//ret |= FillMatWithFix(fixPointType);//固定外层 <=neighborRings, 
	ret |= SolveLU();//ret |= SolveLSCG();
	ret |= UpdateMesh();
	return ret;
}
