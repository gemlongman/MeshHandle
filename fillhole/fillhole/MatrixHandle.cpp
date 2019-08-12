#include "MatrixHandle.h"
#include <cstdio>
#include <iostream>

using namespace Eigen;
using namespace std;

void test() {
	MatrixXd m = MatrixXd::Random(3, 3);              //随机生成3*3的double型矩阵
	m = (m + MatrixXd::Constant(3, 3, 1.2)) * 50;      //MatrixXd::Constant(3,3,1.2)表示生成3*3的double型矩阵，该矩阵所有元素均为1.2
	cout << "m =" << endl << m << endl;
	VectorXd v(3);        // 定义v为3*1的double型向量
	v << 1, 2, 3;         // 向量赋值
	cout << "m * v =" << endl << m * v << endl;

}
