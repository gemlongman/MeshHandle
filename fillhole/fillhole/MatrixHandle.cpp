#include "MatrixHandle.h"
#include <cstdio>
#include <iostream>

using namespace Eigen;
using namespace std;

void test() {
	MatrixXd m = MatrixXd::Random(3, 3);              //�������3*3��double�;���
	m = (m + MatrixXd::Constant(3, 3, 1.2)) * 50;      //MatrixXd::Constant(3,3,1.2)��ʾ����3*3��double�;��󣬸þ�������Ԫ�ؾ�Ϊ1.2
	cout << "m =" << endl << m << endl;
	VectorXd v(3);        // ����vΪ3*1��double������
	v << 1, 2, 3;         // ������ֵ
	cout << "m * v =" << endl << m * v << endl;

}
