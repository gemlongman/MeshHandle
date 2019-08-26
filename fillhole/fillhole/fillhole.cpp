// fillhole.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"

//#include "HoleFiller.h"
#define  _CRT_SECURE_NO_WARNINGS
#include "InheritOMesh.h"
#include "HoleMesh.h"

#include <iostream>
#include <cstdio>
#include <vector>

using namespace std;

///////////////////////

string caseBall = "G:\\gyGit\\data\\hole\\sphere_hole.obj";
string caseBunnyOneSmall = "G:\\gyGit\\data\\hole\\bunny_a_small_hole.obj";// 
string caseBunnyOneBig = "G:\\gyGit\\data\\hole\\bunny_a_big_hole.obj";
string caseBunny = "G:\\gyGit\\data\\hole\\bunny_with_holes.obj";

string outFile0 = "G:\\gyGit\\data\\hole\\out_temp0.obj";
string outFile1 = "G:\\gyGit\\data\\hole\\out_temp1.obj";
string outFile2 = "G:\\gyGit\\data\\hole\\out_temp2.obj";
string outFile3 = "G:\\gyGit\\data\\hole\\out_temp3.obj";

void ot() {
	vector<int> arr={ 1,2,3,4,5 };
	vector<int>::iterator pre, cur, next;
	pre = prev(arr.begin());
	cur = arr.begin();
	next = std::next(arr.end());
	cout << " " << *pre << " " << *cur << " " << *next;
	GetPrevNextIter(arr, pre, cur, next);
	cout << " " << *pre << " " << *cur << " "<< *next;

	cur = prev(arr.end());
	GetPrevNextIter(arr, pre, cur, next);
	cout << " " << *pre << " " << *cur << " " << *next;
}

void useMyHoleFiller()
{
	MeshHandle meshh_handle;
	meshh_handle.Init(caseBunny);
	meshh_handle.FindHoles();
	meshh_handle.FillHoles();
	meshh_handle.WriteMeshFile(outFile3);
}

void useMyHoleMesh()
{
	HoleMeshHandle meshh_handle;
	meshh_handle.Init(caseBall);
	meshh_handle.ExtractHolesRegin(1,3.0);//0 is boundary
	meshh_handle.FillHoles(0);
	meshh_handle.SmoothHoles(1);
	meshh_handle.MergerBack(true);
	meshh_handle.WriteMeshFile(outFile3);
}

void useHoleFiller()
{
	//////////////////////////////////
	//HoleFiller hFiller;
	//hFiller.Run("G:\\gyGit\\data\\hole\\off\\bunny_with_holes.off", "out.off");

	//MeshHandle mesh_save_as;
	//mesh_save_as.Save_file_as("out.off", outFile2);
}

int main(int, char ** argv)
{
	//ot();


	useMyHoleMesh();

	//useMyHoleFiller();

	//useHoleFiller();
	//getchar();
	return 0;
}

