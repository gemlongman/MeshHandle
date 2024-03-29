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

///ot
/*
*/
const int MaxN = 10;

int mapNum[MaxN][MaxN];
bool visited[MaxN];
int numArr[MaxN];

void init() {
	for (int i = 0; i < MaxN; i++) {
		numArr[i] = i;
		visited[i] = false;
	}
	mapNum[1][2] = 1;    mapNum[2][1] = 1;
	mapNum[1][5] = 1;    mapNum[5][1] = 1;
	mapNum[2][3] = 1;    mapNum[3][2] = 1;
	mapNum[2][4] = 1;    mapNum[4][2] = 1;
	mapNum[5][6] = 1;    mapNum[6][5] = 1;
	mapNum[6][9] = 1;    mapNum[9][6] = 1;
	mapNum[5][7] = 1;    mapNum[7][5] = 1;
	mapNum[5][8] = 1;    mapNum[8][5] = 1;
}

void printmap(int start) {
	cout << numArr[start] << endl;
	visited[start] = true;

	vector<int> neibors;
	int top = 0;
	neibors.push_back(start);
	//while(!neibors.empty()){

	while (top < neibors.size()) {
		int cur = neibors[top];
		//neibors.front();
		//neibors.pop_front();
		top++;
		for (int i = 0; i < MaxN; i++) {
			if (!visited[i] && mapNum[cur][i]) {
				cout << numArr[i] << endl;
				neibors.push_back(i);
				visited[i] = true;
			}
		}
	}

}


void ot() {
	init();
	printmap(1);
	return;
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
	ot();


	//useMyHoleMesh();

	//useMyHoleFiller();

	//useHoleFiller();
	//getchar();
	return 0;
}

