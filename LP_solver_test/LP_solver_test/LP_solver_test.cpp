#include"stdafx.h"//プリコンパイル済みヘッダ
#include<fstream>
#include<iostream>
#include<math.h>
#include"LP_SOLVER01.h"
using namespace std;

int main(){
	int i, j, k;

	double objective[2] = { -2.0, -3.0 };
	double leftSide[6] = { 1.0, 2.0, 1.0, 1.0, 3.0, 1.0 };
	double rightSide[3] = { 14.0, 8.0, 18.0 };
	double optimumSolution[2];


	LP_SOLVER01(2, 3, 0, objective, leftSide, rightSide, optimumSolution);
	cout << "x1 = " << optimumSolution[0] << "  x2 = " << optimumSolution[1] << endl;
	cout << "最大値は " << objective[0] * optimumSolution[0] + objective[1] * optimumSolution[1] << endl;

	return 0;
}