#include"stdafx.h"//プリコンパイル済みヘッダ
#include<fstream>
#include<iostream>
#include<math.h>
#include"LP_SOLVER04.h"
using namespace std;

int main(){
	int i, j, k;

	double objective[3] = { 1.0, 2.0 ,3.0 };
	double leftSide[9] = { 1.0, 0.0,1.0 , 2.0, 1.0,2.0, 3.0, 1.0, 2.0 };
	double rightSide[3] = { 2.0, 5.0, 6.0 };
	double optimumSolution[3];
	double result;

	result = LP_SOLVER04(sizeof objective / sizeof(double), sizeof rightSide / sizeof(double), 0, objective, leftSide, rightSide, optimumSolution);
	for (i = 0; i < sizeof rightSide / sizeof(double); i++){
		cout << "x"<<i+1 <<" = " << optimumSolution[i] << " ," ;
	}
	cout << endl;
	cout << "最適値は " << result << endl;

	return 0;
}