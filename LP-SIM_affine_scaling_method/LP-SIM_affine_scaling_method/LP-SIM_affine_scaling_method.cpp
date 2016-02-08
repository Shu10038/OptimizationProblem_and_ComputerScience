#include"stdafx.h"
#include<fstream>
#include<iostream>
#include<math.h>
#include"LP_SOLVER01.h"
using namespace std;

int main(){
	int i, j, k;
	
	int size = 2;
	double objective[2] = { 1.0,1.0 };
	double leftSide[4] = { 2.0, 1.0, 1.0, 3.0 };
	double rightSide[2] = { 4.0, 5.0 };
	double optimumSolution[2];

	LP_SOLVER01(size, objective, leftSide, rightSide, optimumSolution);	
	cout << "x1 = " << optimumSolution[0] << "  x2 = " << optimumSolution[1] << endl;
	cout << "Å¬’l‚Í " << optimumSolution[0] + optimumSolution[1] << endl;

	return 0;
}