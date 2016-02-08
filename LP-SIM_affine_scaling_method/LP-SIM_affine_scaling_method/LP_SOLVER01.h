#include"stdafx.h"
#include<stdlib.h>
#include<iostream>
#include<math.h>
using namespace std;

#define DEGREE 5

/**********　行列の掛け算　**********/
int calcMatrixMultiplication(
	int leftColumn,		// 縦（列、column）；左から掛ける行列の縦の長さ
	int leftRow,		// 横（行、row）   ；左から掛ける行列の横の長さ
	int rightRow,		// 横（行、row）   ；右から掛ける行列の横の長さ
	double *left,		// 左から掛ける行列
	double *right,		// 右から掛ける行列
	double solution[])	// 結果を格納する配列
{
	int i, j, k;
	double term;

	//掛け算
	for (i = 0; i < leftColumn; i++) {
		for (j = 0; j < rightRow; j++) {
			term = 0.0;
			for (k = 0; k < leftRow; k++){
				term = term + left[i*leftRow + k] * right[(rightRow* k) + j];
			}
			solution[rightRow *i + j] = term;
		}
	}

	return 0;
}

/**********　方程式の計算　**********/
double solvEquation(
	int size,			// 変数の数
	double *A,			// 正方行列（連立方程式の行列Aの部分）
	double *B,			// 縦ベクトル（連立方程式のBの部分）
	double *solution)	// 結果を返す配列
{
	int i, j, k, n, m;
	double *augmentedMatrix;
	double pivo, mult;

	augmentedMatrix = (double *)malloc(sizeof(double)*(size + 1)*size);
	if (augmentedMatrix == NULL) exit(0);

	//拡大係数行列生成
	for (i = 0; i < size; ++i){
		for (j = 0; j < size; ++j){
			augmentedMatrix[i*(size + 1) + j] = A[i*size + j];
		}
	}
	for (i = 0; i < size; i++){
		augmentedMatrix[(size + 1)*(i + 1) - 1] = B[i];
	}

	//ガウスの消去法の計算
	for (i = 0; i < size; i++)
	{
		// 正規化(最初の項の係数を１にする)
		pivo = augmentedMatrix[(size + 1)*i + i];
		for (j = 0; j < size + 1; ++j)
		{
			augmentedMatrix[(size + 1)*i + j] = (1 / pivo) * augmentedMatrix[(size + 1)*i + j];
		}

		// 逆三角形の行列を生成する
		// 現在の行より下の行のi列目が0になるように引き算をする
		for (k = i + 1; k < size; ++k)
		{
			mult = augmentedMatrix[(size + 1)*k + i];
			for (n = i; n < size + 1; ++n)
			{
				augmentedMatrix[(size + 1)*k + n] -= mult * augmentedMatrix[(size + 1)*i + n];
			}
		}
	}

	// 後進代入(各変数が独立した形になるように変形する)
	for (i = size - 1; i > 0; --i)
	{
		for (k = i - 1; k >= 0; --k)
		{
			mult = augmentedMatrix[(size + 1)*k + i];
			for (n = i; n < size + 1; ++n)
			{
				augmentedMatrix[(size + 1)*k + n] -= mult * augmentedMatrix[(size + 1)*i + n];
			}
		}
	}

	//計算結果の取り出し
	for (i = 0; i < size; i++){
		solution[i] = augmentedMatrix[(size + 1) *(i + 1) - 1];
	}

	return 0;
}

/**********　メイン関数　**********/
int LP_SOLVER01(int size, double objective[], double leftSide[], double rightSide[], double optimumSolution[]){

	int i, j, k;//カウンタ

	double *A;
	A = (double *)malloc(sizeof(double)*(size + size)*(size + size));

	for (i = 0; i < size; ++i){
		for (j = 0; j < 2 * size; ++j){
			if (j < size){ A[(2 * size)*i + j] = leftSide[size*i + j]; }
			else if (2 * size*i + j - size *(1 + 2 * i) == i){ A[2 * size*i + j] = 1.0; }
			else{ A[2 * size*i + j] = 0.0; }
		}
	}

	for (i = 0; i < size; ++i){
		rightSide[i] = rightSide[i] - 1.0;
	}
	double *xk;
	xk = (double *)malloc(sizeof(double)*(size + size));

	solvEquation(size, leftSide, rightSide, xk);///////////
	for (i = 0; i < size; ++i){
		xk[size + i] = 1.0;
	}


	for (i = 0; i < size; i++){
		xk[size + i] = 1.0;
	}

	double A_tra[8];	//転置行列
	for (i = 0; i < 2 * size; i++){//転置行列の計算
		for (j = 0; j < size; j++){
			A_tra[j + 2 * i] = A[i + j * 4];
		}
	}

	double *C;
	C = (double *)malloc(sizeof(double)*(size + size));
	for (i = 0; i < 2 * size; i++){//転置行列の計算
		if (i<size){
			C[i] = objective[i];
		}
		else{
			C[i] = 0.0;
		}

	}

	//結果や途中家計算を格納するための配列
	double AX2At[4];//方程式の左辺　
	double AX2c[2];//方程式の右辺

	double diagMatrix[16];//対角行列を格納する
	double diagMatrixSqua[16];//対角行列を2乗して格納する

	double solu1[8];//途中計算
	double solu2[2];//途中計算
	double solu3[4];//途中計算
	double Z[4];

	double dx[4];//更新則
	double z_nol = 0;
	double z1 = 0;

	double upsilon = 1.0e-15;//終了条件　ε
	double alpha = 0.8;//更新則に掛かるゲイン

	while (1){
		//左辺
		//対角行列生成
		for (i = 0; i < 4; ++i){
			for (j = 0; j < 4; ++j){
				diagMatrix[i * 4 + j] = 0.0;
				if (j == i)diagMatrix[i * 4 + j] = xk[i];
			}
		}
		calcMatrixMultiplication(4, 4, 4, diagMatrix, diagMatrix, diagMatrixSqua);//対角行列の平方
		calcMatrixMultiplication(2, 4, 4, A, diagMatrixSqua, solu1);//対角行列の平方
		calcMatrixMultiplication(2, 4, 2, solu1, A_tra, AX2At);

		//右辺
		calcMatrixMultiplication(2, 4, 1, solu1, C, AX2c);
		solvEquation(2, AX2At, AX2c, solu2);
		calcMatrixMultiplication(4, 2, 1, A_tra, solu2, solu3);
		for (i = 0; i < 4; ++i){
			Z[i] = C[i] - solu3[i];
		}

		for (i = 0; i < 4; ++i){
			z_nol += (diagMatrix[4 * i + i] * Z[i])*(diagMatrix[4 * i + i] * Z[i]);
		}
		z_nol = sqrt(z_nol);

		//変化量を計算する
		for (i = 0; i < 4; ++i){
			dx[i] = diagMatrixSqua[4 * i + i] * Z[i] / z_nol;
		}

		//ｘの値を更新する
		for (i = 0; i < 4; ++i){
			xk[i] += alpha * dx[i];
		}

		//判定式の計算
		for (i = 0; i < 4; ++i){
			z1 +=  dx[i] * dx[i];
		}
		if (upsilon>z1)break;
		z1 = 0.0;
	}

	//計算結果を返す（返す用の配列に渡す）
	for (i = 0; i < size; ++i){
		optimumSolution[i] = xk[i];
	}


	return 0;
}
