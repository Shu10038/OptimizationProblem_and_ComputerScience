#include"stdafx.h"
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include<math.h>
using namespace std;


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
	int i, j, k, n;
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
int main(){
	ofstream fout("output.txt");
	if (!fout){
		cout << "This file can not open!!";
		return 1;
	}

	int i, j;//カウンタ

	//初期値の設定

	double A[8] = { 2.0, 1.0, 1.0, 0.0,
		1.0, 3.0, 0.0, 1.0 };

	double A_tra[8];	//転置行列
	for (i = 0; i < 4; i++){//転置行列の計算
		for (j = 0; j < 2; j++){
			A_tra[j + 2 * i] = A[i + j * 4];
		}
	}
	double C[4] = { 4.0, 5.0, 0.0, 0.0 };//制約条件の式の右辺
	double xk[4] = { 1.0, 1.0, 1.0, 1.0 };//初期値（初期内点は計算して出した）

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

	double upsilon = 1.0e-20;//終了条件　ε
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

		for (i = 0; i < 4; ++i){
			dx[i] = diagMatrixSqua[4 * i + i] * Z[i] / z_nol;
		}

		//k+1の点の更新
		for (i = 0; i < 4; ++i){
			xk[i] += +alpha * dx[i];
		}

		//判定条件の計算
		for (i = 0; i < 4; ++i){
			z1 += z1 + dx[i] * dx[i];
		}
		if (upsilon>z1)break;
		z1 = 0.0;
	}
	cout << "x1 = " << xk[0] << "  X2 = " << xk[1] << endl;
	cout << "最大値は " << xk[0] + xk[1] << endl;

	return 0;
}
