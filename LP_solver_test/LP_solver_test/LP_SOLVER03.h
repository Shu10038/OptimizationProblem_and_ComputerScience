/*
アフィン変換法を用いた線形計画法のソルバ
This is a solver of linear programming .It's used the affine scaling method.
’16/01/29　S.OHSAWA

’16/02/13 修正
*/

#include"stdafx.h"
#include<stdlib.h>
#include<iostream>
#include<math.h>
using namespace std;


/******************************　
行列の掛け算　
******************************/
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

/******************************
方程式の計算　
******************************/
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

/******************************　　
線形計画法のソルバ
******************************/
double LP_SOLVER03(
	int numbVar,			//変数の数 
	int numbLim,			//制約条件の数
	int numbEqu,			//等式条件
	double objective[],		//目的関数の係数　例)3x_1+2x_2　だったら{3.0,2.0}
	double leftSide[],		//制約条件の左側．例)2x_1+4x_2≦4 と x_1+3x_2≦5　だったら{2.0,4.0,1.0,3.0}
	double rightSide[],		//制約条件の右側．
	double optimumSolution[])
{

	int i, j, k;//カウンタ


	double *A;
	A = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbLim));

	//■STEP1；標準化
	//スラック変数を挿入して標準化する
	for (i = 0; i < numbLim; ++i){
		for (j = 0; j < numbVar + numbLim; ++j){
			if (j < numbVar){ A[(numbVar + numbLim)*i + j] = leftSide[numbVar*i + j]; }
			else if (j - numbVar == i){ A[(numbVar + numbLim)*i + j] = 1.0; }
			else{ A[(numbVar + numbLim)*i + j] = 0.0; }
		}
	}



	//	for (i = 0; i < (numbVar + numbLim)*(numbLim); ++i)cout << A[i] << " ";////////////
	//	cout << endl;//////////////


	//■STEP2；初期値の生成
	
	double *xk;//探索する内点（ｋステップ目のｘの意味）
	xk = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	double buff =1.0;
	
	for (i = 0; i < numbVar + numbLim; ++i){
		if (numbEqu == 1 && i == 0){

			for (j = 0; j < numbVar; ++j){
				if (leftSide[numbVar *i + j] != 0.0){
					xk[i] = rightSide[j];
					buff = rightSide[j];
					break;
				}
			}
		}
		xk[i] = buff;
	}


	//スラック変数以外に１を代入して，内点の初期値を仮決定する
	for (i = 0; i < numbLim; ++i){
		for (j = 0; j <numbVar; ++j){
			rightSide[i] -= leftSide[numbVar*i + j] * xk[j];
		}
	}
		for (i = 0; i < 2; ++i)cout << rightSide[i] << " ";////////////
		cout << endl;/////

	//対角行列生成
	double *diag;//探索する内点（ｋステップ目のｘの意味）
	diag = (double *)malloc(sizeof(double)*numbLim*numbLim);
	for (i = 0; i < numbLim; ++i){
		for (j = 0; j <numbLim; ++j){
			if (i == j)diag[numbLim*i + j] = 1.0;
			else diag[numbLim*i + j] = 0.0;
		}
	}


	solvEquation(numbLim, diag, rightSide, xk);//スラック変数の初期値を求解
//		for (i = 0; i < (numbVar + numbLim); ++i)cout << xk[i] << " ";////////////

	//スラック変数の初期値を然るべき場所に移動する
	for (i = numbLim - 1; i >= 0; --i){
		xk[numbVar + i] = xk[i];
	}
	//初期値生成
	for (i = 0; i < numbVar; ++i){
		xk[i] = 1.0;
	}
	xk[i] = buff;
	for (i = 0; i < (numbVar + numbLim); ++i)cout << xk[i] << " ";////////////
	cout << endl;
	//等式制約の判定
	/*
	if (numbEqu != 0){
		for (i = 0; i < numbEqu; ++i){
			for (j = 0; j < numbVar; ++j){
				if (A[(numbVar + numbLim)*i + j] != 0.0){
					A[(numbVar + numbLim)*i + numbVar + j] = 0.0;//等式条件部分のスラック変数をゼロにする
					xk[j] = xk[numbVar + j ];
					xk[numbVar + j + 1] = 0.0;
					break;
				}
			}
		}
	}*/
	//検査++++
//	for (i = 0; i < (numbVar + numbLim); ++i)cout << xk[i] << " ";////////////
//	cout << endl;
//	for (i = 0; i < (numbVar + numbLim)*(numbLim); ++i)cout << A[i] << " ";////////////
//		cout << endl;//////////////

	//■STEP3；他の係数の準備
	double *A_tra;	//転置行列
	A_tra = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbLim));
	for (i = 0; i < numbVar + numbLim; i++){//転置行列の計算
		for (j = 0; j < numbLim; j++){
			A_tra[numbLim*i + j] = A[(numbVar + numbLim)*j + i];
		}
	}
	//	for (i = 0; i < (numbVar + numbLim)*numbLim; ++i)cout << A_tra[i] << " ";////////////


	double *C;
	C = (double *)malloc(sizeof(double)*(numbVar + numbVar));
	for (i = 0; i < numbVar + numbVar; i++){
		if (i<numbVar){
			C[i] = objective[i];
		}
		else{
			C[i] = 0.0;
		}

	}

	//結果や途中家計算を格納するための配列
	double *AX2At;//方程式の左辺　
	AX2At = (double *)malloc(sizeof(double)* numbLim* numbLim);
	double *AX2c;//方程式の右辺
	AX2c = (double *)malloc(sizeof(double)*(numbVar + numbLim));

	double *diagMatrix;//対角行列を格納する
	diagMatrix = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));
	double *diagMatrixSqua;//対角行列を2乗して格納する
	diagMatrixSqua = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));

	double *solu1;//途中計算
	solu1 = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbLim));
	double *solu2;//途中計算
	solu2 = (double *)malloc(sizeof(double)*numbLim);
	double *solu3;//途中計算
	solu3 = (double *)malloc(sizeof(double)*(numbVar + numbLim));

	double *Z;
	Z = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));

	double *dx;//更新則
	dx = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	double z_nol = 0;
	double z1 = 0;

	double upsilon = 1.0e-18;//終了条件　ε
	double alpha = 0.5;//更新則のゲイン

	double result = 0.0;

	//■STEP4；計算
	while (1){

		//左辺
		//対角行列生成
		for (i = 0; i < numbVar + numbLim; ++i){
			for (j = 0; j < numbVar + numbLim; ++j){
				diagMatrix[i * (numbVar + numbLim) + j] = 0.0;
				if (j == i)diagMatrix[i * (numbVar + numbLim) + j] = xk[i];
				//				cout << diagMatrix[i * (numbVar + numbLim) + j] << " ";/////////////
			}
		}

		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, numbVar + numbLim, diagMatrix, diagMatrix, diagMatrixSqua);//対角行列の平方
		calcMatrixMultiplication(numbLim, numbVar + numbLim, numbVar + numbLim, A, diagMatrixSqua, solu1);//
//				for (i = 0; i < 4; ++i)cout <<solu1[i] << " ";////////////
//		for (i = 0; i < (numbVar + numbLim)*numbLim; ++i)cout << A[i] << " ";////////////
		calcMatrixMultiplication(numbLim, numbVar + numbLim, numbLim, solu1, A_tra, AX2At);
			//	for (i = 0; i < 4; ++i)cout << AX2At[i] << " ";////////////
//				cout << endl;

		//右辺
		calcMatrixMultiplication(numbLim, numbVar + numbLim, 1, solu1, C, AX2c);
		solvEquation(numbLim, AX2At, AX2c, solu2);
//				for (i = 0; i < 4; ++i)cout << solu2[i] << " ";////////////ok
//				cout << endl;//////////////

		calcMatrixMultiplication(numbVar + numbLim, numbLim, 1, A_tra, solu2, solu3);

		//		for (i = 0; i < numbVar + numbLim; ++i)cout <<AX2c[i]<<" ";//////////////
		//		cout << endl;//////////////

		for (i = 0; i < numbVar + numbLim; ++i){
			Z[i] = C[i] - solu3[i];
		}
//				for (i = 0; i < numbVar + numbLim; ++i)cout << Z[i] << " ";////////////
//				cout << "  "<<endl;//////////////

		for (i = 0; i < numbVar + numbLim; ++i){
			z_nol += (diagMatrix[(numbVar + numbLim) * i + i] * Z[i])*(diagMatrix[(numbVar + numbLim) * i + i] * Z[i]);

		}
		z_nol = sqrt(z_nol);
		//		cout << z_nol << endl;

		//更新則を計算する
		for (i = 0; i < numbVar + numbLim; ++i){
			dx[i] = diagMatrixSqua[(numbVar + numbLim) * i + i] * Z[i] / z_nol;
		}
//				for (i = 0; i < numbVar + numbLim; ++i)cout << dx[i] << " ";////////////
//						cout << endl;//////////////

		//ｘの値を更新する
		for (i = 0; i < numbVar + numbLim; ++i){
			xk[i] += alpha * dx[i];
		}
				for (i = 0; i < numbVar + numbLim; ++i)cout << xk[i] << " ";////////////
				cout << endl;//////////////
		//		if (xk[0] <= 0.00353){////////
		//			cout << "Stop" << endl;////////
		//		}//////////////


		//判定式の計算
		for (i = 0; i < numbVar + numbLim; ++i){
			z1 += dx[i] * dx[i];
		}
		//		cout << z1 << endl;
		//		for (i = 75; i < 100; ++i)cout << dx[i] << " ";////////////
		//		cout << endl;//////////////
		if (upsilon>z1)break;
		z1 = 0.0;
	}

	//●STEP5；計算結果を返す（返す用の配列に渡す）
	
	for (i = 0; i < numbVar; ++i){
	optimumSolution[i] = xk[i];
	}
	
/*	for (i = 0; i < numbVar; ++i){
		result += optimumSolution[i] * objective[i];
	}
	*/
	return result;
}
