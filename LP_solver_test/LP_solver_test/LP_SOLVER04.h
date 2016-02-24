/*
インフィージブル内点法を用いた線形計画法のソルバ
This is a solver of linear programming .It's used the infeasible method.
’16/01/29　S.OHSAWA

’16/02/13 修正
*/

#include"stdafx.h"
#include<stdlib.h>
#include<iostream>
#include<math.h>
using namespace std;


/******************************　
対角の逆行列
******************************/
int diagInv(
	int numb,		// 要素数
	double solution[])	// 入力される行列，兼，結果を格納する配列
{
	int i, j;

	for (i = 0; i < numb; i++) {
		for (j = 0; j < numb; j++) {
			if (i == j) solution[numb*i + j] = 1 / solution[numb*i + j];
		}
	}

	return 0;
}


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

/******************************　　
線形計画法のソルバ
******************************/
double LP_SOLVER04(
	int numbVar,			//変数の数 
	int numbLim,			//制約条件の数
	int numbEqu,			//等式条件
	double objective[],		//目的関数の係数　例)3x_1+2x_2　だったら{3.0,2.0}
	double leftSide[],		//制約条件の左側．例)2x_1+4x_2≦4 と x_1+3x_2≦5　だったら{2.0,4.0,1.0,3.0}
	double rightSide[],		//制約条件の右側．
	double optimumSolution[])
{

	int i, j;//カウンタ
	double result = 0.0;

	double *A;
	A = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbLim));

	double alpha_P = 0.5;
	double alpha_D = 0.5;
	double deci = 0.0;

	//■STEP1；標準化
	//スラック変数を挿入して標準化する
	for (i = 0; i < numbLim; ++i){
		for (j = 0; j < numbVar + numbLim; ++j){
			if (j < numbVar){ A[(numbVar + numbLim)*i + j] = leftSide[numbVar*i + j]; }
			else if (j - numbVar == i){ A[(numbVar + numbLim)*i + j] = 1.0; }
			else{ A[(numbVar + numbLim)*i + j] = 0.0; }
		}
	}


	double *A_tra;	//転置行列
	A_tra = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbLim));
	for (i = 0; i < numbVar + numbLim; i++){//転置行列の計算
		for (j = 0; j < numbLim; j++){
			A_tra[numbLim*i + j] = A[(numbVar + numbLim)*j + i];
		}
	}

	double *C;	//
	C = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	for (i = 0; i < numbVar + numbLim; i++){
		if (i < numbVar){ C[i] = -objective[i]; }
		else{ C[i] = 0.0; }
	}

	//対角行列生成
	double *diag;//探索する内点（ｋステップ目のｘの意味）
	diag = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));
	for (i = 0; i < numbVar + numbLim; ++i){
		for (j = 0; j <numbVar + numbLim; ++j){
			if (i == j)diag[(numbVar + numbLim)*i + j] = 1.0;
			else diag[(numbVar + numbLim)*i + j] = 0.0;
		}
	}

	//■STEP2；初期値の生成
	double *dx; dx = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	double *xk; xk = (double *)malloc(sizeof(double)*(numbVar + numbLim));//探索する内点（ｋステップ目のｘの意味）
	double buff = 1.0;
	for (i = 0; i < numbVar + numbLim; ++i){
		xk[i] = 1.0;
	}

	double *dz; dz = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	double *zk; zk = (double *)malloc(sizeof(double)*(numbVar + numbLim));//探索する内点（ｋステップ目のzの意味）
	for (i = 0; i < numbVar + numbLim; ++i){
		zk[i] = 1.0;
	}

	double *dy; dy = (double *)malloc(sizeof(double)*(numbLim));
	double *yk; yk = (double *)malloc(sizeof(double)*(numbLim));//探索する内点（ｋステップ目のzの意味）
	for (i = 0; i < numbLim; ++i){
		yk[i] = 0.0;
	}

	//	double *mu; mu = (double *)malloc(sizeof(double)*(numbVar + numbLim));



	double *X_diag; X_diag = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));//Xの対角行列
	double *Z_diag; Z_diag = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));//Xの対角行列




	//⊿yを求めるために使う配列
	double *dy_L_t1; dy_L_t1 = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));
	double *dy_L_t2; dy_L_t2 = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));

	double *Axb; Axb = (double *)malloc(sizeof(double)*(numbLim)*(numbLim));
	double *Ayzc; Ayzc = (double *)malloc(sizeof(double)*(2 * numbLim));
	double *X0; X0 = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));
	double *MAT2; MAT2 = (double *)malloc(sizeof(double)*numbLim*numbLim);


	double *LEFT; LEFT = (double *)malloc(sizeof(double)*(numbLim)*(numbLim));
	double *RIGHT; RIGHT = (double *)malloc(sizeof(double)*(numbLim));

	double *Xkz; Xkz = (double *)malloc(sizeof(double)*(numbLim));

	double *dy_R_1st_t1; dy_R_1st_t1 = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));
	double *dy_R_1st_t2; dy_R_1st_t2 = (double *)malloc(sizeof(double)*(numbVar + numbLim));


	double *dy_Aterm; dy_Aterm = (double *)malloc(sizeof(double)*(numbLim)*(numbVar + numbLim));
	double *dy_term; dy_term = (double *)malloc(sizeof(double)*(numbLim));
	double *dy_1st; dy_1st = (double *)malloc(sizeof(double)*(numbLim));

	double *dy_3rd; dy_3rd = (double *)malloc(sizeof(double)*(numbLim));


	//⊿zを求めるために使う変数
	double *dz_1st; dz_1st = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	double *dz_2nd; dz_2nd = (double *)malloc(sizeof(double)*(numbVar + numbLim));


	//⊿xを求めるために使う変数
	double *dx_1st_t1; dx_1st_t1 = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));
	double *dx_1st_t2; dx_1st_t2 = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	double *dx_2nd; dx_2nd = (double *)malloc(sizeof(double)*(numbVar + numbLim));

	//μの生成
	double mu = 0.0;
	for (i = 0; i < numbVar + numbLim; ++i){
		mu += (xk[i] * zk[i]);
	}

	while (1){
		//μの更新
		mu *= 0.9;

		//Xの対角行列生成
		for (i = 0; i < numbVar + numbLim; ++i){
			for (j = 0; j <numbVar + numbLim; ++j){
				if (i == j)X_diag[(numbVar + numbLim)*i + j] = xk[i];
				else X_diag[(numbVar + numbLim)*i + j] = 0.0;
			}
		}
		//Zの対角行列生成
		for (i = 0; i < numbVar + numbLim; ++i){
			for (j = 0; j <numbVar + numbLim; ++j){
				if (i == j)Z_diag[(numbVar + numbLim)*i + j] = zk[i];
				else Z_diag[(numbVar + numbLim)*i + j] = 0.0;
			}
		}
		//Zの逆行列
		diagInv(numbVar + numbLim, Z_diag);

		//パス追跡法による点の移動
		//Ax-ｂの計算
		calcMatrixMultiplication(numbLim, numbVar + numbLim, 1, A, xk, Axb);
		for (i = 0; i < numbVar + numbLim; ++i){
			Axb[i] -= rightSide[i];
		}

		//A^Ty+z-cの計算
		calcMatrixMultiplication(numbVar + numbLim, numbLim, 1, A_tra, yk, Ayzc);
		for (i = 0; i < numbVar + numbLim; ++i){
			Ayzc[i] = Ayzc[i] + zk[i] - C[i];
		}


		//⊿yを求める；方程式の左辺の計算
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, numbVar + numbLim, A, Z_diag, dy_L_t1);
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, numbVar + numbLim, dy_L_t1, X_diag, dy_L_t2);
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, numbLim, dy_L_t2, A_tra, LEFT);


		//⊿yを求める；方程式の右辺の計算
		//第1項の計算
		calcMatrixMultiplication(numbLim, numbVar + numbLim, numbVar + numbLim, A, Z_diag, dy_R_1st_t1);
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, 1, X_diag, zk, Xkz);
		for (i = 0; i < numbVar + numbLim; ++i){
			dy_R_1st_t2[i] = mu - Xkz[i];
		}


		calcMatrixMultiplication(numbLim, numbVar + numbLim, 1, dy_L_t1, dy_R_1st_t2, dy_1st);

		//第3項の計算 (第2項は Ax-ｂ なので既に計算した)
		calcMatrixMultiplication(numbLim, numbVar + numbLim, 1, dy_L_t2/*左辺の計算を流用する*/, Ayzc, dy_3rd);

		//右辺の計算
		for (i = 0; i < numbLim; ++i){
			RIGHT[i] = -dy_1st[i] - Axb[i] - dy_3rd[i];
		}

		//⊿yの求解
		solvEquation(numbLim, LEFT, RIGHT, dy);


		//⊿zを求める
		calcMatrixMultiplication(numbVar + numbLim, numbLim, 1, A_tra, dy, dz_1st);
		//⊿zの求解
		for (i = 0; i < numbVar + numbLim; ++i){
			dz[i] = -dz_1st[i] - Ayzc[i];
		}


		//⊿xを求める
		//第1項の計算
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, numbVar + numbLim, Z_diag, X_diag, dx_1st_t1);
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, 1, dx_1st_t1, dz, dx_1st_t2);

		//第2項の計算
		for (i = 0; i < numbVar + numbLim; ++i){
			for (j = 0; j < numbVar + numbLim; ++j){
				if (i == j)dx_2nd[i] = mu * Z_diag[(numbVar + numbLim)*i + j];
			}
		}

		//⊿xの求解
		for (i = 0; i < numbVar + numbLim; ++i){
			dx[i] = -dx_1st_t2[i] + dx_2nd[i] - xk[i];
		}


		//検査+++++
		/*
		for (i = 0; i < 4; ++i){
			cout << dx[i] << " ";
		}
		cout << endl;//+++++
		*/

		//x,y,zの値を更新する
		for (i = 0; i < numbVar + numbLim; ++i){
			xk[i] += alpha_P * dx[i];
			zk[i] += alpha_D * dz[i];
		}
		for (i = 0; i < numbLim; ++i){
			yk[i] += alpha_D * dy[i];
		}


		for (i = 0; i < 4; ++i){
			deci += dx[i] * dx[i];
		}
		if (deci < 1.0e-15)break;
		deci = 0.0;

	}

	//■STEP5；計算結果を返す（返す用の配列に渡す）

	for (i = 0; i < numbVar; ++i){
		optimumSolution[i] = xk[i];
		result += xk[i] + optimumSolution[i];
	}


	return result;
}
