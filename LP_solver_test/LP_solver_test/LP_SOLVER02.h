/*
２段シンプレックス法を用いた線形計画法のソルバ
This is a solver of linear programming .It's used the simplex method.
’16/02/19　S.OHSAWA

*/

#include"stdafx.h"
#include<stdlib.h>
#include<iostream>
#include<math.h>
using namespace std;

double first_culcSimplexMethod(
	int numbVar,			//変数の数 
	int numbLim,			//制約条件の数
	int numbEqu,			//等式条件
	double tableau[],		
	double optimumSolution[])
{
	int i, j, k;//カウンタ

	int minColmNumb = 0, minRowNumb;
	double minColu;
	double pivo, mult, minRowVari, terVari = 0.0;

	while (1){
		//最上段の中から最小の列を選ぶ
		minColu = tableau[0];
		minColmNumb = 0;
		for (i = 0; i <numbVar + 2 * numbLim; ++i){//端以外
			if (minColu >= tableau[i]){
				minColu = tableau[i];
				minColmNumb = i;
			}
		}
		if (minColu >= 0)break;

		//検査+++++
		
		for (i = 0; i < 1 + numbLim; ++i){
			for (j = 0; j < numbVar + 2 * numbLim + 1; ++j){
				cout << tableau[(numbVar + 2 * numbLim + 1)*i + j] << " ";
			}
			cout << endl;
		}
		cout << endl;
		

		//定数項の中から最小の列を選ぶ
		minRowVari = tableau[(numbVar + 2 * numbLim + 1) * 2 - 1] / tableau[(numbVar + 2 * numbLim + 1) * 1 + minColmNumb];
		minRowNumb = 1;
		for (i = 0; i < numbLim; ++i){
			if (minRowVari > tableau[(numbVar + 2 * numbLim + 1)*(i + 2) - 1] / tableau[(numbVar + 2 * numbLim + 1)*(i + 1) + minColmNumb])
			{
				minRowVari = tableau[(numbVar + 2 * numbLim + 1)*(i + 2) - 1] / tableau[(numbVar + 2 * numbLim + 1)*(i + 1) + minColmNumb];
				minRowNumb = i + 1;
			}
		}

		// 正規化(最初の項の係数を１にする)
		pivo = tableau[(numbVar + 2 * numbLim + 1)*minRowNumb + minColmNumb];
		for (i = 0; i < numbVar + 2 * numbLim + 1; ++i){
			tableau[(numbVar + 2 * numbLim + 1)*minRowNumb + i] /= pivo;
		}


		//
		for (i = 0; i < 1 + numbLim; ++i){
			if (i != minRowNumb){
				mult = tableau[(numbVar + 2 * numbLim + 1)*i + minColmNumb];
				for (j = 0; j < numbVar + 2 * numbLim + 1; ++j){
					tableau[(numbVar + 2 * numbLim + 1)*i + j] -= mult * tableau[(numbVar + 2 * numbLim + 1)*minRowNumb + j];
				}
			}
		}
		//検査+++++
		for (i = 0; i < 1 + numbLim; ++i){
			for (j = 0; j < numbVar + 2 * numbLim + 1; ++j){
				cout << tableau[(numbVar + 2 * numbLim + 1)*i + j] << " ";
			}
			cout << endl;
		}
		cout << endl;

	}
	for (i = 0; i < numbVar+1; ++i){
		for (j = 0; j < numbLim + 1; ++j){
			if (tableau[(numbVar + 2 * numbLim + 1)*j + i] == 1){
				optimumSolution[i] = tableau[(numbVar + 2 * numbLim + 1)*(j + 1) - 1];
			}
		}
	}

	return 0;
}




/******************************　　
線形計画法のソルバ
******************************/
double LP_SOLVER02(
	int numbVar,			//変数の数 
	int numbLim,			//制約条件の数
	int numbEqu,			//等式条件
	double objective[],		//目的関数の係数　例)3x_1+2x_2　だったら{3.0,2.0}
	double leftSide[],		//制約条件の左側．例)2x_1+4x_2≦4 と x_1+3x_2≦5　だったら{2.0,4.0,1.0,3.0}
	double rightSide[],		//制約条件の右側．
	double optimumSolution[])
{

	int i, j, k;//カウンタ

	double *tabl;
	tabl = (double *)malloc(sizeof(double)*(numbVar + 2 * numbLim + 1)*(numbLim + 1));

	//■STEP1；標準化
	//サープラス変数を挿入してシンプレックス・タブローを作成する
	//目的関数の列の作成
	for (i = 0; i < numbVar + 2 * numbLim; ++i){
		if (i>numbVar + numbLim - 1){
			tabl[i] = 1.0;
		}
		else{ tabl[i] = 0.0; }
	}
	//制約条件の列の作成
	for (i = 0; i < numbLim; ++i){
		for (j = 0; j < numbVar + 2 * numbLim + 1; ++j){
			if (j < numbVar){ tabl[(numbVar + 2 * numbLim + 1)*(i + 1) + j] = leftSide[numbVar*i + j]; }
			else if (j - numbVar == i){ tabl[(numbVar + 2 * numbLim + 1)*(i + 1) + j] = -1.0; }
			else if (j - numbVar - numbLim == i){ tabl[(numbVar + 2 * numbLim + 1)*(i + 1) + j] = 1.0; }
			else if (j - numbVar - 2 * numbLim + 1 == 1){ tabl[(numbVar + 2 * numbLim + 1)*(i + 2) - 1] = rightSide[i]; }
			else{ tabl[(numbVar + 2 * numbLim + 1)*(i + 1) + j] = 0.0; }
		}
	}
	tabl[numbVar + 2 * numbLim] = 0.0;


	//検査+++++

	for (i = 0; i < 1 + numbLim; ++i){
		for (j = 0; j < numbVar + 2 * numbLim + 1; ++j){
			cout << tabl[(numbVar + 2 * numbLim + 1)*i + j] << " ";
		}
		cout << endl;
	}

	//一段目を制約条件で引く
	for (i = 0; i < numbLim; ++i){
		for (j = 0; j < numbVar + 2 * numbLim + 1; ++j){
			tabl[j] -= tabl[(numbVar + 2 * numbLim + 1)*(i + 1) + j];
		}
	}

	first_culcSimplexMethod(numbVar,  numbLim, numbEqu, tabl, optimumSolution);

		//検査+++++

	for (i = 0; i < 1 + numbLim; ++i){
		cout << optimumSolution[i] << endl;
	}

	int minColmNumb = 0, minRowNumb;
	double minColu;
	double pivo, mult, terNumb, terVari = 0.0;


	while (1){
		//
		minColu = tabl[0];
		for (i = 0; i <numbVar + 2 * numbLim; ++i){//端以外
			if (minColu > tabl[i]){
				minColu = tabl[i];
				minColmNumb = i;
			}
		}
		if (minColu > 0)break;
		//
		for (i = 0; i < numbLim; ++i){
			tabl[(numbVar + 2 * numbLim + 1)*(i + 2) - 1] /= tabl[(numbVar + 2 * numbLim + 1)*(i + 1) + minColmNumb];
		}

		//
		terVari = tabl[(numbVar + 2 * numbLim + 1)*(i + 2) - 1];
		minRowNumb = 1;
		for (i = 0; i < numbLim; ++i){
			if (terVari > tabl[(numbVar + 2 * numbLim + 1)*(i + 2) - 1])
			{
				terNumb = tabl[(numbVar + 2 * numbLim + 1)*(i + 2) - 1];
				minRowNumb = i + 1;
			}
		}
		// 正規化(最初の項の係数を１にする)
		pivo = tabl[(numbVar + 2 * numbLim + 1)*minRowNumb + minColmNumb];
		for (i = 0; i < numbVar + 2 * numbLim; ++i){
			tabl[(numbVar + 2 * numbLim + 1)*minRowNumb + i] /= pivo;
		}


		//
		for (i = 0; i < 1 + numbLim; ++i){
			if (i != minRowNumb){
				mult = tabl[(numbVar + 2 * numbLim + 1)*i + minColmNumb];
				for (j = 0; j < numbVar + 2 * numbLim + 1; ++j){
					tabl[(numbVar + 2 * numbLim + 1)*i + j] -= mult * tabl[(numbVar + 2 * numbLim + 1)*minRowNumb + j];
				}
			}
		}


		for (i = 0; i < 1 + numbLim; ++i){
			for (j = 0; j < numbVar + 2 * numbLim + 1; ++j){
				cout << tabl[(numbVar + 2 * numbLim + 1)*i + j] << " ";
			}
			cout << endl;
		}
		cout << endl;

	}


	cout << minColmNumb;
	return 0;
}



