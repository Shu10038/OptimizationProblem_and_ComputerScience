/*
シンプレックス法を用いた線形計画法のソルバ
This is a solver of linear programming .It's used the simplex method.
’16/01/18　S.OHSAWA

*/

/******************************　　
線形計画法のソルバ
******************************/

#include"stdafx.h"
#include<stdlib.h>
#include<iostream>
#include<math.h>
using namespace std;

double LP_SOLVER01(
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
	tabl = (double *)malloc(sizeof(double)*(numbVar + numbLim + 1)*(numbLim + 1));

	//■STEP1；標準化
	//スラック変数を挿入してシンプレックス・タブローを作成する
	//目的関数の列の作成
	for (i = 0; i < numbVar + numbLim + 1; ++i){
		if (i<numbVar){
			tabl[i] = -objective[i];
		}
		else{ tabl[i] = 0.0; }
	}

	//制約条件の列の作成
	for (i = 0; i < numbLim; ++i){
		for (j = 0; j < numbVar + numbLim + 1; ++j){
			if (j < numbVar){ tabl[(numbVar + numbLim + 1)*(i + 1) + j] = leftSide[numbVar*i + j]; }
			else if (j - numbVar == i){ tabl[(numbVar + numbLim + 1)*(i + 1) + j] = 1.0; }
			else if (j == numbVar + numbLim){ tabl[(numbVar + numbLim + 1)*(i + 2) - 1] = rightSide[i]; }
			else{ tabl[(numbVar + numbLim + 1)*(i + 1) + j] = 0.0; }
		}
	}

	int minColmNumb = 0, minRowNumb;
	double minColu;
	double pivo, mult, minRowVari, terVari = 0.0;

	while (1){
		//最上段の中から最小の列を選ぶ
		minColu = tabl[0];
		minColmNumb = 0;
		for (i = 0; i <numbVar + numbLim; ++i){//端以外
			if (minColu > tabl[i]){
				minColu = tabl[i];
				minColmNumb = i;
			}
		}
		if (minColu >= 0)break;

		//定数項の中から最小の列を選ぶ
		minRowVari = tabl[(numbVar + numbLim + 1) * 2 - 1] / tabl[(numbVar + numbLim + 1) * 1 + minColmNumb];
		minRowNumb = 1;
		for (i = 0; i < numbLim; ++i){
			if (minRowVari > tabl[(numbVar + numbLim + 1)*(i + 2) - 1] / tabl[(numbVar + numbLim + 1)*(i + 1) + minColmNumb])
			{
				minRowVari = tabl[(numbVar + numbLim + 1)*(i + 2) - 1] / tabl[(numbVar + numbLim + 1)*(i + 1) + minColmNumb];
				minRowNumb = i + 1;
			}
		}

		// 正規化(最初の項の係数を１にする)
		pivo = tabl[(numbVar + numbLim + 1)*minRowNumb + minColmNumb];
		for (i = 0; i < numbVar + numbLim + 1; ++i){
			tabl[(numbVar + numbLim + 1)*minRowNumb + i] /= pivo;
		}

		//掃出し操作
		for (i = 0; i < 1 + numbLim; ++i){
			if (i != minRowNumb){
				mult = tabl[(numbVar + numbLim + 1)*i + minColmNumb];
				for (j = 0; j < numbVar + numbLim + 1; ++j){
					tabl[(numbVar + numbLim + 1)*i + j] -= mult * tabl[(numbVar + numbLim + 1)*minRowNumb + j];
				}
			}
		}
	}
	int buff = 0;
	for (i = 0; i < numbVar ; ++i){
		for (j = 0; j < numbLim + 1; ++j){
			if (tabl[(numbVar + numbLim + 1)*j+i] == 1){
				buff = j;
			}
		}if (buff != 0){ optimumSolution[i] = tabl[(numbVar + numbLim + 1)*(buff + 1) - 1]; }
		else{ optimumSolution[i] = 0.0; }
	}

	return tabl[(numbVar + numbLim + 1) - 1];
}



