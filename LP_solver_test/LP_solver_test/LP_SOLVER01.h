/*
�V���v���b�N�X�@��p�������`�v��@�̃\���o
This is a solver of linear programming .It's used the simplex method.
�f16/01/18�@S.OHSAWA

*/

/******************************�@�@
���`�v��@�̃\���o
******************************/

#include"stdafx.h"
#include<stdlib.h>
#include<iostream>
#include<math.h>
using namespace std;

double LP_SOLVER01(
	int numbVar,			//�ϐ��̐� 
	int numbLim,			//��������̐�
	int numbEqu,			//��������
	double objective[],		//�ړI�֐��̌W���@��)3x_1+2x_2�@��������{3.0,2.0}
	double leftSide[],		//��������̍����D��)2x_1+4x_2��4 �� x_1+3x_2��5�@��������{2.0,4.0,1.0,3.0}
	double rightSide[],		//��������̉E���D
	double optimumSolution[])
{
	int i, j, k;//�J�E���^

	double *tabl;
	tabl = (double *)malloc(sizeof(double)*(numbVar + numbLim + 1)*(numbLim + 1));

	//��STEP1�G�W����
	//�X���b�N�ϐ���}�����ăV���v���b�N�X�E�^�u���[���쐬����
	//�ړI�֐��̗�̍쐬
	for (i = 0; i < numbVar + numbLim + 1; ++i){
		if (i<numbVar){
			tabl[i] = -objective[i];
		}
		else{ tabl[i] = 0.0; }
	}

	//��������̗�̍쐬
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
		//�ŏ�i�̒�����ŏ��̗��I��
		minColu = tabl[0];
		minColmNumb = 0;
		for (i = 0; i <numbVar + numbLim; ++i){//�[�ȊO
			if (minColu > tabl[i]){
				minColu = tabl[i];
				minColmNumb = i;
			}
		}
		if (minColu >= 0)break;

		//�萔���̒�����ŏ��̗��I��
		minRowVari = tabl[(numbVar + numbLim + 1) * 2 - 1] / tabl[(numbVar + numbLim + 1) * 1 + minColmNumb];
		minRowNumb = 1;
		for (i = 0; i < numbLim; ++i){
			if (minRowVari > tabl[(numbVar + numbLim + 1)*(i + 2) - 1] / tabl[(numbVar + numbLim + 1)*(i + 1) + minColmNumb])
			{
				minRowVari = tabl[(numbVar + numbLim + 1)*(i + 2) - 1] / tabl[(numbVar + numbLim + 1)*(i + 1) + minColmNumb];
				minRowNumb = i + 1;
			}
		}

		// ���K��(�ŏ��̍��̌W�����P�ɂ���)
		pivo = tabl[(numbVar + numbLim + 1)*minRowNumb + minColmNumb];
		for (i = 0; i < numbVar + numbLim + 1; ++i){
			tabl[(numbVar + numbLim + 1)*minRowNumb + i] /= pivo;
		}

		//�|�o������
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



