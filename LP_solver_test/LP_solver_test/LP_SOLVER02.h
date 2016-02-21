/*
�Q�i�V���v���b�N�X�@��p�������`�v��@�̃\���o
This is a solver of linear programming .It's used the simplex method.
�f16/02/19�@S.OHSAWA

*/

#include"stdafx.h"
#include<stdlib.h>
#include<iostream>
#include<math.h>
using namespace std;

double first_culcSimplexMethod(
	int numbVar,			//�ϐ��̐� 
	int numbLim,			//��������̐�
	int numbEqu,			//��������
	double tableau[],		
	double optimumSolution[])
{
	int i, j, k;//�J�E���^

	int minColmNumb = 0, minRowNumb;
	double minColu;
	double pivo, mult, minRowVari, terVari = 0.0;

	while (1){
		//�ŏ�i�̒�����ŏ��̗��I��
		minColu = tableau[0];
		minColmNumb = 0;
		for (i = 0; i <numbVar + 2 * numbLim; ++i){//�[�ȊO
			if (minColu >= tableau[i]){
				minColu = tableau[i];
				minColmNumb = i;
			}
		}
		if (minColu >= 0)break;

		//����+++++
		
		for (i = 0; i < 1 + numbLim; ++i){
			for (j = 0; j < numbVar + 2 * numbLim + 1; ++j){
				cout << tableau[(numbVar + 2 * numbLim + 1)*i + j] << " ";
			}
			cout << endl;
		}
		cout << endl;
		

		//�萔���̒�����ŏ��̗��I��
		minRowVari = tableau[(numbVar + 2 * numbLim + 1) * 2 - 1] / tableau[(numbVar + 2 * numbLim + 1) * 1 + minColmNumb];
		minRowNumb = 1;
		for (i = 0; i < numbLim; ++i){
			if (minRowVari > tableau[(numbVar + 2 * numbLim + 1)*(i + 2) - 1] / tableau[(numbVar + 2 * numbLim + 1)*(i + 1) + minColmNumb])
			{
				minRowVari = tableau[(numbVar + 2 * numbLim + 1)*(i + 2) - 1] / tableau[(numbVar + 2 * numbLim + 1)*(i + 1) + minColmNumb];
				minRowNumb = i + 1;
			}
		}

		// ���K��(�ŏ��̍��̌W�����P�ɂ���)
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
		//����+++++
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




/******************************�@�@
���`�v��@�̃\���o
******************************/
double LP_SOLVER02(
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
	tabl = (double *)malloc(sizeof(double)*(numbVar + 2 * numbLim + 1)*(numbLim + 1));

	//��STEP1�G�W����
	//�T�[�v���X�ϐ���}�����ăV���v���b�N�X�E�^�u���[���쐬����
	//�ړI�֐��̗�̍쐬
	for (i = 0; i < numbVar + 2 * numbLim; ++i){
		if (i>numbVar + numbLim - 1){
			tabl[i] = 1.0;
		}
		else{ tabl[i] = 0.0; }
	}
	//��������̗�̍쐬
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


	//����+++++

	for (i = 0; i < 1 + numbLim; ++i){
		for (j = 0; j < numbVar + 2 * numbLim + 1; ++j){
			cout << tabl[(numbVar + 2 * numbLim + 1)*i + j] << " ";
		}
		cout << endl;
	}

	//��i�ڂ𐧖�����ň���
	for (i = 0; i < numbLim; ++i){
		for (j = 0; j < numbVar + 2 * numbLim + 1; ++j){
			tabl[j] -= tabl[(numbVar + 2 * numbLim + 1)*(i + 1) + j];
		}
	}

	first_culcSimplexMethod(numbVar,  numbLim, numbEqu, tabl, optimumSolution);

		//����+++++

	for (i = 0; i < 1 + numbLim; ++i){
		cout << optimumSolution[i] << endl;
	}

	int minColmNumb = 0, minRowNumb;
	double minColu;
	double pivo, mult, terNumb, terVari = 0.0;


	while (1){
		//
		minColu = tabl[0];
		for (i = 0; i <numbVar + 2 * numbLim; ++i){//�[�ȊO
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
		// ���K��(�ŏ��̍��̌W�����P�ɂ���)
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



