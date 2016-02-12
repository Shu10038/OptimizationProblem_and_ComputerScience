/*
�A�t�B���ϊ��@��p�������`�v��@�̃\���o
This is a solver of linear programming .It's used the affine scaling method.
														�f16/01/29�@S.OHSAWA
*/

#include"stdafx.h"
#include<stdlib.h>
#include<iostream>
#include<math.h>
using namespace std;


/******************************�@
	�s��̊|���Z�@
******************************/
int calcMatrixMultiplication(
	int leftColumn,		// �c�i��Acolumn�j�G������|����s��̏c�̒���
	int leftRow,		// ���i�s�Arow�j   �G������|����s��̉��̒���
	int rightRow,		// ���i�s�Arow�j   �G�E����|����s��̉��̒���
	double *left,		// ������|����s��
	double *right,		// �E����|����s��
	double solution[])	// ���ʂ��i�[����z��
{
	int i, j, k;
	double term;

	//�|���Z
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
	�������̌v�Z�@
 ******************************/
double solvEquation(
	int size,			// �ϐ��̐�
	double *A,			// �����s��i�A���������̍s��A�̕����j
	double *B,			// �c�x�N�g���i�A����������B�̕����j
	double *solution)	// ���ʂ�Ԃ��z��
{
	int i, j, k, n, m;
	double *augmentedMatrix;
	double pivo, mult;

	augmentedMatrix = (double *)malloc(sizeof(double)*(size + 1)*size);
	if (augmentedMatrix == NULL) exit(0);

	//�g��W���s�񐶐�
	for (i = 0; i < size; ++i){
		for (j = 0; j < size; ++j){
			augmentedMatrix[i*(size + 1) + j] = A[i*size + j];
		}
	}
	for (i = 0; i < size; i++){
		augmentedMatrix[(size + 1)*(i + 1) - 1] = B[i];
	}

	//�K�E�X�̏����@�̌v�Z
	for (i = 0; i < size; i++)
	{
		// ���K��(�ŏ��̍��̌W�����P�ɂ���)
		pivo = augmentedMatrix[(size + 1)*i + i];
		for (j = 0; j < size + 1; ++j)
		{
			augmentedMatrix[(size + 1)*i + j] = (1 / pivo) * augmentedMatrix[(size + 1)*i + j];
		}

		// �t�O�p�`�̍s��𐶐�����
		// ���݂̍s��艺�̍s��i��ڂ�0�ɂȂ�悤�Ɉ����Z������
		for (k = i + 1; k < size; ++k)
		{
			mult = augmentedMatrix[(size + 1)*k + i];
			for (n = i; n < size + 1; ++n)
			{
				augmentedMatrix[(size + 1)*k + n] -= mult * augmentedMatrix[(size + 1)*i + n];
			}
		}
	}

	// ��i���(�e�ϐ����Ɨ������`�ɂȂ�悤�ɕό`����)
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

	//�v�Z���ʂ̎��o��
	for (i = 0; i < size; i++){
		solution[i] = augmentedMatrix[(size + 1) *(i + 1) - 1];
	}

	return 0;
}

/******************************�@�@
	���`�v��@�̃\���o
******************************/
int LP_SOLVER02(
	int numbVar,			//�ϐ��̐� 
	int numbLim,			//��������̐�
	double objective[],		//�ړI�֐��̌W���@��)3x_1+2x_2�@��������{3.0,2.0}
	double leftSide[],		//��������̍����D��)2x_1+4x_2��4 �� x_1+3x_2��5�@��������{2.0,4.0,1.0,3.0}
	double rightSide[],		//��������̉E���D
	double optimumSolution[])
{

	int i, j, k;//�J�E���^

	double *A;
	A = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbLim));

	//��STEP1�G�W����
	//�X���b�N�ϐ���}�����ĕW��������
	for (i = 0; i < numbVar; ++i){
		for (j = 0; j < numbVar + numbLim; ++j){
			if (j < numbVar){ A[(2 * numbVar)*i + j] = leftSide[numbVar*i + j]; }
			else if (j - numbVar == i){ A[(numbVar + numbLim)*i + j] = 1.0; }
			else{ A[(numbVar + numbLim)*i + j] = 0.0; }
		}
	}

	//��STEP2�G�����l�̐���
	//�X���b�N�ϐ��ȊO�ɂP�������āC���_�̏����l�������肷��
	for (i = 0; i < numbLim; ++i){
		for (j = 0; j <numbVar; ++j){
			rightSide[i] -= leftSide[numbVar*i + j];
		}
	}

	//�Ίp�s�񐶐�
	double *diag;//�T��������_�i���X�e�b�v�ڂ̂��̈Ӗ��j
	diag = (double *)malloc(sizeof(double)*numbLim*numbLim);
	for (i = 0; i < numbLim; ++i){
		for (j = 0; j <numbLim; ++j){
			if (i == j)diag[numbLim*i + j] = 1.0;
			else diag[numbLim*i + j] = 0.0;
		}
	}

	double *xk;//�T��������_�i���X�e�b�v�ڂ̂��̈Ӗ��j
	xk = (double *)malloc(sizeof(double)*(numbVar + numbLim));

	solvEquation(numbLim, diag, rightSide, xk);

	//�����l�ړ�
	for (i = numbLim-1; i >=0; --i){
		xk[numbVar + i] = xk[i];
	}
	//�����l����
	for (i = 0; i < numbVar; ++i){
		xk[i] = 1.0;
	}

	//��STEP3�G���̌W���̏���
	double *A_tra;	//�]�u�s��
	A_tra = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbLim));
	for (i = 0; i < numbVar + numbLim; i++){//�]�u�s��̌v�Z
		for (j = 0; j < numbLim; j++){
			A_tra[numbLim*i + j] = A[(numbVar + numbLim)*j+i];
		}
	}

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

	//���ʂ�r���ƌv�Z���i�[���邽�߂̔z��
	double *AX2At;//�������̍��Ӂ@
		AX2At = (double *)malloc(sizeof(double)* numbLim* numbLim);
	double *AX2c;//�������̉E��
		AX2c = (double *)malloc(sizeof(double)*(numbVar + numbLim));

	double *diagMatrix;//�Ίp�s����i�[����
		diagMatrix = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));
	double *diagMatrixSqua;//�Ίp�s���2�悵�Ċi�[����
		diagMatrixSqua = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));

	double *solu1;//�r���v�Z
		solu1 = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbLim));
	double *solu2;//�r���v�Z
		solu2 = (double *)malloc(sizeof(double)*numbLim);
	double *solu3;//�r���v�Z
		solu3 = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	
		double *Z;
		Z = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));

	double *dx;//�X�V��
		dx = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	double z_nol = 0;
	double z1 = 0;

	double upsilon = 1.0e-16;//�I�������@��
	double alpha = 0.8;//�X�V���̃Q�C��

	//��STEP4�G���_�@�ɂ����`�v��@
	while (1){
		
		//����
			//�Ίp�s�񐶐�
		for (i = 0; i < numbVar + numbLim; ++i){
			for (j = 0; j < numbVar + numbLim; ++j){
				diagMatrix[i * (numbVar + numbLim) + j] = 0.0;
				if (j == i)diagMatrix[i * (numbVar + numbLim) + j] = xk[i];
			}
		}
		
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, numbVar + numbLim, diagMatrix, diagMatrix, diagMatrixSqua);//�Ίp�s��̕���
		calcMatrixMultiplication(numbLim, numbVar + numbLim, numbVar + numbLim, A, diagMatrixSqua, solu1);//�Ίp�s��̕���
		calcMatrixMultiplication(numbLim, numbVar + numbLim, numbLim, solu1, A_tra, AX2At);

		//�E��
		calcMatrixMultiplication(numbLim, numbVar + numbLim, 1, solu1, C, AX2c);
		solvEquation(numbLim, AX2At, AX2c, solu2);

		calcMatrixMultiplication(numbVar + numbLim, numbLim, 1, A_tra, solu2, solu3);

		for (i = 0; i < numbVar + numbLim; ++i){
			Z[i] = C[i] - solu3[i];
		}

		for (i = 0; i < numbVar + numbLim; ++i){
			z_nol += (diagMatrix[(numbVar + numbLim) * i + i] * Z[i])*(diagMatrix[(numbVar + numbLim) * i + i] * Z[i]);
			
		}
		z_nol = sqrt(z_nol);

		//�X�V�����v�Z����
		for (i = 0; i < numbVar + numbLim; ++i){
			dx[i] = diagMatrixSqua[(numbVar + numbLim) * i + i] * Z[i] / z_nol;
		}

		//���̒l���X�V����
		for (i = 0; i < numbVar + numbLim; ++i){
			xk[i] += alpha * dx[i];
		}

		//���莮�̌v�Z
		for (i = 0; i < numbVar + numbLim; ++i){
			z1 +=  dx[i] * dx[i];
		}
		if (upsilon > z1)break;//�Èȉ��̏ꍇ�̓��[�v�E�o
		z1 = 0.0;
	}

	//��STEP5�G�v�Z���ʂ�Ԃ��i�Ԃ��p�̔z��ɓn���j
	for (i = 0; i < numbVar; ++i){
		optimumSolution[i] = xk[i];
	}

	return 0;
}
