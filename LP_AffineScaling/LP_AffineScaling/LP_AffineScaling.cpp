#include"stdafx.h"
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include<math.h>
using namespace std;


/**********�@�s��̊|���Z�@**********/
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

/**********�@�������̌v�Z�@**********/
double solvEquation(
	int size,			// �ϐ��̐�
	double *A,			// �����s��i�A���������̍s��A�̕����j
	double *B,			// �c�x�N�g���i�A����������B�̕����j
	double *solution)	// ���ʂ�Ԃ��z��
{
	int i, j, k, n;
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

/**********�@���C���֐��@**********/
int main(){
	ofstream fout("output.txt");
	if (!fout){
		cout << "This file can not open!!";
		return 1;
	}

	int i, j;//�J�E���^

	//�����l�̐ݒ�

	double A[8] = { 2.0, 1.0, 1.0, 0.0,
		1.0, 3.0, 0.0, 1.0 };

	double A_tra[8];	//�]�u�s��
	for (i = 0; i < 4; i++){//�]�u�s��̌v�Z
		for (j = 0; j < 2; j++){
			A_tra[j + 2 * i] = A[i + j * 4];
		}
	}
	double C[4] = { 4.0, 5.0, 0.0, 0.0 };//��������̎��̉E��
	double xk[4] = { 1.0, 1.0, 1.0, 1.0 };//�����l�i�������_�͌v�Z���ďo�����j

	//���ʂ�r���ƌv�Z���i�[���邽�߂̔z��
	double AX2At[4];//�������̍��Ӂ@
	double AX2c[2];//�������̉E��

	double diagMatrix[16];//�Ίp�s����i�[����
	double diagMatrixSqua[16];//�Ίp�s���2�悵�Ċi�[����

	double solu1[8];//�r���v�Z
	double solu2[2];//�r���v�Z
	double solu3[4];//�r���v�Z
	double Z[4];

	double dx[4];//�X�V��
	double z_nol = 0;
	double z1 = 0;

	double upsilon = 1.0e-20;//�I�������@��
	double alpha = 0.8;//�X�V���Ɋ|����Q�C��

	while (1){
		//����
		//�Ίp�s�񐶐�
		for (i = 0; i < 4; ++i){
			for (j = 0; j < 4; ++j){
				diagMatrix[i * 4 + j] = 0.0;
				if (j == i)diagMatrix[i * 4 + j] = xk[i];
			}
		}

		calcMatrixMultiplication(4, 4, 4, diagMatrix, diagMatrix, diagMatrixSqua);//�Ίp�s��̕���
		calcMatrixMultiplication(2, 4, 4, A, diagMatrixSqua, solu1);//�Ίp�s��̕���
		calcMatrixMultiplication(2, 4, 2, solu1, A_tra, AX2At);

		//�E��
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

		//k+1�̓_�̍X�V
		for (i = 0; i < 4; ++i){
			xk[i] += +alpha * dx[i];
		}

		//��������̌v�Z
		for (i = 0; i < 4; ++i){
			z1 += z1 + dx[i] * dx[i];
		}
		if (upsilon>z1)break;
		z1 = 0.0;
	}
	cout << "x1 = " << xk[0] << "  X2 = " << xk[1] << endl;
	cout << "�ő�l�� " << xk[0] + xk[1] << endl;

	return 0;
}
