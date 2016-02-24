/*
�C���t�B�[�W�u�����_�@��p�������`�v��@�̃\���o
This is a solver of linear programming .It's used the infeasible method.
�f16/01/29�@S.OHSAWA

�f16/02/13 �C��
*/

#include"stdafx.h"
#include<stdlib.h>
#include<iostream>
#include<math.h>
using namespace std;


/******************************�@
�Ίp�̋t�s��
******************************/
int diagInv(
	int numb,		// �v�f��
	double solution[])	// ���͂����s��C���C���ʂ��i�[����z��
{
	int i, j;

	for (i = 0; i < numb; i++) {
		for (j = 0; j < numb; j++) {
			if (i == j) solution[numb*i + j] = 1 / solution[numb*i + j];
		}
	}

	return 0;
}


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

/******************************�@�@
���`�v��@�̃\���o
******************************/
double LP_SOLVER04(
	int numbVar,			//�ϐ��̐� 
	int numbLim,			//��������̐�
	int numbEqu,			//��������
	double objective[],		//�ړI�֐��̌W���@��)3x_1+2x_2�@��������{3.0,2.0}
	double leftSide[],		//��������̍����D��)2x_1+4x_2��4 �� x_1+3x_2��5�@��������{2.0,4.0,1.0,3.0}
	double rightSide[],		//��������̉E���D
	double optimumSolution[])
{

	int i, j;//�J�E���^
	double result = 0.0;

	double *A;
	A = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbLim));

	double alpha_P = 0.5;
	double alpha_D = 0.5;
	double deci = 0.0;

	//��STEP1�G�W����
	//�X���b�N�ϐ���}�����ĕW��������
	for (i = 0; i < numbLim; ++i){
		for (j = 0; j < numbVar + numbLim; ++j){
			if (j < numbVar){ A[(numbVar + numbLim)*i + j] = leftSide[numbVar*i + j]; }
			else if (j - numbVar == i){ A[(numbVar + numbLim)*i + j] = 1.0; }
			else{ A[(numbVar + numbLim)*i + j] = 0.0; }
		}
	}


	double *A_tra;	//�]�u�s��
	A_tra = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbLim));
	for (i = 0; i < numbVar + numbLim; i++){//�]�u�s��̌v�Z
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

	//�Ίp�s�񐶐�
	double *diag;//�T��������_�i���X�e�b�v�ڂ̂��̈Ӗ��j
	diag = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));
	for (i = 0; i < numbVar + numbLim; ++i){
		for (j = 0; j <numbVar + numbLim; ++j){
			if (i == j)diag[(numbVar + numbLim)*i + j] = 1.0;
			else diag[(numbVar + numbLim)*i + j] = 0.0;
		}
	}

	//��STEP2�G�����l�̐���
	double *dx; dx = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	double *xk; xk = (double *)malloc(sizeof(double)*(numbVar + numbLim));//�T��������_�i���X�e�b�v�ڂ̂��̈Ӗ��j
	double buff = 1.0;
	for (i = 0; i < numbVar + numbLim; ++i){
		xk[i] = 1.0;
	}

	double *dz; dz = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	double *zk; zk = (double *)malloc(sizeof(double)*(numbVar + numbLim));//�T��������_�i���X�e�b�v�ڂ�z�̈Ӗ��j
	for (i = 0; i < numbVar + numbLim; ++i){
		zk[i] = 1.0;
	}

	double *dy; dy = (double *)malloc(sizeof(double)*(numbLim));
	double *yk; yk = (double *)malloc(sizeof(double)*(numbLim));//�T��������_�i���X�e�b�v�ڂ�z�̈Ӗ��j
	for (i = 0; i < numbLim; ++i){
		yk[i] = 0.0;
	}

	//	double *mu; mu = (double *)malloc(sizeof(double)*(numbVar + numbLim));



	double *X_diag; X_diag = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));//X�̑Ίp�s��
	double *Z_diag; Z_diag = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));//X�̑Ίp�s��




	//��y�����߂邽�߂Ɏg���z��
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


	//��z�����߂邽�߂Ɏg���ϐ�
	double *dz_1st; dz_1st = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	double *dz_2nd; dz_2nd = (double *)malloc(sizeof(double)*(numbVar + numbLim));


	//��x�����߂邽�߂Ɏg���ϐ�
	double *dx_1st_t1; dx_1st_t1 = (double *)malloc(sizeof(double)*(numbVar + numbLim)*(numbVar + numbLim));
	double *dx_1st_t2; dx_1st_t2 = (double *)malloc(sizeof(double)*(numbVar + numbLim));
	double *dx_2nd; dx_2nd = (double *)malloc(sizeof(double)*(numbVar + numbLim));

	//�ʂ̐���
	double mu = 0.0;
	for (i = 0; i < numbVar + numbLim; ++i){
		mu += (xk[i] * zk[i]);
	}

	while (1){
		//�ʂ̍X�V
		mu *= 0.9;

		//X�̑Ίp�s�񐶐�
		for (i = 0; i < numbVar + numbLim; ++i){
			for (j = 0; j <numbVar + numbLim; ++j){
				if (i == j)X_diag[(numbVar + numbLim)*i + j] = xk[i];
				else X_diag[(numbVar + numbLim)*i + j] = 0.0;
			}
		}
		//Z�̑Ίp�s�񐶐�
		for (i = 0; i < numbVar + numbLim; ++i){
			for (j = 0; j <numbVar + numbLim; ++j){
				if (i == j)Z_diag[(numbVar + numbLim)*i + j] = zk[i];
				else Z_diag[(numbVar + numbLim)*i + j] = 0.0;
			}
		}
		//Z�̋t�s��
		diagInv(numbVar + numbLim, Z_diag);

		//�p�X�ǐՖ@�ɂ��_�̈ړ�
		//Ax-���̌v�Z
		calcMatrixMultiplication(numbLim, numbVar + numbLim, 1, A, xk, Axb);
		for (i = 0; i < numbVar + numbLim; ++i){
			Axb[i] -= rightSide[i];
		}

		//A^Ty+z-c�̌v�Z
		calcMatrixMultiplication(numbVar + numbLim, numbLim, 1, A_tra, yk, Ayzc);
		for (i = 0; i < numbVar + numbLim; ++i){
			Ayzc[i] = Ayzc[i] + zk[i] - C[i];
		}


		//��y�����߂�G�������̍��ӂ̌v�Z
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, numbVar + numbLim, A, Z_diag, dy_L_t1);
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, numbVar + numbLim, dy_L_t1, X_diag, dy_L_t2);
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, numbLim, dy_L_t2, A_tra, LEFT);


		//��y�����߂�G�������̉E�ӂ̌v�Z
		//��1���̌v�Z
		calcMatrixMultiplication(numbLim, numbVar + numbLim, numbVar + numbLim, A, Z_diag, dy_R_1st_t1);
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, 1, X_diag, zk, Xkz);
		for (i = 0; i < numbVar + numbLim; ++i){
			dy_R_1st_t2[i] = mu - Xkz[i];
		}


		calcMatrixMultiplication(numbLim, numbVar + numbLim, 1, dy_L_t1, dy_R_1st_t2, dy_1st);

		//��3���̌v�Z (��2���� Ax-�� �Ȃ̂Ŋ��Ɍv�Z����)
		calcMatrixMultiplication(numbLim, numbVar + numbLim, 1, dy_L_t2/*���ӂ̌v�Z�𗬗p����*/, Ayzc, dy_3rd);

		//�E�ӂ̌v�Z
		for (i = 0; i < numbLim; ++i){
			RIGHT[i] = -dy_1st[i] - Axb[i] - dy_3rd[i];
		}

		//��y�̋���
		solvEquation(numbLim, LEFT, RIGHT, dy);


		//��z�����߂�
		calcMatrixMultiplication(numbVar + numbLim, numbLim, 1, A_tra, dy, dz_1st);
		//��z�̋���
		for (i = 0; i < numbVar + numbLim; ++i){
			dz[i] = -dz_1st[i] - Ayzc[i];
		}


		//��x�����߂�
		//��1���̌v�Z
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, numbVar + numbLim, Z_diag, X_diag, dx_1st_t1);
		calcMatrixMultiplication(numbVar + numbLim, numbVar + numbLim, 1, dx_1st_t1, dz, dx_1st_t2);

		//��2���̌v�Z
		for (i = 0; i < numbVar + numbLim; ++i){
			for (j = 0; j < numbVar + numbLim; ++j){
				if (i == j)dx_2nd[i] = mu * Z_diag[(numbVar + numbLim)*i + j];
			}
		}

		//��x�̋���
		for (i = 0; i < numbVar + numbLim; ++i){
			dx[i] = -dx_1st_t2[i] + dx_2nd[i] - xk[i];
		}


		//����+++++
		/*
		for (i = 0; i < 4; ++i){
			cout << dx[i] << " ";
		}
		cout << endl;//+++++
		*/

		//x,y,z�̒l���X�V����
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

	//��STEP5�G�v�Z���ʂ�Ԃ��i�Ԃ��p�̔z��ɓn���j

	for (i = 0; i < numbVar; ++i){
		optimumSolution[i] = xk[i];
		result += xk[i] + optimumSolution[i];
	}


	return result;
}
