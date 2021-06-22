#ifndef QEA_H__
#define QEA_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#define SWAP_STRATEGY 0 // �Ό����헪���s�����ǂ���(Yes: 1, No: 0)�C���F�̒����Œ�̏ꍇ�̂ݎg�p��
#define GRAY_CODE 1 // 2�i10�i�ϊ����O���C�R�[�h�����ŕϊ����邩�ǂ���(Yes: 1, No: 0)

#define CROSSOVER 1 // 1�_�������s�����ǂ���(Yes: 1, No: 0)
#define MUTATION 1 // �ˑR�ψق��s�����ǂ���(Yes: 1, No: 0)
#define ELITE_FLAG 0 // �G���[�g�ۑ����s�����ǂ���(Yes: 1, No: 0)
#define UPDATE_MODE 1 // �A�b�v�f�[�g�̕��@(0: QEA�C1: QEA��)

#define NUM_GENE 10 // �̐�(�Ό����헪������̓s����������]�܂���)
#define VARIABLE_LENGTH 0 // �l�X�Ȑ��F�̒��̌̂𐶐�(Yes: 1�CNo:0)�CNo�̏ꍇ�͍ő���F�̒��ŌŒ�
#define NUM_QBIT 5 // 1�̈�`�q��������l�͈̔͂�����(0�`15�������ꍇ��4bit�K�v�Ȃ̂�NUM_QBIT=4�ƂȂ�)
#define LENGTH_GENE 20  // �ő���F�̒�
#define MAX_GENERATION 100 // �ő吢�㐔
#define CROSSOVER_RATE 0.8 // �����̋N����m��
#define MUTATION_RATE 0.01 // �ˑR�ψق̋N����m��

#define NUM_TRIAL 10 // ���s��
#define RESULT_FILE_NAME  "result.txt" // ���ʂ��o�͂���t�@�C����

double del_theta_list[8] = {0.0, 0.0, 0.01, 0.0, -0.01, 0.0, 0.0, 0.0}; // ���̒l�~�΂����ۂ̕ω���

/* �摜�����p�ϐ��錾*/
#define INPUT_FNAME "input.pgm" // ���̓t�@�C����
#define TARGET_FNAME "teach.pgm" // �ڕW�t�@�C����
#define OUTPUT_FNAME "output" // �o�̓t�@�C����
#define MAX_IMAGESIZE   1000 /* �z�肷��c�E���̍ő��f��    */
#define MAX_BRIGHTNESS  255 /* �z�肷��ő�K���l            */
#define GRAYLEVEL       256 /* �z�肷��K����(=�ő�K���l+1) */
#define MAX_FILENAME    256 /* �z�肷��t�@�C�����̍ő咷    */
#define MAX_BUFFERSIZE  256 /* ���p����o�b�t�@�ő咷        */
unsigned char image1[MAX_IMAGESIZE][MAX_IMAGESIZE], image2[MAX_IMAGESIZE][MAX_IMAGESIZE]; // �摜�p�z��1�C�摜�p�z��2
int width, height; // �摜�̉���f���C�c��f��
unsigned char original[MAX_IMAGESIZE][MAX_IMAGESIZE]; // ���摜
unsigned char target[MAX_IMAGESIZE][MAX_IMAGESIZE]; // �ڕW�摜

// ��`�q���
struct q_gene{
	int length; // ���F�̂̒���
	int b[LENGTH_GENE * NUM_QBIT]; // ���݂܂ł̍œK��
	int p[LENGTH_GENE * NUM_QBIT]; // ���̒l(0�܂���1)
	double q[LENGTH_GENE * NUM_QBIT]; // �m���U��(�ʑ��l)
	double fitness; // ���݂̓K���x�̒l
	double best_fitness; // ���݂܂łɍł����������K���x�̒l
};

q_gene *qg;
q_gene elite;
double grobal_best_fitness; // �S�̒��ɂ�����K���x�̍ō��l
int grobal_best_solution[LENGTH_GENE * NUM_QBIT]; // ����̉��̒l
int grobal_best_length; // ����̉��̒���


void allocate_memory(void); // �������̊��蓖��
void initialize(void); // ������
void observation(void); // �ϑ�
void evaluate(void); // �]��
void swap_answer(void); // �Ό����헪���s���֐�
void my_swap(int swap1, int swap2); // �ŗǉ�������������
void update(int generation); // �Ƃ̍X�V
double return_del_theta(int i, int j); // ���Ƃ̒l��Ԃ��֐�
double return_del_theta_Improve_QEA(int i, int j, int generation); // ���Ƃ̒l��Ԃ��֐�(Improve_QEA)
void disp_result(int trial);  // ���ʂ���ʂɏo��
void f_write(int trial); // ���ʂ��t�@�C���ɏ�������
void singlepoint_crossover(void); // 1�_����
void execute_crossover(int num1, int num2); // ���������s����֐�
void execute_crossover2(int num1, int num2); // ���������s����֐�(���F�̒����ς̏ꍇ)
void mutation(void); // �ˑR�ψ�
void sort_by_fitness(void); // �K���x�̍������Ƀ\�[�g
void elite_strategy(void); // �G���[�g�ۑ�
void trans_2_to_10(int *sol_binary, int *sol_decimal); // 2�i����10�i�֕ϊ�
void load_learning_image(void); // �w�K�p�摜��ǂݍ���
void save_image(int trial, int generation); // �摜��ۑ�����
void all_gene(int trial, int generation);


#endif


