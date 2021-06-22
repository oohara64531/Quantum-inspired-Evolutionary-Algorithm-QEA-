#ifndef QEA_H__
#define QEA_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#define SWAP_STRATEGY 0 // 対交換戦略を行うかどうか(Yes: 1, No: 0)，染色体長が固定の場合のみ使用可
#define GRAY_CODE 1 // 2進10進変換をグレイコード方式で変換するかどうか(Yes: 1, No: 0)

#define CROSSOVER 1 // 1点交叉を行うかどうか(Yes: 1, No: 0)
#define MUTATION 1 // 突然変異を行うかどうか(Yes: 1, No: 0)
#define ELITE_FLAG 0 // エリート保存を行うかどうか(Yes: 1, No: 0)
#define UPDATE_MODE 1 // アップデートの方法(0: QEA，1: QEA改)

#define NUM_GENE 10 // 個体数(対交換戦略や交叉の都合上偶数が望ましい)
#define VARIABLE_LENGTH 0 // 様々な染色体長の個体を生成(Yes: 1，No:0)，Noの場合は最大染色体長で固定
#define NUM_QBIT 5 // 1つの遺伝子が扱える値の範囲を決定(0～15を扱う場合は4bit必要なのでNUM_QBIT=4となる)
#define LENGTH_GENE 20  // 最大染色体長
#define MAX_GENERATION 100 // 最大世代数
#define CROSSOVER_RATE 0.8 // 交叉の起こる確率
#define MUTATION_RATE 0.01 // 突然変異の起こる確率

#define NUM_TRIAL 10 // 試行回数
#define RESULT_FILE_NAME  "result.txt" // 結果を出力するファイル名

double del_theta_list[8] = {0.0, 0.0, 0.01, 0.0, -0.01, 0.0, 0.0, 0.0}; // 左の値×πが実際の変化量

/* 画像処理用変数宣言*/
#define INPUT_FNAME "input.pgm" // 入力ファイル名
#define TARGET_FNAME "teach.pgm" // 目標ファイル名
#define OUTPUT_FNAME "output" // 出力ファイル名
#define MAX_IMAGESIZE   1000 /* 想定する縦・横の最大画素数    */
#define MAX_BRIGHTNESS  255 /* 想定する最大階調値            */
#define GRAYLEVEL       256 /* 想定する階調数(=最大階調値+1) */
#define MAX_FILENAME    256 /* 想定するファイル名の最大長    */
#define MAX_BUFFERSIZE  256 /* 利用するバッファ最大長        */
unsigned char image1[MAX_IMAGESIZE][MAX_IMAGESIZE], image2[MAX_IMAGESIZE][MAX_IMAGESIZE]; // 画像用配列1，画像用配列2
int width, height; // 画像の横画素数，縦画素数
unsigned char original[MAX_IMAGESIZE][MAX_IMAGESIZE]; // 元画像
unsigned char target[MAX_IMAGESIZE][MAX_IMAGESIZE]; // 目標画像

// 遺伝子情報
struct q_gene{
	int length; // 染色体の長さ
	int b[LENGTH_GENE * NUM_QBIT]; // 現在までの最適解
	int p[LENGTH_GENE * NUM_QBIT]; // 解の値(0または1)
	double q[LENGTH_GENE * NUM_QBIT]; // 確率振幅(位相値)
	double fitness; // 現在の適応度の値
	double best_fitness; // 現在までに最も高かった適応度の値
};

q_gene *qg;
q_gene elite;
double grobal_best_fitness; // 全個体中における適応度の最高値
int grobal_best_solution[LENGTH_GENE * NUM_QBIT]; // 同上の解の値
int grobal_best_length; // 同上の解の長さ


void allocate_memory(void); // メモリの割り当て
void initialize(void); // 初期化
void observation(void); // 観測
void evaluate(void); // 評価
void swap_answer(void); // 対交換戦略を行う関数
void my_swap(int swap1, int swap2); // 最良怪情報を交換する
void update(int generation); // θの更新
double return_del_theta(int i, int j); // ⊿θの値を返す関数
double return_del_theta_Improve_QEA(int i, int j, int generation); // ⊿θの値を返す関数(Improve_QEA)
void disp_result(int trial);  // 結果を画面に出力
void f_write(int trial); // 結果をファイルに書き込む
void singlepoint_crossover(void); // 1点交叉
void execute_crossover(int num1, int num2); // 交叉を実行する関数
void execute_crossover2(int num1, int num2); // 交叉を実行する関数(染色体長が可変の場合)
void mutation(void); // 突然変異
void sort_by_fitness(void); // 適応度の高い順にソート
void elite_strategy(void); // エリート保存
void trans_2_to_10(int *sol_binary, int *sol_decimal); // 2進から10進へ変換
void load_learning_image(void); // 学習用画像を読み込む
void save_image(int trial, int generation); // 画像を保存する
void all_gene(int trial, int generation);


#endif


