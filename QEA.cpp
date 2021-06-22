#include "QEA.h"
#include "mt.h"
#include "pgm.h"
#include "fltr.h"

int main(void)
{
	int trial, generation;
	int i;
	double ave[MAX_GENERATION];
	FILE *fp;
	for(i=0; i<MAX_GENERATION; i++) ave[i] = 0.0;

	allocate_memory();
	load_learning_image();


	for(trial=1; trial<=NUM_TRIAL; trial++){

		printf("%d\n",trial);
		initialize();
		for(generation = 1; generation <= MAX_GENERATION; generation++){

			observation();
			evaluate();
			all_gene(trial, generation);
			fp=fopen("qea.txt","a");
			fprintf(fp,"%d\t%d\t%f\n",trial,generation-1,grobal_best_fitness);
			fclose(fp);
			ave[generation-1] += grobal_best_fitness;
			if(ELITE_FLAG) sort_by_fitness();
			if(SWAP_STRATEGY && !VARIABLE_LENGTH) swap_answer();
			update(generation);
			if(CROSSOVER) singlepoint_crossover();
			if(MUTATION) mutation();
			if(ELITE_FLAG) elite_strategy();
			}
		disp_result(trial);
		f_write(trial);
		save_image(trial, MAX_GENERATION);
		}

	fp=fopen("qea_ave.txt","a");
	for(i=0; i<MAX_GENERATION; i++) fprintf(fp,"%d\t%f\n",i,ave[i]/(double)NUM_TRIAL);
	fclose(fp);

	return 0;

}


void allocate_memory(void){

	qg = (q_gene *) calloc(NUM_GENE, sizeof(q_gene));

}

void initialize(void){
	int i,j;
	double pi_div_4;

	pi_div_4 = atan(1.0);
	grobal_best_fitness = 100000000;
	for(i=0; i<LENGTH_GENE * NUM_QBIT; i++) grobal_best_solution[i] = 0;

	for(i=0; i<NUM_GENE; i++){
		qg[i].best_fitness = 100000000;
		if(VARIABLE_LENGTH) qg[i].length = (genrand_int32() % LENGTH_GENE) + 1;
		else qg[i].length = LENGTH_GENE;
		for(j=0; j<qg[i].length * NUM_QBIT; j++){
			qg[i].b[j] = qg[i].p[j] = 0;
			qg[i].q[j] = pi_div_4;
			}
		}
}


void observation(void){
	int i,j;

	for(i=0; i<NUM_GENE; i++){
		for(j=0; j<qg[i].length * NUM_QBIT; j++){
			if(genrand_real1() <= pow(sin(qg[i].q[j]), 2.0)) qg[i].p[j] = 1;
			else qg[i].p[j] = 0;
			}
		}
}


void evaluate(void){
	int i, j, k;
	int x, y;
	int tmp_index;
	double tmp_min=100000000;
	int solution_decimal[LENGTH_GENE];

	// �R�X�g�̌v�Z
	for(i=0; i<NUM_GENE; i++){
		trans_2_to_10(qg[i].p, solution_decimal);

		for(y=0; y<height; y++){
			for(x=0; x<width; x++){
				image1[y][x] = original[y][x];
				}
			}
		for(j=0; j<qg[i].length; j++){
            /* �t�B���^�����O�����s */
            filtering(solution_decimal[j]);
			for(y=0; y<height; y++){
				for(x=0; x<width; x++){
					image1[y][x] = image2[y][x];
					}
				}
			}

		qg[i].fitness = 0.0;
		for(y=0; y<height; y++){
			for(x=0; x<width; x++){
				qg[i].fitness += (pow(((double)target[y][x] - image2[y][x]), 2.0) / (double)(width * height * MAX_BRIGHTNESS));
				}
			}

		// ���菬�����덷�����������Ƃ�
		if(qg[i].fitness < qg[i].best_fitness){
			qg[i].best_fitness = qg[i].fitness;
			memcpy(qg[i].b, qg[i].p, sizeof(qg[i].p));
			}
		}

	// �S�̂ōł��덷�̏������̂����߂�
	for(i=0; i<NUM_GENE; i++){
		if(qg[i].best_fitness < tmp_min){
			tmp_min = qg[i].best_fitness;
			tmp_index = i;
			}
		}
	if(tmp_min < grobal_best_fitness){
		grobal_best_fitness = tmp_min;
		grobal_best_length = qg[tmp_index].length;
		memcpy(grobal_best_solution, qg[tmp_index].b, sizeof(qg[tmp_index].b));
		}
}


void update(int generation){
	int i, j;

	for(i=0; i<NUM_GENE; i++){
		for(j=0; j<qg[i].length * NUM_QBIT; j++){
			if(UPDATE_MODE == 0) qg[i].q[j] += return_del_theta(i,j);
			else if(UPDATE_MODE ==1) qg[i].q[j] += return_del_theta_Improve_QEA(i,j,generation);
			}
		}
}


double return_del_theta(int i, int j){
	int index=0;
	double pi;

	pi = atan(1.0) * 4.0;

	if(qg[i].p[j] == 0){
		if(qg[i].b[j] == 0){
			if(qg[i].fitness < qg[i].best_fitness) index = 0;
			else index = 1;
		}
		else{
			if(qg[i].fitness < qg[i].best_fitness) index = 2;
			else index = 3;
		}
	}
	else{
		if(qg[i].b[j] == 0){
			if(qg[i].fitness < qg[i].best_fitness) index = 4;
			else index = 5;
		}
		else{
			if(qg[i].fitness < qg[i].best_fitness) index = 6;
			else index = 7;
		}
	}

	return del_theta_list[index] * pi;
}


double return_del_theta_Improve_QEA(int i, int j, int generation){
	double sign;
	double pi;
	double k=0.5;

	pi = atan(1.0) * 4.0;

	if(qg[i].p[j] == 0){
		if(qg[i].b[j] == 0){
			if(qg[i].fitness < qg[i].best_fitness) sign = 0.0;
			else sign = 0.0;
			}
		else{
			if(qg[i].fitness < qg[i].best_fitness){
				if(cos(qg[i].q[j]) * sin(qg[i].q[j]) > 0.0) sign = 1.0;
				else if(cos(qg[i].q[j]) * sin(qg[i].q[j]) < 0.0) sign = -1.0;
				else if(cos(qg[i].q[j]) == 0.0) sign = 0.0;
				else sign = 1.0;
				}
			else{
				sign=0.0;
				}
			}
		}
	else{
		if(qg[i].b[j] == 0){
			if(qg[i].fitness < qg[i].best_fitness){
				if(cos(qg[i].q[j]) * sin(qg[i].q[j]) > 0.0) sign = -1.0;
				else if(cos(qg[i].q[j]) * sin(qg[i].q[j]) < 0.0) sign = 1.0;
				else if(cos(qg[i].q[j]) == 0.0) sign = 1.0;
				else sign = 0.0;
				}
			else{
				sign=0.0;
				}
			}
		else{
			if(qg[i].fitness < qg[i].best_fitness) sign = 0.0;
			else sign = 0.0;
			}
		}

	return sign * 0.01 * pi * (1 - (k * (double)generation/(MAX_GENERATION + 1)));
}


void swap_answer(void){
	int i;
	int flag[NUM_GENE];
	int swap1, swap2;

	for(i=0; i<NUM_GENE; i++) flag[i] = 0;

	while(1){
		swap1 = genrand_int32() % NUM_GENE;
		if(!flag[swap1]){
			while(1){
				swap2 = genrand_int32() % NUM_GENE;
				if(swap1 != swap2 && !flag[swap2]){
					my_swap(swap1, swap2);
					flag[swap1] = flag[swap2] = 1;
					break;
					}
				}
			}

		for(i=0; i<NUM_GENE && flag[i]; i++);
		if(i == NUM_GENE) break;

		}
}


void my_swap(int swap1, int swap2){
	double bf_tmp;
	int b_tmp[LENGTH_GENE * NUM_QBIT];

	bf_tmp = qg[swap1].best_fitness;
	qg[swap1].best_fitness = qg[swap2].best_fitness;
	qg[swap2].best_fitness = bf_tmp;

	memcpy(b_tmp, qg[swap1].b, sizeof(qg[swap1].b));
	memcpy(qg[swap1].b, qg[swap2].b, sizeof(qg[swap2].b));
	memcpy(qg[swap2].b, b_tmp, sizeof(b_tmp));

}


void disp_result(int trial){
	int i;
	int solution_decimal[LENGTH_GENE];

	printf("trial=%d best_fitness=%f\n", trial, grobal_best_fitness);
	trans_2_to_10(grobal_best_solution, solution_decimal);
	for(i=0; i<grobal_best_length; i++) if(solution_decimal[i]>0 && solution_decimal[i]<17) printf("%d ", solution_decimal[i]);
	putchar('\n');

}


void f_write(int trial){
	int i;
	int solution_decimal[LENGTH_GENE];
	FILE *fp;

	if((fp = fopen(RESULT_FILE_NAME, "a")) == NULL){
		printf("開けませーン\n");	
		exit(1);
		}

	fprintf(fp, "%d\t%f\t", trial, grobal_best_fitness);
	trans_2_to_10(grobal_best_solution, solution_decimal);
	for(i=0; i<grobal_best_length; i++) if(solution_decimal[i]>0 && solution_decimal[i]<17) fprintf(fp, "%d ", solution_decimal[i]);
	fputc('\n', fp);

	fclose(fp);

}

#if(GRAY_CODE)
void trans_2_to_10(int *sol_binary, int *sol_decimal){
	int i, j, k;
	int val, n, parity;

	for(i=0; i<LENGTH_GENE; i++){
		for(j=NUM_QBIT*(i+1)-1, n=1, val=0; j>=NUM_QBIT*i; j--,n*=2){
			for(k=NUM_QBIT*i, parity=0; k<=j; k++) parity ^= sol_binary[k];
			val += (n * parity);
			}
		sol_decimal[i] = val;
		}
}
#else
void trans_2_to_10(int *sol_binary, int *sol_decimal){
	int i, j;
	int val, n;

	for(i=0; i<LENGTH_GENE; i++){
		for(j=NUM_QBIT*(i+1)-1, n=1, val=0; j>=NUM_QBIT*i; j--,n*=2) val += (n * sol_binary[j]);
		sol_decimal[i] = val;
		}
}
#endif

void singlepoint_crossover(void){
	int i;
	int num[NUM_GENE];
	int tmp;

	// �����_���ɑg�ݍ��킹��
	for(i=0; i<NUM_GENE; i++) num[i] = -1;
	for(i=0; i<NUM_GENE; i++){
		while(1){
			tmp = genrand_int32() % NUM_GENE;
			if(num[tmp] == -1){
				num[tmp] = i;
				break;
				}
			}
		}

	for(i=0; i<NUM_GENE; i+=2){
		if(genrand_real1() <= CROSSOVER_RATE){
			if(!VARIABLE_LENGTH) execute_crossover(num[i], num[i+1]);
			else if(qg[num[i]].length > 1 && qg[num[i+1]].length > 1) execute_crossover2(num[i], num[i+1]);
			}
		}

}


void execute_crossover(int num1, int num2){
	int i;
	int point;
	q_gene tmp;

	point = (genrand_int32() % (LENGTH_GENE - 1)) + 1; // 1�`LENGTH_GENE-1�͈̔͂Ō����ʒu�����߂�
	tmp = qg[num2];
	for(i=point * NUM_QBIT; i<LENGTH_GENE * NUM_QBIT; i++){
		qg[num2].q[i] = qg[num1].q[i];
		qg[num1].q[i] = tmp.q[i];
		}
}


void execute_crossover2(int num1, int num2){
	int i, j;
	int point1, point2;
	int tmp_length;
	q_gene tmp;

	point1 = (genrand_int32() % (qg[num1].length - 1)) + 1; // 1�`length-1�͈̔͂Ō����ʒu�����߂�
	point2 = (genrand_int32() % (qg[num2].length  - 1)) + 1; // 1�`length-1�͈̔͂Ō����ʒu�����߂�

	tmp_length = qg[num1].length;
	qg[num1].length = point1 + qg[num2].length - point2;
	qg[num2].length = point2 + tmp_length - point1;
	if(qg[num1].length > LENGTH_GENE) qg[num1].length = LENGTH_GENE;
	if(qg[num2].length > LENGTH_GENE) qg[num2].length = LENGTH_GENE;
	tmp = qg[num2];
	for(i=point2 * NUM_QBIT, j=0; i<qg[num2].length * NUM_QBIT; i++, j++)	qg[num2].q[i] = qg[num1].q[point1 * NUM_QBIT + j];
	for(i=point1 * NUM_QBIT, j=0; i<qg[num1].length * NUM_QBIT; i++, j++)	qg[num1].q[i] = tmp.q[point2 * NUM_QBIT + j];

}


void mutation(void){
	int i;
	int point;

	for(i=0; i<NUM_GENE; i++){
		if(genrand_real1() <= MUTATION_RATE){
			point = genrand_int32() % (qg[i].length * NUM_QBIT); // 0�`LENGTH_GENE*NUM_QBIT-1�͈̔͂œˑR�ψق̈ʒu�����߂�
			qg[i].q[point] = atan(1.0) * 2.0 - qg[i].q[point]; // �m���𔽓]������
			}
		}
}


void sort_by_fitness(void){
	int i, j, tmp_i;
	int sorted_index[NUM_GENE];
	q_gene *tmp;

	tmp = (q_gene *) calloc(NUM_GENE, sizeof(q_gene));

    for(i=0; i<NUM_GENE; i++) sorted_index[i] = i;

    for(i=0; i<NUM_GENE; i++){
      for(j=i+1; j<NUM_GENE; j++){
        if(qg[sorted_index[i]].fitness < qg[sorted_index[j]].fitness){
          tmp_i = sorted_index[i];
          sorted_index[i] = sorted_index[j];
          sorted_index[j] = tmp_i;
				}
			}
		}

    memcpy(tmp, qg, _msize(qg));
	memset(qg, 0, _msize(qg));
    for(i=0; i<NUM_GENE; i++) qg[i] = tmp[sorted_index[i]];

	elite = qg[NUM_GENE-1]; // �G���[�g�̂��ۑ�
	free(tmp);
}


void elite_strategy(void){
	qg[0] = elite;
}


void load_learning_image(void){
	int x, y;

	// ���͉摜�̓ǂݍ���
	load_image_data(INPUT_FNAME);
	for(y=0; y<height; y++){
		for(x=0; x<width; x++){
			original[y][x] = image1[y][x];
			}
		}
	// �ڕW�摜�̓ǂݍ���
	load_image_data(TARGET_FNAME);
	for(y=0; y<height; y++){
		for(x=0; x<width; x++){
			target[y][x] = image1[y][x];
			}
		}
}


void save_image(int trial, int generation){
	int i, j, x, y;
	int solution_decimal[LENGTH_GENE];
	char fname[100];

	trans_2_to_10(grobal_best_solution, solution_decimal);
	for(y=0; y<height; y++){
		for(x=0; x<width; x++){
			image1[y][x] = original[y][x];
			}
		}
	for(j=0; j<grobal_best_length; j++){
        /* �t�B���^�����O�����s */
        filtering(solution_decimal[j]);
		for(y=0; y<height; y++){
			for(x=0; x<width; x++){
				image1[y][x] = image2[y][x];
				}
			}
		}

	sprintf(fname, "%s_%d_%d.pgm", OUTPUT_FNAME, trial, generation);
	save_image_data(fname);

}


void all_gene(int trial, int generation){
	int i, j;
	FILE *fp;
	int solution_decimal;

	fp=fopen("qea_all.txt","a");
	for(i=0; i<NUM_GENE; i++){
		//trans_2_to_10(qg[i].p, &solution_decimal);
		fprintf(fp,"%d\t%d\t%d\t%f\t",trial, generation-1, i, qg[i].fitness);
		for(j=0; j<LENGTH_GENE * NUM_QBIT; j++) fprintf(fp, "%d", qg[i].p[j]);
		fprintf(fp, "\n");
	}
	fclose(fp);

}
