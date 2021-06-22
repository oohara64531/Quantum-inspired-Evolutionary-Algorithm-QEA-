#pragma warning(disable: 4996)
#include "QEA.h"

/* �֐��̃v���g�^�C�v�錾 */
void load_image_data(char *fname); /* �摜�ǂݍ��ݗp�֐� */
void save_image_data(char *fname); /* �摜�������ݗp�֐� */

void load_image_data(const char *fname)
/* pgm �摜�C����f���C�c��f���̃f�[�^���t�@�C������ǂݍ��݁C*/
/* image1[ ][ ]�Cwidth�Cheight �ɂ��ꂼ��������D            */
{
    char buffer[MAX_BUFFERSIZE];  /* �f�[�^�ǂݍ��ݗp��ƕϐ�  */
    FILE *fp;                     /* �t�@�C���|�C���^          */
    int max_gray;                 /* �ő�K���l                */
    int x, y;                     /* ���[�v�ϐ�                */

    /* ���̓t�@�C���̃I�[�v�� */
    fp = fopen( fname, "rb" );
    if ( NULL == fp ){
        printf("���̖��O�̃t�@�C���͑��݂��܂���D\n");
        exit(1);
    }
    /* �t�@�C���^�C�v(=P5)�̊m�F */
    fgets( buffer, MAX_BUFFERSIZE, fp );
    if ( buffer[0] != 'P' || buffer[1] != '5' ){
        printf("�t�@�C���̃t�H�[�}�b�g�� P5 �Ƃ͈قȂ�܂��D\n");
        exit(1);
    }
    /* width, height �̑���i#����n�܂�R�����g�͓ǂݔ�΂��j */
    width = 0;
    height = 0;
    while ( width == 0 || height == 0 ){
        fgets( buffer, MAX_BUFFERSIZE, fp );
        if ( buffer[0] != '#' ){
            sscanf( buffer, "%d %d", &width, &height );
        }
    }
    /* max_gray �̑���i#����n�܂�R�����g�͓ǂݔ�΂��j */
    max_gray = 0;
    while ( max_gray == 0 ){
        fgets( buffer, MAX_BUFFERSIZE, fp );
        if ( buffer[0] != '#' ){
            sscanf( buffer, "%d", &max_gray );
        }
    }
    /* �p�����[�^�̉�ʂւ̕\�� */
    if ( width > MAX_IMAGESIZE || height > MAX_IMAGESIZE ){
        printf("�z��l %d x %d �𒴂��Ă��܂��D\n", 
            MAX_IMAGESIZE, MAX_IMAGESIZE);
        printf("�������������ȉ摜���g���ĉ������D\n");
        exit(1);
    }
    if ( max_gray != MAX_BRIGHTNESS ){
        printf("�ő�K���l���s�K�؂ł��D\n");
        exit(1);
    }
    /* �摜�f�[�^��ǂݍ���ŉ摜�p�z��ɑ������ */
    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            image1[y][x] = (unsigned char)fgetc( fp );
        }
    }
    fclose(fp);
}

void save_image_data(char *fname)
/* image2[ ][ ], width, height �̃f�[�^���C���ꂼ�� pgm �摜�C*/
/* ����f���C�c��f���Ƃ��ăt�@�C���ɕۑ�����D               */
{
    FILE *fp;                     /* �t�@�C���|�C���^         */
    int x, y;                     /* ���[�v�ϐ�               */

    /* �o�̓t�@�C���̃I�[�v�� */
    fp = fopen(fname, "wb");
    /* �t�@�C�����ʎq "P5" ��擪�ɏo�͂��� */
    fputs( "P5\n", fp );
    /* # �Ŏn�܂�R�����g�s�i�ȗ��\�j */
    fputs( "# Created by Image Processing\n", fp );
    /* �摜�̉����C�c���̏o�� */
    fprintf( fp, "%d %d\n", width, height );
    /* �ő�K���l�̏o�� */
    fprintf( fp, "%d\n", MAX_BRIGHTNESS );
    /* �摜�f�[�^�̏o�� */
    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            fputc( image2[y][x], fp );
        }
    }
    fclose(fp);
}
