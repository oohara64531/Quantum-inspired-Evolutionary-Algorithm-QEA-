#pragma warning(disable: 4996)
#include "QEA.h"

/* 関数のプロトタイプ宣言 */
void load_image_data(char *fname); /* 画像読み込み用関数 */
void save_image_data(char *fname); /* 画像書き込み用関数 */

void load_image_data(const char *fname)
/* pgm 画像，横画素数，縦画素数のデータをファイルから読み込み，*/
/* image1[ ][ ]，width，height にそれぞれ代入する．            */
{
    char buffer[MAX_BUFFERSIZE];  /* データ読み込み用作業変数  */
    FILE *fp;                     /* ファイルポインタ          */
    int max_gray;                 /* 最大階調値                */
    int x, y;                     /* ループ変数                */

    /* 入力ファイルのオープン */
    fp = fopen( fname, "rb" );
    if ( NULL == fp ){
        printf("その名前のファイルは存在しません．\n");
        exit(1);
    }
    /* ファイルタイプ(=P5)の確認 */
    fgets( buffer, MAX_BUFFERSIZE, fp );
    if ( buffer[0] != 'P' || buffer[1] != '5' ){
        printf("ファイルのフォーマットが P5 とは異なります．\n");
        exit(1);
    }
    /* width, height の代入（#から始まるコメントは読み飛ばす） */
    width = 0;
    height = 0;
    while ( width == 0 || height == 0 ){
        fgets( buffer, MAX_BUFFERSIZE, fp );
        if ( buffer[0] != '#' ){
            sscanf( buffer, "%d %d", &width, &height );
        }
    }
    /* max_gray の代入（#から始まるコメントは読み飛ばす） */
    max_gray = 0;
    while ( max_gray == 0 ){
        fgets( buffer, MAX_BUFFERSIZE, fp );
        if ( buffer[0] != '#' ){
            sscanf( buffer, "%d", &max_gray );
        }
    }
    /* パラメータの画面への表示 */
    if ( width > MAX_IMAGESIZE || height > MAX_IMAGESIZE ){
        printf("想定値 %d x %d を超えています．\n", 
            MAX_IMAGESIZE, MAX_IMAGESIZE);
        printf("もう少し小さな画像を使って下さい．\n");
        exit(1);
    }
    if ( max_gray != MAX_BRIGHTNESS ){
        printf("最大階調値が不適切です．\n");
        exit(1);
    }
    /* 画像データを読み込んで画像用配列に代入する */
    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            image1[y][x] = (unsigned char)fgetc( fp );
        }
    }
    fclose(fp);
}

void save_image_data(char *fname)
/* image2[ ][ ], width, height のデータを，それぞれ pgm 画像，*/
/* 横画素数，縦画素数としてファイルに保存する．               */
{
    FILE *fp;                     /* ファイルポインタ         */
    int x, y;                     /* ループ変数               */

    /* 出力ファイルのオープン */
    fp = fopen(fname, "wb");
    /* ファイル識別子 "P5" を先頭に出力する */
    fputs( "P5\n", fp );
    /* # で始まるコメント行（省略可能） */
    fputs( "# Created by Image Processing\n", fp );
    /* 画像の横幅，縦幅の出力 */
    fprintf( fp, "%d %d\n", width, height );
    /* 最大階調値の出力 */
    fprintf( fp, "%d\n", MAX_BRIGHTNESS );
    /* 画像データの出力 */
    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            fputc( image2[y][x], fp );
        }
    }
    fclose(fp);
}
