/*                                                       */
/* 「進化的画像処理」昭晃堂　掲載プログラム              */
/*                                                       */
/* 画像のフィルタリングプログラムのヘッダファイル fltr.h */

#define MAX_FILTER_NUMBER  16      /* フィルタの最大番号 */

/* 関数のプロトタイプ宣言 */
void filtering( int number );
/* number:フィルタの種類，img_number:原画像の番号 */

unsigned char _image1( int y, int x )
/* image1[y][x] の値を返す                               */
/* 画像の外周の画素についても3×3近傍のフィルタリングを  */
/* 実行するためのものである．                            */
/* (x,y) が必ず範囲内であることが既知な場合は使用しない．*/
{
    if ( x < 0 ) x = 0; 
    else if ( x >= width ) x = width - 1;
    if ( y < 0 ) y = 0;
    else if ( y >= height ) y = height - 1;
    return image1[y][x];
}

void _thresholding( int n )
/* 固定しきい値で2値化 (0またはMAX_BRIGHTNESS(=255)) する */
{
    int threshold; /* 2値化のしきい値 */
    int x, y;      /* ループ変数       */

    switch( n ){
        case 1 : threshold =  50;    break;
        case 2 : threshold = 100;    break;
        case 3 : threshold = 150;    break;
        case 4 : threshold = 200;    break;
    }
    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            if ( image1[y][x] >= threshold )
                image2[y][x] = MAX_BRIGHTNESS; 
            else image2[y][x] = 0;
        }
    }
}

void _mean( )
/* 3×3近傍内の平均値に置き換えて平滑化する */
{
    int x, y, i, j; /* ループ変数   */
    int sum;        /* 階調値の合計 */

    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            sum = 0;
            for ( i = -1; i <= 1; i ++ ){
                for ( j = -1; j <= 1; j ++ ){
                    sum = sum + _image1(y+i, x+j);
                }
            }
            image2[y][x] = (unsigned char)( sum / 9.0 );
        }
    }
}

void _min( )
/* 3×3近傍内の最小値に置き換える */
{
    int x, y, i, j; /* ループ変数       */
    int min;        /* 近傍の最小階調値 */

    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            min = MAX_BRIGHTNESS;
            for ( i = -1; i <= 1; i ++ ){
                for ( j = -1; j <= 1; j ++ ){
                    if ( _image1( y+i, x+j ) < min ){
                        min = _image1( y+i, x+j );
                    }
                }
            }
            image2[y][x] = min;
        }
    }
}

void _max( )
/* 3×3近傍内の最大値に置き換える */
{
    int x, y, i, j; /* ループ変数       */
    int max;        /* 近傍の最大階調値 */

    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            max = 0;
            for ( i = -1; i <= 1; i ++ ){
                for ( j = -1; j <= 1; j ++ ){
                    if ( _image1( y+i, x+j ) > max ){
                        max = _image1( y+i, x+j );
                    }
                }
            }
            image2[y][x] = max;
        }
    }
}

void _sobel( int n )
/* Sobelフィルタによる1次微分値にする */
/* n=1：横方向， n=2：縦方向， n=3：微分の大きさ */
{
    int x, y; /* ループ変数 */
    double horizontal, vertical; /* 横方向微分値，縦方向微分値 */
    int strength; /* 微分の大きさ */

    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            horizontal = - _image1(y-1,x-1) - 2 * _image1(y,x-1)
                      - _image1(y+1,x-1) + _image1(y-1,x+1) 
                       + 2 * _image1(y,x+1) + _image1(y+1,x+1);
            vertical   = - _image1(y-1,x-1) - 2 * _image1(y-1,x) 
                      - _image1(y-1,x+1) + _image1(y+1,x-1)
                      + 2 * _image1(y+1,x) + _image1(y+1,x+1);
            switch ( n ){
                case 1: strength = 
                            (int)( horizontal );
                    break;
                case 2: strength =
                            (int)( vertical );
                    break;
                case 3: strength = (int)( sqrt( horizontal * horizontal + 
                                           vertical * vertical ) );
            }
            if ( strength < 0 ) strength = 0;
            else if ( strength > MAX_BRIGHTNESS )
                strength = MAX_BRIGHTNESS;
            image2[y][x] = strength;
        }
    }
}

void _lightedge( )
/* ラプラシアンによる2次微分値にする */
{
    int x, y;     /* ループ変数   */
    int strength; /* 微分の大きさ */

    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            strength = (int)(
                 - _image1(y-1,x-1) - _image1(y-1,x) - _image1(y-1,x+1)
                 - _image1(y  ,x-1) + 8 * _image1(y, x) - _image1(y,x+1)
                 - _image1(y+1,x-1) - _image1(y+1,x) - _image1(y+1,x+1) );
            if ( strength < 0 ) strength = 0;
            else if ( strength > MAX_BRIGHTNESS )
                strength = MAX_BRIGHTNESS;
            image2[y][x] = strength;
        }
    }
}

void _darkedge( )
/* ラプラシアンによる2次微分値 + MAX_BRIGHTNESS にする */
{
    int x, y;     /* ループ変数   */
    int strength; /* 微分の大きさ */

    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            strength = (int)(
                 - _image1(y-1,x-1) - _image1(y-1,x) - _image1(y-1,x+1)
                 - _image1(y  ,x-1) + 8 * _image1(y, x) - _image1(y,x+1)
                 - _image1(y+1,x-1) - _image1(y+1,x) - _image1(y+1,x+1) );
            strength = strength + MAX_BRIGHTNESS;
            if ( strength < 0 ) strength = 0;
            else if ( strength > MAX_BRIGHTNESS )
                strength = MAX_BRIGHTNESS;
            image2[y][x] = strength;
        }
    }
}

void _lightpixel( )
/* 平均階調より暗い画素 → 0 */
{
    int x, y;       /* ループ変数 */
    double average; /* 平均階調値 */

    average = 0.0;
    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            average = average + (double)image1[y][x] / width / height;
        }
    }
    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            if ( image1[y][x] < (int)average ) image2[y][x] = 0;
            else image2[y][x] = image1[y][x];
        }
    }
}

void _darkpixel( )
/* 平均階調より明るい画素 → MAX_BRIGHTNESS */
{
    int x, y;       /* ループ変数 */
    double average; /* 平均階調値 */

    average = 0.0;
    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            average = average + (double)image1[y][x] / width / height;
        }
    }
    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            if ( image1[y][x] >= (int)average )
                image2[y][x] = MAX_BRIGHTNESS;
            else image2[y][x] = image1[y][x];
        }
    }
}

void _inversion( )
/* 反転処理 */
{
    int x, y; /* ループ変数 */

    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            image2[y][x] = MAX_BRIGHTNESS - image1[y][x];
        }
    }
}

void _subtract(void)
/* 原画像（=original[img_num-1][][] ）との差分の絶対値にする */
{
    int x, y; /* ループ変数 */

    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            image2[y][x] = abs( original[y][x] - image1[y][x] );
        }
    }
}

void copy(void){
	int x, y;

	for(y=0; y<height; y++){
		for(x=0; x<width; x++){
			image2[y][x] = image1[y][x];
			}
		}

}

void filtering( int number )
/* 大域変数で宣言されている横width 縦height の画像 データ */
/* image1, image2 を用いて画像のフィルタリングを実行するプログラム */
/* 常に image1 から image2 への変換を行うものとする． */
/* int number : 適用するフィルタの番号(1,2,...,MAX_FILTER_NUMBER) */
/* int img_number : 対象とする原画像の番号(1,2,...,MAX_LEARN) */
{
    switch( number ){
        case  1: case  2: case  3: case  4: _thresholding( number );
            break;
        case  5: _mean( );  break;
        case  6: _min( );  break;
        case  7: _max( );  break;
        case  8: case  9: case 10: _sobel( number-7 );  break;
        case 11: _lightedge( );  break;
        case 12: _darkedge( );  break;
        case 13: _lightpixel( );  break;
        case 14: _darkpixel( );  break;
        case 15: _inversion( );  break;
        case 16: _subtract();  break;
		default: copy();
    }
}
