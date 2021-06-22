/*                                                       */
/* �u�i���I�摜�����v���W���@�f�ڃv���O����              */
/*                                                       */
/* �摜�̃t�B���^�����O�v���O�����̃w�b�_�t�@�C�� fltr.h */

#define MAX_FILTER_NUMBER  16      /* �t�B���^�̍ő�ԍ� */

/* �֐��̃v���g�^�C�v�錾 */
void filtering( int number );
/* number:�t�B���^�̎�ށCimg_number:���摜�̔ԍ� */

unsigned char _image1( int y, int x )
/* image1[y][x] �̒l��Ԃ�                               */
/* �摜�̊O���̉�f�ɂ��Ă�3�~3�ߖT�̃t�B���^�����O��  */
/* ���s���邽�߂̂��̂ł���D                            */
/* (x,y) ���K���͈͓��ł��邱�Ƃ����m�ȏꍇ�͎g�p���Ȃ��D*/
{
    if ( x < 0 ) x = 0; 
    else if ( x >= width ) x = width - 1;
    if ( y < 0 ) y = 0;
    else if ( y >= height ) y = height - 1;
    return image1[y][x];
}

void _thresholding( int n )
/* �Œ肵�����l��2�l�� (0�܂���MAX_BRIGHTNESS(=255)) ���� */
{
    int threshold; /* 2�l���̂������l */
    int x, y;      /* ���[�v�ϐ�       */

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
/* 3�~3�ߖT���̕��ϒl�ɒu�������ĕ��������� */
{
    int x, y, i, j; /* ���[�v�ϐ�   */
    int sum;        /* �K���l�̍��v */

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
/* 3�~3�ߖT���̍ŏ��l�ɒu�������� */
{
    int x, y, i, j; /* ���[�v�ϐ�       */
    int min;        /* �ߖT�̍ŏ��K���l */

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
/* 3�~3�ߖT���̍ő�l�ɒu�������� */
{
    int x, y, i, j; /* ���[�v�ϐ�       */
    int max;        /* �ߖT�̍ő�K���l */

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
/* Sobel�t�B���^�ɂ��1�������l�ɂ��� */
/* n=1�F�������C n=2�F�c�����C n=3�F�����̑傫�� */
{
    int x, y; /* ���[�v�ϐ� */
    double horizontal, vertical; /* �����������l�C�c���������l */
    int strength; /* �����̑傫�� */

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
/* ���v���V�A���ɂ��2�������l�ɂ��� */
{
    int x, y;     /* ���[�v�ϐ�   */
    int strength; /* �����̑傫�� */

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
/* ���v���V�A���ɂ��2�������l + MAX_BRIGHTNESS �ɂ��� */
{
    int x, y;     /* ���[�v�ϐ�   */
    int strength; /* �����̑傫�� */

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
/* ���ϊK�����Â���f �� 0 */
{
    int x, y;       /* ���[�v�ϐ� */
    double average; /* ���ϊK���l */

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
/* ���ϊK����薾�邢��f �� MAX_BRIGHTNESS */
{
    int x, y;       /* ���[�v�ϐ� */
    double average; /* ���ϊK���l */

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
/* ���]���� */
{
    int x, y; /* ���[�v�ϐ� */

    for ( y = 0; y < height; y ++ ){
        for ( x = 0; x < width; x ++ ){
            image2[y][x] = MAX_BRIGHTNESS - image1[y][x];
        }
    }
}

void _subtract(void)
/* ���摜�i=original[img_num-1][][] �j�Ƃ̍����̐�Βl�ɂ��� */
{
    int x, y; /* ���[�v�ϐ� */

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
/* ���ϐ��Ő錾����Ă��鉡width �cheight �̉摜 �f�[�^ */
/* image1, image2 ��p���ĉ摜�̃t�B���^�����O�����s����v���O���� */
/* ��� image1 ���� image2 �ւ̕ϊ����s�����̂Ƃ���D */
/* int number : �K�p����t�B���^�̔ԍ�(1,2,...,MAX_FILTER_NUMBER) */
/* int img_number : �ΏۂƂ��錴�摜�̔ԍ�(1,2,...,MAX_LEARN) */
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
