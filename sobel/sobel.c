/*
 *  Sobel Algorithm
 */

#include <stdio.h>
#include <math.h>
#include <string.h>

#define THRESHOLD_HIGH 40
#define THRESHOLD_LOW 20

int pic[256][256];
int outpicx[256][256];
int outpicy[256][256];
int maskx[3][3] = {{-1,0,1}, {-2,0,2}, {-1,0,1}};
int masky[3][3] = {{1,2,1}, {0,0,0}, {-1,-2,-1}};
double ival[256][256], maxival;

void  readHeader(FILE *ifp);

int main(int argc, const char * argv[])
{
    int i, j, p, q, mr, sum1, sum2;
    FILE *ifp, *opf1, *opf2, *opf3;

    //  Static filenames
    ifp = fopen("face05.pgm","rb");
    opf1 = fopen("mag.pgm","wb");
    opf2 = fopen("low.pgm","wb");
    opf3 = fopen("high.pgm","wb");

    readHeader(ifp);

    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
            pic[i][j] = getc (ifp);
            pic[i][j] &= 0377;
        }
    }


    fclose(ifp);

    mr = 1;
    for (i = mr; i < 256 - mr; i++)
    {
        for (j = mr; j < 256 - mr; j++)
        {
            sum1 = 0;
            sum2 = 0;

            for (p = -mr; p <= mr; p++)
            {
                for (q = -mr; q <= mr; q++)
                {
                    sum1 += pic[i+p][j+q] * maskx[p+mr][q+mr];
                    sum2 += pic[i+p][j+q] * masky[p+mr][q+mr];
                }
            }

            outpicx[i][j] = sum1;
            outpicy[i][j] = sum2;
        }
    }

    maxival = 0;
    for (i = mr; i < 256 - mr; i++)
    {
        for (j = mr; j< 256 - mr; j++)
        {
            ival[i][j] = sqrt((double)((outpicx[i][j] * outpicx[i][j]) + (outpicy[i][j] * outpicy[i][j])));

            if (ival[i][j] > maxival)
            {
                maxival = ival[i][j];
            }
        }
    }

    // printf("High threshold: %d\nLow threshold: %d\n", THRESHOLD_HIGH, THRESHOLD_LOW); // bug testing

    fprintf(opf1, "P5\n256 256\n255\n");
    fprintf(opf2, "P5\n256 256\n255\n");
    fprintf(opf3, "P5\n256 256\n255\n");

    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            int mag = (ival[i][j] / maxival) * 255;
            int low = mag > THRESHOLD_LOW ? 255 : 0;
            int high = mag > THRESHOLD_HIGH ? 255 : 0;

            fprintf(opf1,"%c",(char)mag);
            fprintf(opf2,"%c", (char)low);
            fprintf(opf3,"%c", (char)high);
        }
    }

    fclose(opf1);
    fclose(opf2);
    fclose(opf3);
}

//  Fix for .pgm "shifting" right.
void readHeader(FILE * ifp)
{
    int headerLines = 3;
    char line[256];

    while(headerLines > 0)
    {
        fgets(line, 256, ifp);

        //  Ignore empty lines and comments.
        if((strlen(line) <= 1) || line[0] == '#')
        {
            continue;
        }

        headerLines -= 1;
    }
}
