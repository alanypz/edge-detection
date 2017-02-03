/*
 *  Canny Algorithm
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define PICSIZE 256
#define MAXMASK 100
#define true 1
#define false 0

int    pic[PICSIZE][PICSIZE];
double outpicp[PICSIZE][PICSIZE];
double outpicq[PICSIZE][PICSIZE];
double maskp[MAXMASK][MAXMASK];
double maskq[MAXMASK][MAXMASK];
double ival[PICSIZE][PICSIZE];
double mag[PICSIZE][PICSIZE];
double peaks[PICSIZE][PICSIZE];
double thresh[PICSIZE][PICSIZE];
int histogram[PICSIZE];

void  readHeader(FILE *ifp);

int main(int argc, const char * argv[])
{
    int i, j, p, q, low, high, mr, centerp, centerq, sum, quit, moretodo;
    double  threshold, sigma, maskval, sump, sumq, mapival, slope;

    //  Static filenames
    FILE *ifp, *ofp1, *ofp2, *ofp3;
    ifp = fopen("garb34.pgm","rb");
    ofp1 = fopen("mag.pgm","wb");
    ofp2 = fopen("peaks.pgm","wb");
    ofp3 = fopen("thresh.pgm","wb");

    //  Fix .pgm "shift" to the right.
    readHeader(ifp);

    threshold = .02;
    sigma = 1.0;

    printf("\nSetting sigma value to default: %lf\n", sigma);
    // printf("Enter threshold percentage, as decimal:\n");
    // scanf("%lf", &threshold); //  Values for best outputs < .14
    printf("Setting threshold value to: %lf\n", threshold);

    mr = (int)(sigma * 3);
    centerp = (MAXMASK / 2);
    centerq = (MAXMASK / 2);

    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            pic[i][j] = getc(ifp);
            pic[i][j] &= 0377;
        }
    }

    fclose(ifp);

    //  Canny Algorithm – Part 1: Compute Gradient Magnitude
    for (p = -mr; p <= mr; p++)
    {
        for (q = -mr; q <= mr; q++)
        {
            maskval = exp(-((p*p+q*q) / (2*sigma*sigma)));
            maskp[q+centerq][p+centerp] = -p * maskval;
            maskq[q+centerq][p+centerp] = -q * maskval;
        }
    }

    for (i = mr; i <= 256 - mr; i++)
    {
        for (j = mr; j <= 256 - mr; j++)
        {
            sump = 0;
            sumq = 0;

            for (p = -mr; p <= mr; p++)
            {
                for (q = -mr; q<= mr; q++)
                {
                    sump += pic[i+p][j+q] * maskp[p+centerq][q+centerp];
                    sumq += pic[i+p][j+q] * maskq[p+centerq][q+centerp];
                }
            }
            outpicp[i][j] = sump;
            outpicq[i][j] = sumq;
        }
    }

    mapival = 0;
    for (i = mr; i < 256 - mr; i++)
    {
        for (j = mr; j < 256 - mr; j++)
        {
            ival[i][j] = sqrt((double)((outpicp[i][j] * outpicp[i][j]) + (outpicq[i][j] * outpicq[i][j])));

            if (ival[i][j] > mapival)
            {
                mapival = ival[i][j];
            }
        }
    }

    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            mag[i][j] = (ival[i][j] / mapival) * 255;
        }
    }

    //  Canny Algorithm – Part 2: Compute Peaks
    for(i = mr; i < 256 - mr; i++)
    {
       for(j = mr; j < 256 - mr; j++)
       {
            peaks[i][j] = 0;

            if(outpicp[i][j] == 0.0)
            {
                outpicp[i][j] = 0.0001;
            }

            slope = outpicq[i][j]/outpicp[i][j];

            if( (slope <= 0.4142)&&(slope > -0.4142))
            {

                if((mag[i][j] > mag[i][j-1])&& (mag[i][j] > mag[i][j+1]))
                {
                    peaks[i][j] = 255;
                }
            }
            else if( (slope <= 2.4142)&&(slope > 0.4142))
            {
                if((mag[i][j] > mag[i+1][j+1])&&(mag[i][j] > mag[i-1][j-1]))
                {
                    peaks[i][j] = 255;
                }
            }
            else if( (slope <= -0.4142)&&(slope > -2.4142))
            {
                if((mag[i][j] > mag[i-1][j+1])&&(mag[i][j] > mag[i+1][j-1]))
                {
                    peaks[i][j] = 255;
                }
           }
           else
           {
                if((mag[i][j] > mag[i-1][j])&&(mag[i][j] > mag[i+1][j]))
                {
                    peaks[i][j] = 255;
                }
            }
        }
    }

    fprintf(ofp1, "P5\n256 256\n255\n");
    fprintf(ofp2, "P5\n256 256\n255\n");
    fprintf(ofp3, "P5\n256 256\n255\n");

    //  Output magnitude (mag.pgm) and peaks (peaks.pgm).
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            fprintf(ofp1,"%c",(char)((int)(mag[i][j])));
            fprintf(ofp2, "%c", (char)((int)(peaks[i][j])));
        }
    }

    //  Canny Aglorithm – Part 3: Automatically Compute High and Low
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            histogram[(int)mag[i][j]%256]++;
        }
    }

    quit = threshold * PICSIZE * PICSIZE;
    sum = 0;
    high = 255;

    for (i = 256; i > 0; i--)
    {
        sum += histogram[i];

        if (sum > quit)
        {
            high = i;
            break;
        }
    }

    low = high * .35;

    // printf("quit: %d\n", quit); // bug testing
   printf("\nHigh threshold: %d\nLow threshold: %d\n", high, low);

    // Canny Algorithm – Part 4: Double Thresholding
    for(i = 0; i < 256; i++)
    {
        for(j = 0; j < 256; j++)
        {
            if (peaks[i][j] == 255)
            {
                if (ival[i][j] > high)
                {
                    peaks[i][j] = 0;
                    thresh[i][j] = 255;
                }
                else if (ival[i][j] < low)
                {
                    peaks[i][j] = 0;
                    thresh[i][j] = 0;
                }
            }
        }
    }

    do
    {
        moretodo = false;

        for(i = 0; i < 256; i++)
        {
            for(j = 0 ; j < 256; j++)
            {
                if (peaks[i][j] == 255)
                {
                    for (p = -1;p <= 1; p++)
                    {
                        for (q = -1; q <= 1; q++)
                        {
                            if (thresh[i+p][j+q] == 255)
                            {
                                peaks[i][j] = 0;
                                thresh[i][j] = 255;
                                moretodo = true;
                             }
                        }
                    }
                }
            }
        }
    }
    while (moretodo);

    //  Output result of double-thresholding (thresh.pgm)
    for (i = 0; i < 256; i++)
    {
        for (j = 0;j < 256; j++)
        {
            fprintf(ofp3, "%c", (char)thresh[i][j]);
        }
    }

    fclose(ofp1);
    fclose(ofp2);
    fclose(ofp3);

    return 0;
}

//  Fix for .pgm "shifting" right.
void readHeader(FILE * ifp)
{
    int headerLines = 3;
    char line[PICSIZE];

    while(headerLines > 0)
    {
        fgets(line, PICSIZE, ifp);

        //  Ignore empty lines and comments.
        if((strlen(line) <= 1) || line[0] == '#')
        {
            continue;
        }

        headerLines -= 1;
    }
}
