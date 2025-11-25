// gcc sobel.c -o sobel
// ./sobel garb34.pgm sobelmag.pgm sobelout1.pgm sobelout2.pgm 110 40

#include <stdio.h> /* Sobel.c */
#include <math.h>
#include <stdlib.h>

int pic[256][256];
int outpicx[256][256];
int outpicy[256][256];
int maskx[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
int masky[3][3] = {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};
double ival[256][256], maxival;

void main(argc, argv) int argc;
char **argv;
{
    int i, j, p, q, mr, sum1, sum2;
    double thresholdHI, thresholdLO;
    FILE *fo1, *fo2, *fo3, *fp1, *fopen();
    char *foobar;

    if (argc < 7)
    {
        printf("Usage: %s garb34.pgm sobelmag.pgm sobelout1.pgm sobelout2.pgm HI LO\n", argv[0]);
        return;
    }

    argc--;
    argv++;
    foobar = *argv;
    fp1 = fopen(foobar, "rb");
    if (fp1 == NULL)
    {
        printf("File not found: %s\n", foobar);
        return;
    }

    argc--;
    argv++;
    foobar = *argv;
    fo1 = fopen(foobar, "wb");
    fprintf(fo1, "P5\n256 256\n255\n");

    argc--;
    argv++;
    foobar = *argv;
    fo2 = fopen(foobar, "wb");
    fprintf(fo2, "P5\n256 256\n255\n");

    argc--;
    argv++;
    foobar = *argv;
    fo3 = fopen(foobar, "wb");
    fprintf(fo3, "P5\n256 256\n255\n");

    argc--;
    argv++;
    foobar = *argv;
    thresholdHI = atof(foobar);

    argc--;
    argv++;
    foobar = *argv;
    thresholdLO = atof(foobar);

    // read image
    char header[100];
    fgets(header, sizeof(header), fp1); // "P5\n"
    fgets(header, sizeof(header), fp1); // "256 256\n"
    fgets(header, sizeof(header), fp1); // "255\n"

    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            pic[i][j] = getc(fp1);
            pic[i][j] &= 0377;
        }
    }

    // scanning convolutions
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
                    sum1 += pic[i + p][j + q] * maskx[p + mr][q + mr];
                    sum2 += pic[i + p][j + q] * masky[p + mr][q + mr];
                }
            }
            outpicx[i][j] = sum1;
            outpicy[i][j] = sum2;
        }
    }

    // gradient magnitude + normalization
    maxival = 0;
    for (i = mr; i < 256 - mr; i++)
    {
        for (j = mr; j < 256 - mr; j++)
        {
            ival[i][j] = sqrt((double)((outpicx[i][j] * outpicx[i][j]) +
                                       (outpicy[i][j] * outpicy[i][j])));
            if (ival[i][j] > maxival)
                maxival = ival[i][j];
        }
    }

    // output magnitude image into sobelmag.pgm
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            ival[i][j] = (ival[i][j] / maxival) * 255;
            fprintf(fo1, "%c", (char)((int)(ival[i][j])));
        }
    }

    // output LO threshold image into sobelout1.pgm
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            if (ival[i][j] > thresholdLO)
                fprintf(fo2, "%c", (char)255);
            else
                fprintf(fo2, "%c", (char)0);
        }
    }

    // output HI threshold image into sobelout2.pgm
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            if (ival[i][j] > thresholdHI)
                fprintf(fo3, "%c", (char)255);
            else
                fprintf(fo3, "%c", (char)0);
        }
    }

}
