// Jake Weber

// Expected command line format:
//     gcc marrh.c -o marrh
//     ./marrh.exe garb34.pgm magnitude.pgm peaks.pgm Final Edges.pgm 1.0 5.0
// Input: input.pgm
// Output: magnitude.pgm (contains the image of magnitudes)
// Output: peaks.pgm (contains the image of peaks determineed from the non-maxima suppression algorithm)
// Output: Final Edges.pgm (contains the fully computed edge detected image)
// Parameter: sigma (used as the unit of standard deviation in pixels)
// Parameter: percent (used as the percent of pixels marked as HI automatically)

#include <stdio.h>     /* Marr-Hildreth + Canny-like NMS + Hysteresis */
#include <math.h>
#include <stdlib.h>

#define PICSIZE 256
#define MAXMASK 100

int pic[PICSIZE][PICSIZE];
double xconv[PICSIZE][PICSIZE] = {0};
double yconv[PICSIZE][PICSIZE] = {0};
double ival[PICSIZE][PICSIZE];
int cand[PICSIZE][PICSIZE];   // candidate peaks from NMS
int peaks[PICSIZE][PICSIZE];  // hysteresis peaks
int final[PICSIZE][PICSIZE];  // final edge map
int histogram[256] = {0};

double xmask[MAXMASK][MAXMASK];
double ymask[MAXMASK][MAXMASK];

int main(int argc, char **argv)
{
    int i, j, p, q, mr, centx, centy, HI, LO, areaOfTops = 0;
    double sig, sumX, sumY, maxival, slope, percent;
    FILE *fp1, *fo1, *fo2, *fo3;
    char *foobar;

    if (argc < 7) {
        printf("Usage: %s input.pgm magnitude.pgm peaks.pgm Final Edges.pgm sigma percent\n", argv[0]);
        return 1;
    }

    foobar = argv[1];
    fp1 = fopen(foobar, "rb");
    if (fp1 == NULL) { printf("File not found: %s\n", foobar); return 1; }

    foobar = argv[2];
    fo1 = fopen(foobar, "wb");
    fprintf(fo1, "P5\n256 256\n255\n");

    foobar = argv[3];
    fo2 = fopen(foobar, "wb");
    fprintf(fo2, "P5\n256 256\n255\n");

    foobar = argv[4];
    fo3 = fopen(foobar, "wb");
    fprintf(fo3, "P5\n256 256\n255\n");

    sig = atof(argv[5]);
    mr = (int)(sig * 3);

    percent = atof(argv[6]);

    // read image
    char header[100];
    fgets(header, sizeof(header), fp1); // "P5\n"
    fgets(header, sizeof(header), fp1); // "256 256\n"
    fgets(header, sizeof(header), fp1); // "255\n"

    for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
                {
                  pic[i][j]  =  getc (fp1);
                  pic[i][j]  &= 0377;
                }
        }

    // build derivative of Gaussian masks
    for (p = -mr; p <= mr; p++) {
        for (q = -mr; q <= mr; q++) {
            xmask[p + mr][q + mr] = p * exp(-(p*p + q*q) / (2 * sig * sig));
            ymask[p + mr][q + mr] = q * exp(-(p*p + q*q) / (2 * sig * sig));
        }
    }

    // convolution
    for (i = mr; i < 256 - mr; i++) {
        for (j = mr; j < 256 - mr; j++) {
            sumX = 0;
            sumY = 0;
            for (p = -mr; p <= mr; p++) {
                for (q = -mr; q <= mr; q++) {
                    sumX += pic[i + p][j + q] * xmask[p + mr][q + mr];
                    sumY += pic[i + p][j + q] * ymask[p + mr][q + mr];
                }
            }
            xconv[i][j] = sumX;
            yconv[i][j] = sumY;
        }
    }
   
    // gradient magnitude + normalization
    maxival = 0;
    for (i = mr; i < 256 - mr; i++) {
        for (j = mr; j < 256 - mr; j++) {
            ival[i][j] = sqrt(xconv[i][j]*xconv[i][j] + yconv[i][j]*yconv[i][j]);
            if (ival[i][j] > maxival) maxival = ival[i][j];
        }
    }

    // ival only contains gradient values from mr to (255 - mr)
    for (i = mr; i < 256 - mr; i++) {
        for (j = mr; j < 256 - mr; j++) {
            ival[i][j] = (ival[i][j] / maxival) * 255;
        }
    }

    // write magnitude map
    for (i = 0; i < 256; i++) {
        for (j = 0; j < 256; j++) {
            fprintf(fo1, "%c", (char)((int)(ival[i][j])));
        }
    }

    // non-maximum suppression
    for (i = mr; i < 256 - mr; i++) {
        for (j = mr; j < 256 - mr; j++) {
            double mag = ival[i][j];
            if (xconv[i][j] == 0.0) xconv[i][j] = 0.00001;

            slope = yconv[i][j] / xconv[i][j];

            if ((slope <= 0.4142) && (slope > -0.4142)) {
                if ((mag > ival[i-1][j]) && (mag > ival[i+1][j]))
                    cand[i][j] = 255;
            }
            else if ((slope <= 2.4142) && (slope > 0.4142)) {
                if ((mag > ival[i+1][j+1]) && (mag > ival[i-1][j-1]))
                    cand[i][j] = 255;
            }
            else if ((slope <= -0.4142) && (slope > -2.4142)) {
                if ((mag > ival[i-1][j+1]) && (mag > ival[i+1][j-1]))
                    cand[i][j] = 255;
            }
            else if ((slope >= 2.4142) || (slope <= -2.4142)) {
                if ((mag > ival[i][j-1]) && (mag > ival[i][j+1]))
                    cand[i][j] = 255;
            }
        }
    }

    // write peaks map
    for (i = 0; i < 256; i++) {
        for (j = 0; j < 256; j++) {
            fprintf(fo2, "%c", (char)((int)(cand[i][j])));
        }
    }

    // Compute Histogram of scaled magnitudes
    for (i = 0; i < PICSIZE; i++) {
        for (j = 0; j < PICSIZE; j++) {
            histogram[(int)ival[i][j]]++;
        }
    }

    // Automatically Get HI
    int cutoff = (int)(percent * 0.01 * PICSIZE * PICSIZE);
    for (i = 255; i >= 1; i--) {
        areaOfTops += histogram[i];
        if (areaOfTops > cutoff) {
            break;
        }
    }
    HI = i;
    LO = 0.35 * HI;

    // Hysteresis thresholding
    for (i = 0; i < 256; i++) {
        for (j = 0; j < 256; j++) {
            if (cand[i][j] == 255) {
                if (ival[i][j] >= HI) {
                    peaks[i][j] = 0;
                    final[i][j] = 255;
                }
                else if (ival[i][j] < LO) {
                    peaks[i][j] = 0;
                    final[i][j] = 0;
                }
                else {
                    peaks[i][j] = 255;
                    final[i][j] = 0;
                }
            }
        }
    }

    int moretodo = 1;
    while (moretodo == 1) {
        moretodo = 0;
        for (i = 1; i < 255; i++) {
            for (j = 1; j < 255; j++) {
                if (peaks[i][j] == 255) {
                    for (p = -1; p <= 1; p++) {
                        for (q = -1; q <= 1; q++) {
                            if (final[i+p][j+q] == 255) {
                                peaks[i][j] = 0;
                                final[i][j] = 255;
                                moretodo = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    // write Final Edges map
    for (i = 0; i < 256; i++) {
        for (j = 0; j < 256; j++) {
            fprintf(fo3, "%c", (char)((int)(final[i][j])));
        }
    }

    fclose(fp1);
    fclose(fo1);
    fclose(fo2);
    fclose(fo3);
    return 0;
}
