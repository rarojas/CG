#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "projection.h"
#define  MAX_NUMBER 255
#define  WIDTH 1920
#define  HEIGTH 1080
#define MAGIC_NUMBER "P3"
#define FILE_NAME "projection.ppm"



void bezier (int x[4], int y[4],Pixel** matrix) {
    int i;
    double t;

    for (t = 0.0; t < 1.0; t += 0.0005) {
	     double xt = pow (1-t, 3) * x[0] + 3 * t * pow (1-t, 2) * x[1] + 3 * pow (t, 2) * (1-t) * x[2] + pow (t, 3) * x[3];
       double yt = pow (1-t, 3) * y[0] + 3 * t * pow (1-t, 2) * y[1] + 3 * pow (t, 2) * (1-t) * y[2] + pow (t, 3) * y[3];
	     putpixel(xt, yt,matrix, 0,HEIGTH,WIDTH);
    }
    int j;
    for (i = 0; i < 3; i++){
      Line line;
      line.a.x = x[i];
      line.a.y = y[i];
      line.b.x = x[i + 1];
      line.b.y = y[i + 1];
      drawLine(matrix,line,HEIGTH,WIDTH);
    }
}

int main() {
    int x[4], y[4];
    int i;
    printf ("Enter the x- and y-coordinates of the four control points.\n");
    for (i = 0; i < 4; i++)
	     scanf ("%d,%d", &x[i], &y[i]);
    Pixel** matrix = raster(HEIGTH,WIDTH);
    bezier (x, y,matrix);

    char* str = "img/bezier.ppm";
    FILE *fp = createFile(str,HEIGTH,WIDTH,MAX_NUMBER,MAGIC_NUMBER);
    saveMatrix(fp, matrix,HEIGTH,WIDTH);
    free(matrix);
}
