#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "projection.h"
#define  MAX_NUMBER 255
#define MAGIC_NUMBER "P3"
#define FILE_NAME "projection.ppm"
#define  WIDTH 1920
#define  HEIGTH 1080


double getH1(double s) { return (2.0 * pow(s,3)) - (3.0 * s * s) + 1; }

double getH2(double s) { return (-2.0 * pow(s,3)) + (3.0 * s * s); }

double getH3(double s) { return pow(s,3) - (2.0 * s * s) + s; }

double getH4(double s) { return pow(s,3) - (s * s) ; }

void hermite(int x[4], int y[4], Pixel** matrix){
  double s;
  Line line;
  line.a.x = x[0];
  line.a.y = y[0];
  line.b.x = x[0] + x[2];
  line.b.y = y[0] + y[2];
  drawLine(matrix, line, HEIGTH, WIDTH);

  line.a.x = x[1];
  line.a.y = y[1];
  line.b.x = x[1] + x[3];
  line.b.y = y[1] + y[3];
  drawLine(matrix, line, HEIGTH, WIDTH);
  for (s = 0.0; s < 1.0; s += 0.0001) {
    double h1 = getH1(s);
    double h2 = getH2(s);
    double h3 = getH3(s);
    double h4 = getH4(s);
    double xt = (h1 * x[0]) + (h2 * x[1]) + (h3 * x[2]) + (h4 * x[3]);
    double yt = (h1 * y[0]) + (h2 * y[1]) + (h3 * y[2]) + (h4 * y[3]);
    putpixel(xt, yt, matrix, 0,HEIGTH,WIDTH);
 }
}

int main() {
    int x[4], y[4];
    int i;
    printf ("Enter the x- and y-coordinates \n");
    for (i = 0; i < 2; i++)
	     scanf ("%d,%d", &x[i], &y[i]);
    printf("Enter the tangents at p1,p4\n");
    for (i = 2; i < 4; i++)
       scanf ("%d,%d", &x[i], &y[i]);
    Pixel** matrix = raster(HEIGTH,WIDTH);
    hermite(x, y, matrix);
    char* str = "img/hermite.ppm";
    FILE *fp = createFile(str,HEIGTH,WIDTH,MAX_NUMBER,MAGIC_NUMBER);
    saveMatrix(fp, matrix,HEIGTH,WIDTH);
    free(matrix);
}
