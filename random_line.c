#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lines.h"
#define MAX_NUMBER 255
#define NUMBER_LINES 1000
#define WIDTH 1920
#define HEIGTH 1080
#define MAGIC_NUMBER "P3"
#define FILE_NAME "lines.ppm"


int main() {
	int i;
	Line * lines = (Line * ) malloc(NUMBER_LINES * sizeof(Line));
	Pixel **m = (Pixel **) malloc(WIDTH * sizeof(Pixel *));
	for(i = 0;i < WIDTH;i++){
		m[i] = (Pixel *) malloc(HEIGTH * sizeof(Pixel));
	}

	for(i = 0; i < NUMBER_LINES; i ++) {
	 	generateRandomLine(&lines[i]);
	}
	raster(m);
	for(i = 0; i < NUMBER_LINES;i++)
		drawLine(m,lines[i]);
	FILE *fp = createFile();
	saveMatrix(fp, m);
	fclose(fp);
	return(0);
}

void generateRandomLine(Line* line){
	generateRandomColor(&(line->color));
	generateRandomPoint(&(line->a));
	generateRandomPoint(&(line->b));
}

void raster(Pixel **matrix){
	int i, j;
	for(i = 0; i < WIDTH;  i++) {
		for(j = 0; j < HEIGTH; j++){
			matrix[i][j].red = 255;
			matrix[i][j].green = 255;
			matrix[i][j].blue = 255;
		}
	}
}

void saveMatrix(FILE *fp, Pixel** matrix){
	int row, col;
	for(col = 0; col < HEIGTH; col++){
		for(row = 0; row < WIDTH;  row++) {
			fprintf(fp,"%i %i %i ",matrix[row][col].red, matrix[row][col].green, matrix[row][col].blue);
	 	}
	 	fputs("\n", fp);
	}
}

FILE* createFile() {
	FILE *fp;
	fp = fopen(FILE_NAME, "w+" );
	fprintf(fp, "%s\n# %s\n%i %i\n%i\n",MAGIC_NUMBER, FILE_NAME, WIDTH,HEIGTH, MAX_NUMBER);
	return fp;
}

void drawLine(Pixel **matrix, Line line){
		int dx =  line.b.x - line.a.x;
		int dy =  line.b.y - line.a.y;
		int dxAbs = abs(dx);
		int dyAbs = abs(dy);
		int steps = dxAbs > dyAbs ? dxAbs: dyAbs;
		float xIncrement =  dx/(float)steps;
		float yIncrement =  dy/(float)steps;
		float x = (float)line.a.x;
		float y = (float)line.a.y;
		int v ;
		for(v = 0;v < steps; v++){
			 x += xIncrement;
			 y += yIncrement;

			 int xCord = round(x);
			 int yCord = round(y);

			 matrix[xCord][yCord].red = line.color.red;
			 matrix[xCord][yCord].green = line.color.green;
			 matrix[xCord][yCord].blue = line.color.blue;
		}
}

void generateRandomColor(Color *color){
	color->red = rand() % MAX_NUMBER;
	color->green = rand() % MAX_NUMBER;
	color->blue = rand() % MAX_NUMBER;
}

void generateRandomPoint(Point *point){
	point->x = rand() % (WIDTH-1);
	point->y = rand() % (HEIGTH-1);
}
