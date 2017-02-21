#ifndef LINE_H_
#define LINE_H_

typedef struct {
	int red;
	int green;
	int blue;
} Pixel;

typedef struct {
   int  x;
   int 	y;
} Point;

typedef struct {
	int red;
	int green;
	int blue;
} Color;

typedef struct{
	Point a;
	Point b;
	Color color;
} Line;

void generateRandomColor(Color *color);
void drawLine(Pixel **matrix, Line line);
void generateRandomPoint(Point *point);
FILE * createFile();
void saveMatrix(FILE *fp, Pixel** matrix);
void raster(Pixel** matrix);
void generateRandomLine(Line* line);

#endif
