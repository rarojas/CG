#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include "projection.h"

#define  MAX_NUMBER 255
#define  WIDTH 720
#define  HEIGTH 480
#define MAGIC_NUMBER "P3"
#define FILE_NAME "projection.ppm"

typedef struct {
	int red;
	int green;
	int blue;
} Pixel;



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


typedef struct{
	Point a;
	Point b;
  Point c;
  Line* lines;
	Color color;
} Triangle;



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

FILE* createFile(char*  filename) {
	FILE *fp;
	fp = fopen(filename, "w+" );
	fprintf(fp, "%s\n# %s\n%i %i\n%i\n",MAGIC_NUMBER, filename, WIDTH,HEIGTH, MAX_NUMBER);
	return fp;
}

void multiply (double a[4][4],double b[4], double *resultmatrix)
{   int i,j;
    for(i = 0; i < 4; i++) {
      double result  = 0;
      for(j = 0; j < 4; j++)
      	result += a[i][j] * b[j];
      resultmatrix[i] = result;
    }
}

void drawLine(Pixel **matrix, Line line){
		int dx =  line.b.x - line.a.x;
		int dy =  line.b.y - line.a.y;
		int dxAbs = abs(dx);
		int dyAbs = abs(dy);
		int steps = dxAbs > dyAbs ? dxAbs: dyAbs;
		double xIncrement =  dx/(double)steps;
		double yIncrement =  dy/(double)steps;
		double x = (double)line.a.x;
		double y = (double)line.a.y;
		int v;
		for(v = 0;v < steps; v++){
				x += xIncrement;
				y += yIncrement;
				int xCord = round(x);
				int yCord = round(y);
				if(yCord >= HEIGTH && yCord <= 0 && xCord >= WIDTH && xCord <= 0)
			 		return;
				matrix[xCord][yCord].red = line.color.red;
				matrix[xCord][yCord].green = line.color.green;
				matrix[xCord][yCord].blue = line.color.blue;
		}
}


char** str_split(char* a_str, const char a_delim){
    char** result    = 0;
    size_t count     = 0;
    char* tmp        = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;
    while (*tmp){
        if (a_delim == *tmp) {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }
    count += last_comma < (a_str + strlen(a_str) - 1);
    count++;
    result = malloc(sizeof(char*) * count);
    if (result)
    {
        size_t idx  = 0;
        char* token = strtok(a_str, delim);
        while (token){
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }
    return result;
}



Line rotateLine(double angle, Line line){
	double matrixRotation[4][4] = {
		{	cos(angle), 		0,  	-sin(angle), 0},
		{	0	,							1, 		0 , 0 },
		{ sin(angle) , 	0,    cos(angle), 0},
		{	0	,							0, 		0 , 1},
	 	};

		double aMatrix[4] = { line.a.x , line.a.y, line.a.z, 1 };
		multiply(matrixRotation,aMatrix, (double *)&aMatrix);
		line.a.x = aMatrix[0];
		line.a.y = aMatrix[1];
		line.a.z = aMatrix[2];
		double bMatrix[4] = { line.b.x , line.b.y, line.b.z , 1};
		multiply(matrixRotation,bMatrix, (double *)&bMatrix);
		line.b.x = bMatrix[0];
		line.b.y = bMatrix[1];
		line.b.z = bMatrix[2];
		return line;
}


Line translateLine(Line line, double x,double y,double z){
	double scale = 100;
		double matrixTranlate[4][4] = {
			{	scale, 		0,  	0,   x },
			{	0	,		scale, 		0,   y},
			{ 0 , 	0,    scale ,  z},
			{ 0 , 	0,    0    ,   	 1}
		 };
		double aMatrix[4] = { line.a.x , line.a.y, line.a.z ,1};
		multiply(matrixTranlate,aMatrix, (double *)&aMatrix);
		line.a.x = aMatrix[0];
		line.a.y = aMatrix[1];
		line.a.z = aMatrix[2];
		double bMatrix[4] = { line.b.x , line.b.y, line.b.z, 1};
		multiply(matrixTranlate,bMatrix, (double *)&bMatrix);
		line.b.x = bMatrix[0];
		line.b.y = bMatrix[1];
		line.b.z = bMatrix[2];
		return line;
}

Line projection(Line line){
  int f = 300;
  double matrixProjection[4][4] = {
			{	f, 		0,   300 , 	0},
			{	0	,		f ,  300 ,  0},
		 	{ 0 , 	0,   f ,	0 },
			{ 0 , 	0,   0 , 	f },
	 };
  double aMatrix[4] = { line.a.x / line.a.z, line.a.y / line.a.z, 1 , 1 / line.a.z };
  double bMatrix[4] = { line.b.x / line.b.z, line.b.y / line.b.z, 1, 1 / line.b.z};

  double* resultA = malloc(sizeof(double) * 4);
  multiply(matrixProjection,aMatrix, resultA);

  double* resultB = malloc(sizeof(double) * 4);
  multiply(matrixProjection,bMatrix, resultB);

  Line projection;
  projection.a.x = resultA[0];
  projection.a.y = resultA[1];

  projection.b.x = resultB[0];
  projection.b.y = resultB[1];
  projection.color.red = 0;
  projection.color.green = 0;
  projection.color.blue = 0;
	//printf("%f,%f %f,%f\n",projection.a.x,projection.a.y,projection.b.x,projection.b.y);

  free(resultA);
  free(resultB);
  return projection;
}

int getSize(char** ch) {
      int length = 0;
      while (*(ch + length)) {
        length++;
      }
      return length;
}

void createFrame(Triangle* triangles,int noFrame,int numberObjects,int length) {
		int  i ,j,index ;
		double angle = noFrame * 0.1;
		Pixel **m = (Pixel **) malloc(WIDTH * sizeof(Pixel *));
		for(i = 0;i < WIDTH;i++) {
			m[i] = (Pixel *) malloc(HEIGTH * sizeof(Pixel));
		}
		printf("angle %f cos %f \n",angle,cos(angle));

		raster(m);

		int translateZ = 300;
		vector3 camera;
		camera.x = 0;
		camera.y = 0;
		camera.z = 1;


		for(j = 0;j < numberObjects;j++) {
			Line vector0 =	 translateLine(rotateLine(angle,triangles[j].lines[0]),0,0,translateZ);
			Line vector1 =	 translateLine(rotateLine(angle,triangles[j].lines[1]),0,0,translateZ);
			Line vector2 =	 translateLine(rotateLine(angle,triangles[j].lines[2]),0,0,translateZ);

			vector3 aux0 = minus(vector1.a, vector0.a);
			vector3 aux1 = minus(vector2.a, vector0.a);

			vector3 normal = crossProduct(aux0, aux1);
			normal = normalize(normal);
			double dotCamera = dotProduct(normal,camera);
			dotCamera =  max(0.00001, dotCamera);
			if(dotCamera == 0.00001)
				continue;


			for(i = 0;  i < length - 1; i++){
				Line line = rotateLine(angle,triangles[j].lines[i]);
				line = translateLine(line,0,0,translateZ);
				line = projection(line);
				vector3 aux0 = minus(vector1.a, vector0.a);
				vector3 aux1 = minus(vector2.a, vector0.a);

				vector3 normal = crossProduct(aux0, aux1);
				normal = normalize(normal);
				double dotCamera = dotProduct(normal,camera);
				dotCamera =  max(0.00001, dotCamera);
				if(dotCamera == 0.00001)
					continue;


				drawLine(m,line);
			}
		}
		char str[80];
		sprintf(str,"img/image%d.ppm", noFrame);
		FILE *fp = createFile(str);
		saveMatrix(fp, m);
		free(m);
		fclose(fp);
}

int main() {
	int i,j;
	FILE *file;
	file = fopen("cow.raw", "r" );
  if(file == NULL){
    printf("Hubo un error al abrir el archivo\n" );
    return -1;
  }
  Triangle* triangle = malloc(sizeof(Triangle) * 0);
	int length;
  char *line = NULL;
  size_t len = 0;
  char read;
  char** tokens;
	int numberObjects = 0, index = 0;
  while ((read = getline(&line, &len, file)) != -1) {
      tokens = str_split(line,' ');
      length = getSize(tokens) / 3;
			triangle = (Triangle *)realloc(triangle ,sizeof(Triangle) * ++numberObjects);
      triangle[index].lines = malloc(sizeof(Line) * length);
      if (tokens){
        for (i = 0, j = 0;*(tokens + i); j++) {
            triangle[index].lines[j].a.x = atof(*(tokens + i));
						i++;
            triangle[index].lines[j].a.y = atof(*(tokens + i));
						i++;
            triangle[index].lines[j].a.z = atof(*(tokens + i));
						i++;
      	}
				for(i = 0,j = 1;i < length; i++,j++){
	        	if(j == length)
							j = 0;
	          triangle[index].lines[i].b.x = triangle[index].lines[j].a.x;
	          triangle[index].lines[i].b.y = triangle[index].lines[j].a.y;
	          triangle[index].lines[i].b.z = triangle[index].lines[j].a.z;
	      }
    	}
			index++;
		}

		int totalFrames = 10;

		int noFrame = 0;
		for(noFrame = 0; noFrame < totalFrames; noFrame++){
				createFrame(triangle,noFrame,numberObjects,length);
		}


  return 0;
}
