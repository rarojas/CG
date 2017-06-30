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
		Pixel **m = raster(WIDTH,HEIGTH);
		printf("angle %f cos %f \n",angle,cos(angle));


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
				line = projection(line,WIDTH,HEIGTH);
				vector3 aux0 = minus(vector1.a, vector0.a);
				vector3 aux1 = minus(vector2.a, vector0.a);

				vector3 normal = crossProduct(aux0, aux1);
				normal = normalize(normal);
				double dotCamera = dotProduct(normal,camera);
				dotCamera =  max(0.00001, dotCamera);
				if(dotCamera == 0.00001)
					continue;

				drawLine(m,line,HEIGTH,WIDTH);
			}
		}
		char str[80];
		sprintf(str,"img/image%d.ppm", noFrame);
		FILE *fp = createFile(str,HEIGTH,WIDTH,MAX_NUMBER,MAGIC_NUMBER);
		saveMatrix(fp, m,HEIGTH,WIDTH);
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
