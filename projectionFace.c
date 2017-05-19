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
    if (result) {
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
	double cosine = cos(angle);
	double sinus = sin(angle);
	double matrixRotation[4][4] = {
		{	cosine, 	0,  sinus, 0},
		{	0	,				1, 	 0 ,    0 },
		{ -sinus ,   0,  cosine, 0},
		{	0	,				0, 	0 ,     1},
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
		double angle = noFrame * 0.3;
		Pixel** m = raster(HEIGTH,WIDTH);

		double **depthBuffer = (double **) malloc(WIDTH * sizeof(double *));
		for(i = 0;i < WIDTH;i++) {
			depthBuffer[i] = (double *) malloc(HEIGTH * sizeof(double));
		}

		for(i = 0;i < WIDTH;i++)
			for(j = 0;j < HEIGTH;j++)
				depthBuffer[i][j] = INFINITY;

		int translateZ = 400;
		vector3 camera;

		vector3 light;
		light.x = 0;
		light.y = 0;
		light.z = 1;



		for(j = 0;j < numberObjects;j++) {
				Line vector0 =	translateLine(rotateLine(angle,scale(triangles[j].lines[0], 100.0)),0,0,translateZ);
				Line vector1 =	translateLine(rotateLine(angle,scale(triangles[j].lines[1], 100.0)),0,0,translateZ);
				Line vector2 =	translateLine(rotateLine(angle,scale(triangles[j].lines[2], 100.0)),0,0,translateZ);

				vector3 aux0 = minus(vector1.a, vector0.a);
				vector3 aux1 = minus(vector2.a, vector0.a);

				vector3 normal = crossProduct(aux0, aux1);
				normal = normalize(normal);

				camera.x = vector0.a.x;
				camera.y = vector0.a.y;
				camera.z = -vector0.a.z;
				camera =  normalize(camera);

				double dotCamera = dotProduct(normal,camera);
				dotCamera =  max(0, dotCamera);
				if(dotCamera == 0) continue;

				 Line v0, v1, v2;
				 v0 = projection(vector0,HEIGTH,WIDTH);
				 v1 = projection(vector1,HEIGTH,WIDTH);
				 v2 = projection(vector2,HEIGTH,WIDTH);

				 double xmin = min3(v0.a.x, v1.a.x, v2.a.x);
				 double ymin = min3(v0.a.y, v1.a.y, v2.a.y);
				 double xmax = max3(v0.a.x, v1.a.x, v2.a.x);
				 double ymax = max3(v0.a.y, v1.a.y, v2.a.y);

				 if (xmin > HEIGTH - 1 || xmax < 0 || ymin > WIDTH - 1 || ymax < 0)
				 		continue;

				 uint32_t x0 = max((int32_t) 0, (int32_t)(floor(xmin)));
				 uint32_t x1 = min((int32_t)(WIDTH - 1), (int32_t)(floor(xmax)));
				 uint32_t y0 = max((int32_t) 0, (int32_t)(floor(ymin)));
				 uint32_t y1 = min((int32_t)(HEIGTH - 1), (int32_t)(floor(ymax)));
				 double area = edgeFunction(v0.a, v1.a, v2.a);

			 for (uint32_t y = y0; y <= y1; ++y) {
						for (uint32_t x = x0; x <= x1; ++x) {
							Point sample;
							sample.x = x + 0.5;
							sample.y = y + 0.5;

							double w0 = edgeFunction(v1.a, v2.a, sample);
							double w1 = edgeFunction(v2.a, v0.a, sample);
							double w2 = edgeFunction(v0.a, v1.a, sample);

							if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
	                    w0 /= area;
	                    w1 /= area;
	                    w2 /= area;
											float oneOverZ = v0.a.z * w0 + v1.a.z * w1 + v2.a.z * w2;
	                    float z = 1 / oneOverZ;

											if (z < depthBuffer[x][y]) {
												depthBuffer[x][y] = z;

												float px = (v0.a.x/ -v0.a.z) * w0 + (v1.a.x/ -v1.a.z) * w1 + (v2.a.x / -v2.a.z) * w2;
												float py = (v0.a.y/ -v0.a.z) * w0 + (v1.a.y/ -v1.a.z) * w1 + (v2.a.y / -v2.a.z) * w2;

												vector3 pt;
												pt.x = px * z;
												pt.y = py * z;
												pt.z =  -z;
												pt = normalize(pt);

												float dotLigth = dotProduct(normal, pt);
												float diffuseLight =  max(0.0, dotLigth);

												float specularLight = 0;
												if (diffuseLight <= 0) specularLight = 0;
												float ambientLight = 0.1;

												float totalLight = ambientLight + diffuseLight + specularLight ;
												//totalLight = min(1, totalLight);

												m[x][y].red = min((int)200 * totalLight, 255);
											 	m[x][y].green = min((int)200 * totalLight, 255);
											 	m[x][y].blue = min((int)200 * totalLight, 255);
									}
							}
					}
				}
		}
		char str[80];
		printf("Guardando imagen %d\n", noFrame);
		sprintf(str,"img/image%d.ppm", noFrame);
		FILE *fp = createFile(str,HEIGTH,WIDTH,MAX_NUMBER,MAGIC_NUMBER);
		saveMatrix(fp, m,HEIGTH,WIDTH);
		free(m);
		free(depthBuffer);
		fclose(fp);
}

int main() {
	int i,j;
	FILE *file;
	file = fopen("cube.raw", "r" );
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
