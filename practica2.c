#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include "projection.h"


#define  MAX_NUMBER 255
#define  WIDTH 480
#define  HEIGTH 480
#define MAGIC_NUMBER "P3"
#define FILE_NAME "projection.ppm"
#define FRAMES 36


const float nearClippingPLane = 1;
const float farClippingPLane = 100;
float focalLength = 20;
float filmApertureWidth = 0.980;
float filmApertureHeight = 0.735;

vector3 materialColor = { 43, 143, 219 };
vector3 ld =  { 211, 10, 10  };
vector3 background = { 100, 100, 100 };
vector3 lookAt = { 0 ,  0,  5};
vector3 camera = { 0,   0,  0};


float materialShininess = 30;

float ks = 0.1;
float kd = 1;
float ka = 0.01;
float ambientLight = 0.1;




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
		{	cosine, 	0,   -sinus,  0},
		{	0	,				1, 	 0 ,    0 },
		{ sinus ,  0,  cosine, 0},
		{	0	,				0, 	  0 ,     1},
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

void createFrame(Triangle* triangles,int noFrame,int numberObjects,int length,
  float l, float r, float t, float b, float near) {
    int  i ,j,index;
    Pixel** m = rasterBackground(HEIGTH,WIDTH,background);
    double **depthBuffer = createDeepBuffer(WIDTH, HEIGTH, INFINITY);
    double angle = 2 * M_PI * (float)noFrame / (float)FRAMES;
    vector3 light =  { cos(angle) * 3 ,   sin(angle) * 3,     sin(angle) * 2.5  - 2.5   };
    light =  (vector3){ 2, -2, -2};




		for(j = 0;j < numberObjects;j++) {
        Line vector0 =	rotateLineY(triangles[j].lines[0], angle);
        Line vector1 =	rotateLineY(triangles[j].lines[1], angle);
        Line vector2 =	rotateLineY(triangles[j].lines[2], angle);

        vector0 =	getPointToCamera(camera, lookAt, vector0);
        vector1 =	getPointToCamera(camera, lookAt, vector1);
        vector2 =	getPointToCamera(camera, lookAt, vector2);


				vector3 normal = getNormal(vector0.a,vector1.a,vector2.a);

        vector3 viewPointSurface = { vector0.a.x, vector0.a.y, vector0.a.z };
        viewPointSurface = normalize(viewPointSurface);
        float dotCamera = max(0.0, dotProduct(viewPointSurface,normal));
        if(dotCamera == 0.0) continue;

        vector3 v0, v1, v2;
        v0 = projectionFocal(vector0,HEIGTH,WIDTH,near,l, r, t, b, near);
        v1 = projectionFocal(vector1,HEIGTH,WIDTH,near,l, r, t, b, near);
        v2 = projectionFocal(vector2,HEIGTH,WIDTH,near,l, r, t, b, near);


        			float xmin = min3(v0.x, v1.x, v2.x);
        			float ymin = min3(v0.y, v1.y, v2.y);
        			float xmax = max3(v0.x, v1.x, v2.x);
        			float ymax = max3(v0.y, v1.y, v2.y);

        			if (xmin > HEIGTH - 1 || xmax < 0 || ymin > WIDTH - 1 || ymax < 0) continue;

        			uint32_t x0 = max((int32_t) 0, (int32_t)(floor(xmin)));
        			uint32_t x1 = min((int32_t)(WIDTH - 1), (int32_t)(floor(xmax)));
        			uint32_t y0 = max((int32_t) 0, (int32_t)(floor(ymin)));
        			uint32_t y1 = min((int32_t)(HEIGTH - 1), (int32_t)(floor(ymax)));
        			float area = edge(v0, v1, v2);

        			 for (uint32_t y = y0; y <= y1; ++y) {
        						for (uint32_t x = x0; x <= x1; ++x) {
        							vector3 sample = {x + 0.5,y + 0.5, 0};

        							float w0 = edge(v1, v2, sample);
        							float w1 = edge(v2, v0, sample);
        							float w2 = edge(v0, v1, sample);

        							if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
        	               w0 /= area;
        	               w1 /= area;
        	               w2 /= area;
                         float oneOverZ = v0.z * w0 + v1.z * w1 + v2.z * w2;
        	               float z = 1 / oneOverZ;

											if (z < depthBuffer[x][y]) {
												depthBuffer[x][y] = z;

												float px = (vector0.a.x/ vector0.a.z) * w0 + (vector1.a.x/ vector1.a.z) * w1 + (vector2.a.x / vector2.a.z) * w2;
												float py = (vector0.a.y/ vector0.a.z) * w0 + (vector1.a.y/ vector1.a.z) * w1 + (vector2.a.y / vector2.a.z) * w2;

												vector3 pt  = { -px * z, -py * z, z};
                        vector3 viewPoint = normalize(pt);
                        vector3 lightVector = minusVector3(pt, light);
                        float distanceFromLigth = distance(lightVector);
                        float attenuation = min(1.0 ,1 / pow(distanceFromLigth * ka + 1, 2));
                        vector3 normalLight = normalize(lightVector);


                        float dotLigth = dotProduct(normalLight, normal);
                        float diffuseLight =  max(0.0, dotLigth);
                        float specularCoefficient = 0;
                        if (diffuseLight > 0) {
                            float factorReflection = 2 * dotProduct(normalize(lightVector), normal);
                            vector3 tmpNormal = { factorReflection *  normal.x, factorReflection *  normal.y, factorReflection *  normal.z };
                            vector3 reflect = minusVector3( lightVector,tmpNormal);
                            reflect = normalize(reflect);
                            float cosAngle = max(0.0, dotProduct(viewPoint, reflect));
                            specularCoefficient = pow(cosAngle, materialShininess);
                        }

                        float finalColorCoff = ambientLight + (attenuation * (diffuseLight * kd + (specularCoefficient * ks)));
                        vector3 finalColor = { finalColorCoff, finalColorCoff, finalColorCoff};
                        finalColor.x *= materialColor.x;
                        finalColor.y *= materialColor.y;
                        finalColor.z *= materialColor.z;
												m[x][y].red =  min(finalColor.x, MAX_NUMBER);
											 	m[x][y].green = min(finalColor.y, MAX_NUMBER);
											 	m[x][y].blue = min(finalColor.z, MAX_NUMBER);

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
      triangle[index].normal = getNormal(triangle[index].lines[0].a,triangle[index].lines[1].a, triangle[index].lines[2].a);
			index++;
		}

    float t, b, l, r;

    computeScreenCoordinates(
           filmApertureWidth, filmApertureHeight,
           WIDTH, HEIGTH,
           1,
           nearClippingPLane,
           focalLength,
           &t, &b, &l, &r);

		int noFrame = 0;
		for(noFrame = 0; noFrame < FRAMES; noFrame++){
        createFrame(triangle,noFrame,numberObjects,3,l, r, t, b, nearClippingPLane);
		}

  return 0;
}
