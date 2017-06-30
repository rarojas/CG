#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "projection.h"

#define NI 4
#define NJ 3
#define RESOLUTIONI 100 * NI
#define RESOLUTIONJ 100 * NJ
#define  MAX_NUMBER 255
#define  WIDTH 720
#define  HEIGTH 480
#define MAGIC_NUMBER "P3"
#define FILE_NAME "projection.ppm"


Point inp[NI+1][NJ+1];
Point outp[RESOLUTIONI][RESOLUTIONJ];



double BezierBlend(int k, double mu, int n){
   int nn,kn,nkn;
   double blend=1;

   nn = n;
   kn = k;
   nkn = n - k;

   while (nn >= 1) {
      blend *= nn;
      nn--;
      if (kn > 1) {
         blend /= (double)kn;
         kn--;
      }
      if (nkn > 1) {
         blend /= (double)nkn;
         nkn--;
      }
   }
   if (k > 0)
      blend *= pow(mu,(double)k);
   if (n-k > 0)
      blend *= pow(1-mu,(double)(n-k));

   return(blend);
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


void createFrame(Triangle* triangles,int noFrame,int numberObjects,int length) {
		int  i ,j,index ;
		Pixel** m = raster(HEIGTH,WIDTH);

		double **depthBuffer = (double **) malloc(WIDTH * sizeof(double *));
		for(i = 0;i < WIDTH;i++) {
			depthBuffer[i] = (double *) malloc(HEIGTH * sizeof(double));
		}

		for(i = 0;i < WIDTH;i++)
			for(j = 0;j < HEIGTH;j++)
				depthBuffer[i][j] = INFINITY;

		double translateZ = 300;
    double translateX = 0;
    double translateY = -150;
    double angle = 0.1 * noFrame;
    double scaleFactor = 200;
		vector3 camera;

		vector3 light;
		light.x = 0;
		light.y = 0;
		light.z = 1;



		for(j = 0;j < numberObjects;j++) {
				Line vector0 =	translateLine(rotateLine(angle,scale(triangles[j].lines[0], scaleFactor)),translateX,translateY,translateZ);
				Line vector1 =	translateLine(rotateLine(angle,scale(triangles[j].lines[1], scaleFactor)),translateX,translateY,translateZ);
				Line vector2 =	translateLine(rotateLine(angle,scale(triangles[j].lines[2], scaleFactor)),translateX,translateY,translateZ);

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
		sprintf(str,"img/bezierSurface%d.ppm", noFrame);
		FILE *fp = createFile(str,HEIGTH,WIDTH,MAX_NUMBER,MAGIC_NUMBER);
		saveMatrix(fp, m,HEIGTH,WIDTH);
		free(m);
		free(depthBuffer);
		fclose(fp);
}




int main(argc,argv)
int argc;
char **argv;
{
   int i,j,ki,kj;
   double mui,muj,bi,bj;

   srandom(1111);
   for (i=0;i<=NI;i++) {
      for (j=0;j<=NJ;j++) {
         inp[i][j].x = i;
         inp[i][j].y = j;
         inp[i][j].z = (random() % 10000) / 5000.0 + 1;
      }
   }

   for (i=0;i<RESOLUTIONI;i++) {
      mui = i / (double)(RESOLUTIONI - 1);
      for (j=0;j<RESOLUTIONJ;j++) {
         muj = j / (double)(RESOLUTIONJ - 1);
         outp[i][j].x = 0;
         outp[i][j].y = 0;
         outp[i][j].z = 0;
         for (ki=0; ki <= NI; ki++) {
            bi = BezierBlend(ki, mui, NI);
            for (kj=0;kj<=NJ;kj++) {
               bj = BezierBlend(kj,muj,NJ);
               outp[i][j].x += (inp[ki][kj].x * bi * bj);
               outp[i][j].y += (inp[ki][kj].y * bi * bj);
               outp[i][j].z += (inp[ki][kj].z * bi * bj);
            }
         }
      }
   }

   Triangle* triangle = malloc(sizeof(Triangle) * 0);
   int numberObjects = 0, index = 0;

   for (i=0;i<RESOLUTIONI - 1;i++) {
      for (j=0;j<RESOLUTIONJ - 1;j++) {

        triangle = (Triangle *)realloc(triangle ,sizeof(Triangle) * ++numberObjects);
        triangle[index].lines = malloc(sizeof(Line) * 3);

        triangle[index].lines[0].a.x = outp[i][j].x;
        triangle[index].lines[0].a.y = outp[i][j].y;
        triangle[index].lines[0].a.z = outp[i][j].z;

        triangle[index].lines[1].a.x = outp[i][j + 1].x;
        triangle[index].lines[1].a.y = outp[i][j + 1].y;
        triangle[index].lines[1].a.z = outp[i][j + 1].z;

        triangle[index].lines[2].a.x = outp[i + 1][j + 1].x;
        triangle[index].lines[2].a.y = outp[i + 1][j + 1].y;
        triangle[index].lines[2].a.z = outp[i + 1][j + 1].z;

        index++;

        triangle = (Triangle *)realloc(triangle ,sizeof(Triangle) * ++numberObjects);
        triangle[index].lines = malloc(sizeof(Line) * 3);

        triangle[index].lines[0].a.x = outp[i + 1][j + 1].x;
        triangle[index].lines[0].a.y = outp[i + 1][j + 1].y;
        triangle[index].lines[0].a.z = outp[i + 1][j + 1].z;

        triangle[index].lines[1].a.x = outp[i + 1][j].x;
        triangle[index].lines[1].a.y = outp[i + 1][j].y;
        triangle[index].lines[1].a.z = outp[i + 1][j].z;

        triangle[index].lines[2].a.x = outp[i][j].x;
        triangle[index].lines[2].a.y = outp[i][j].y;
        triangle[index].lines[2].a.z = outp[i][j].z;

        index++;

      }
   }


   int noFrame = 0;
   for(noFrame = 0; noFrame < 10; noFrame++){
       createFrame(triangle,noFrame,numberObjects,3);
   }


}
