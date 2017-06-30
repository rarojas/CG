#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <OpenCL/opencl.h>
#include <string.h>
#include "projection.h"


#define NI 8
#define NJ 2
#define RESOLUTIONI 15 * NI
#define RESOLUTIONJ 15 * NJ
#define MAX_NUMBER 255
#define WIDTH 480
#define HEIGTH 480
#define MAGIC_NUMBER "P3"
#define FILE_NAME "projection.ppm"
#define FRAMES 10
#define MODE 1


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


float materialShininess = 40;

float ks = 0.1;
float kd = 1;
float ka = 0.05;
float ambientLight = 0.1;


Line rotateLine(Line line){
	double cosine = cos( M_PI );
	double sinus = sin( M_PI);
	double matrixRotation[4][4] = {
    {	cosine	,	 0, 	   sinus,      0   },
    { 0 ,        1 ,    0,   0 },
    { -sinus ,    0 ,    cosine,   0 },
		{	0	,			   0, 	   0,      1    },
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


void createFrame(Triangle* triangles,int noFrame,int numberObjects,int length,
  float l, float r, float t, float b, float near, vector3** vectorNormals
) {
		int  i ,j,index;
		Pixel** m = rasterBackground(HEIGTH,WIDTH,background);
		double **depthBuffer = createDeepBuffer(WIDTH, HEIGTH, INFINITY);
    double angle = 2 * M_PI * (float)noFrame / (float)FRAMES;
    vector3 light =  { cos(angle) * 3 ,   sin(angle) * 3,     sin(angle) * 2.5  - 2.5   };

		for(j = 0;j < numberObjects;j++) {
      Line vector0 =	triangles[j].lines[0];
      Line vector1 =	triangles[j].lines[1];
      Line vector2 =	triangles[j].lines[2];
      vector3 normal = triangles[j].normal;

      vector0 =	getPointToCamera(camera, lookAt, vector0);
      vector1 =	getPointToCamera(camera, lookAt, vector1);
      vector2 =	getPointToCamera(camera, lookAt, vector2);

      normal = getVectorToCamera(camera, lookAt, normal);

      vector3 normal0 = getVectorToCamera(camera, lookAt,vectorNormals[vector0.i][vector0.j]);
      vector3 normal1 = getVectorToCamera(camera, lookAt,vectorNormals[vector1.i][vector1.j]);
      vector3 normal2 = getVectorToCamera(camera, lookAt,vectorNormals[vector2.i][vector2.j]);

      normal0 = normalize(normal0);
      normal1 = normalize(normal1);
      normal2 = normalize(normal2);

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

			#if defined MODE && MODE == 3
			vector3 color[3];


			#endif

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

                      float z = INFINITY;
                      vector3 pt;

                      #if defined MODE && MODE == 0
                        float px = (vector0.a.x / vector0.a.z) * w0 + (vector1.a.x / vector1.a.z) * w1 + (vector2.a.x / vector2.a.z) * w2;
                        float py = (vector0.a.y / vector0.a.z) * w0 + (vector1.a.y / vector1.a.z) * w1 + (vector2.a.y / vector2.a.z) * w2;
  											float oneOverZ = v0.z * w0 + v1.z * w1 + v2.z * w2;
  	                    z = 1 / oneOverZ;
                        pt =  (vector3){ -px * z, -py * z, -z};
                      #endif

                      #if defined MODE && MODE == 1

                        float px = (vector0.a.x / vector0.a.z) * w0 + (vector1.a.x / vector1.a.z) * w1 + (vector2.a.x / vector2.a.z) * w2;
                        float py = (vector0.a.y / vector0.a.z) * w0 + (vector1.a.y / vector1.a.z) * w1 + (vector2.a.y / vector2.a.z) * w2;
                        float oneOverZ = v0.z * w0 + v1.z * w1 + v2.z * w2;
                        z = 1 / oneOverZ;
                        pt = (vector3){ -px * z, -py * z, -z};

                        vector3 normalInterpolated = interpolateVector(normal2, normal1, normal0, w1, w0, w2);
                        normal = normalize(normalInterpolated);


                      #endif

                      #if defined MODE && MODE == 2


                        Vertex vertex = interpolatePhong(normal2, normal1, normal0, vector2.a, vector1.a, vector0.a, w1, w0, w2);
                        normal =  normalize(vertex.m_normal);
                        if(triangles[j].even == 0) {
                          vertex = interpolatePhong(normal0, normal1, normal2, vector0.a, vector1.a, vector2.a, w1, w2, w0);
                          normal =  normalize(vertex.m_normal);
                        }

                        pt = (vector3) { vertex.m_pos.x, vertex.m_pos.y, vertex.m_pos.z};
                        z = pt.z;

                        float px = (vector0.a.x / vector0.a.z) * w0 + (vector1.a.x / vector1.a.z) * w1 + (vector2.a.x / vector2.a.z) * w2;
                        float py = (vector0.a.y / vector0.a.z) * w0 + (vector1.a.y / vector1.a.z) * w1 + (vector2.a.y / vector2.a.z) * w2;
                        float oneOverZ = v0.z * w0 + v1.z * w1 + v2.z * w2;
                        //z = 1 / oneOverZ;
                        //pt = (vector3){ px * z, py * z, z};
                        //pt.z = vertex.m_pos.z;
                      #endif

											#if defined MODE && MODE == 3


											#endif
											if (z < depthBuffer[x][y]) {
												depthBuffer[x][y] = z;

												vector3 viewPoint = normalize(pt);
                        vector3 lightVector = minusVector3(pt, light);
                        double distanceFromLigth = distance(lightVector);
                        double attenuation = min(1.0 ,1 / pow(distanceFromLigth * ka + 1, 2));
                        vector3 normalLight = normalize(lightVector);


                        double dotLigth = dotProduct(normalLight, normal);
                        double diffuseLight =  max(0.0, dotLigth);
                        double specularCoefficient = 0;

                        if (diffuseLight > 0) {
                            double factorReflection = 2 * dotProduct(normalLight, normal);
                            vector3 tmpNormal = { factorReflection *  normal.x, factorReflection *  normal.y, factorReflection *  normal.z };
                            vector3 reflect = minusVector3(tmpNormal, lightVector);
                            reflect = normalize(reflect);
                            double cosAngle = max(0.0, dotProduct(viewPoint, reflect));
                            specularCoefficient = pow(cosAngle, materialShininess);
                        }


                        double finalColorCoff = ambientLight + (attenuation * (diffuseLight * kd + (specularCoefficient * ks)));
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
		sprintf(str,"img/bezierSurface%d.ppm", noFrame);
		FILE *fp = createFile(str,HEIGTH,WIDTH,MAX_NUMBER,MAGIC_NUMBER);
		saveMatrix(fp, m,HEIGTH,WIDTH);
		free(m);
		free(depthBuffer);
		fclose(fp);
}


int main( int argc,char **argv){
   int i,j;

   dispatch_queue_t queue =
   #ifdef USE_CPU
      NULL;
   #else
      gcl_create_dispatch_queue(CL_DEVICE_TYPE_GPU, NULL);
   #endif
      if (queue == NULL) {
          queue = gcl_create_dispatch_queue(CL_DEVICE_TYPE_CPU, NULL);
          fprintf(stderr, "Warning: Running on CPU\n");
      }
      else {
          char name[128];
          cl_device_id gpu = gcl_get_device_id_with_dispatch_queue(queue);
          clGetDeviceInfo(gpu, CL_DEVICE_NAME, 128, name, NULL);
          fprintf(stderr, "Running on GPU %s\n", name);
      }



   Point** outp  = (Point **)malloc(sizeof(Point *) * RESOLUTIONI);
   int aux;
   for(aux = 0; aux < RESOLUTIONI; aux++)
      outp[aux] = (Point *) malloc(sizeof(Point) * RESOLUTIONJ);

   Point inp[NI + 1][NJ + 1] = {
     { {0, 1, 0}, { 0,  1, -1},{ 0,  0, -1} },
     { {0, 1, 0}, { -1, 1, -1},{ -1, 0, -1} },
     { {0, 1, 0}, { -1, 1, 0}, { -1, 0, 0} },
     { {0, 1, 0}, { -1, 1, 1}, { -1, 0, 1} },
     { {0, 1, 0}, { 0,  1, 1}, { 0,  0, 1} },
     { {0, 1, 0}, { 1,  1, 1}, { 1,  0, 1} },
     { {0, 1, 0}, { 1,  1, 0}, { 1,  0, 0} },
     { {0, 1, 0}, { 1,  1,-1}, { 1,  0,-1} },
     { {0, 1, 0}, { 0,  1,-1}, { 0,  0,-1} }
   };

   Point** in_p  = (Point **)malloc(sizeof(Point *) * (NI + 1));
   for(aux = 0; aux <= NI; aux++) in_p[aux] = (Point *) malloc(sizeof(Point) * (NJ + 1));
   for(i = 0; i <= NI; i++) for(j = 0; j <= NJ; j++) in_p[i][j] = inp[i][j];

   getPointFromBezierPoints(outp,RESOLUTIONI,RESOLUTIONJ, in_p, NI,NJ);

   Triangle* triangle = malloc(sizeof(Triangle) * 0);
   int numberObjects = 0, index = 0;

   vector3** verticesNormal  = (vector3**) malloc(sizeof(vector3 *) * RESOLUTIONI );
   for(i = 0; i < RESOLUTIONI * 2; i++)
      verticesNormal[i] = (vector3 * ) malloc(sizeof(vector3) * RESOLUTIONJ * 2);


   for (i = 0;i < RESOLUTIONI - 1;i++) {
      for (j = 0; j < RESOLUTIONJ - 1;j++) {
        triangle = (Triangle *)realloc(triangle ,sizeof(Triangle) * ++numberObjects);
        triangle[index].lines = malloc(sizeof(Line) * 3);
        triangle[index].lines[2].a = (Point){ outp[i][j].x, outp[i][j].y, outp[i][j].z };
        triangle[index].lines[1].a = (Point){ outp[i][j + 1].x, outp[i][j + 1].y, outp[i][j + 1].z };
        triangle[index].lines[0].a = (Point){ outp[i + 1][j + 1].x, outp[i + 1][j + 1].y, outp[i + 1][j + 1].z };

        triangle[index].lines[0].i = triangle[index].lines[1].i = triangle[index].lines[2].i = i;
        triangle[index].lines[0].j = triangle[index].lines[1].j = triangle[index].lines[2].j = j;
        triangle[index].lines[1].j++;
        triangle[index].lines[0].i++;
        triangle[index].lines[0].j++;


        vector3 normal = getNormal(triangle[index].lines[0].a, triangle[index].lines[1].a, triangle[index].lines[2].a);
        triangle[index].normal = normal;

        verticesNormal[i][j] = midVector(verticesNormal[i][j], normal);
        verticesNormal[i][j + 1] = midVector(verticesNormal[i][j + 1], normal);
        verticesNormal[i + 1][j + 1] = midVector(verticesNormal[i + 1][j + 1], normal);


        index++;
        triangle = (Triangle *)realloc(triangle ,sizeof(Triangle) * ++numberObjects);
        triangle[index].lines = malloc(sizeof(Line) * 3);

        triangle[index].lines[0].a = (Point){ outp[i][j].x, outp[i][j].y, outp[i][j].z };
        triangle[index].lines[1].a = (Point){ outp[i + 1][j].x, outp[i + 1][j].y, outp[i + 1][j].z };
        triangle[index].lines[2].a = (Point){ outp[i + 1][j + 1].x, outp[i + 1][j + 1].y, outp[i + 1][j + 1].z };

        triangle[index].lines[0].i = triangle[index].lines[1].i = triangle[index].lines[2].i = i;
        triangle[index].lines[0].j = triangle[index].lines[1].j = triangle[index].lines[2].j = j;
        triangle[index].lines[2].i = ++triangle[index].lines[1].i;
        triangle[index].lines[2].j++;
        triangle[index].even = 1;


        normal = getNormal(triangle[index].lines[0].a, triangle[index].lines[1].a, triangle[index].lines[2].a);
        triangle[index].normal = normal;

        verticesNormal[i][j] = midVector(verticesNormal[i][j], normal);
        verticesNormal[i + 1][j] = midVector(verticesNormal[i  + 1][j], normal);
        verticesNormal[i + 1][j + 1] = midVector(verticesNormal[i + 1][j + 1], normal);

        index++;

      }
   }

   for(j = 0;j < RESOLUTIONJ; j ++){
    vector3 normal = midVector(verticesNormal[0][j],verticesNormal[RESOLUTIONI - 2][j]);
    verticesNormal[0][j] = verticesNormal[RESOLUTIONI - 2][j]  = normal;
   }

   Point inp2[NI+1][NJ+1] = {
     { {0, -1, 0}, { 0,  -1, -1}, { 0,  0, -1} },
     { {0, -1, 0}, { -1, -1, -1}, { -1, 0, -1} },
     { {0, -1, 0}, { -1, -1,  0}, { -1, 0,  0} },
     { {0, -1, 0}, { -1, -1,  1}, { -1, 0,  1} },
     { {0, -1, 0}, { 0,  -1,  1}, { 0,  0,  1} },
     { {0, -1, 0}, { 1,  -1,  1}, { 1,  0,  1} },
     { {0, -1, 0}, { 1,  -1,  0}, { 1,  0,  0} },
     { {0, -1, 0}, { 1,  -1, -1}, { 1,  0, -1} },
     { {0, -1, 0}, { 0,  -1, -1}, { 0,  0, -1} }
   };

   for(i = 0; i <= NI; i++) for(j = 0; j <= NJ; j++) in_p[i][j] = inp2[i][j];
   getPointFromBezierPoints(outp,RESOLUTIONI,RESOLUTIONJ, in_p, NI,NJ);


   for (i = 0;i < RESOLUTIONI - 1;i++) {
      for (j = 0;j < RESOLUTIONJ - 1;j++) {
        triangle = (Triangle *)realloc(triangle ,sizeof(Triangle) * ++numberObjects);
        triangle[index].lines = malloc(sizeof(Line) * 3);

        triangle[index].lines[1].a = (Point){ outp[i + 1][j].x, outp[i + 1][j].y, outp[i + 1][j].z };
        triangle[index].lines[0].a = (Point){ outp[i][j + 1].x, outp[i][j + 1].y, outp[i][j + 1].z };
        triangle[index].lines[2].a = (Point){ outp[i][j].x, outp[i][j].y, outp[i][j].z };

        triangle[index].lines[0].i = triangle[index].lines[1].i = triangle[index].lines[2].i = i;
        triangle[index].lines[0].j = triangle[index].lines[1].j = triangle[index].lines[2].j = j + RESOLUTIONJ;
        triangle[index].lines[0].j++;
        triangle[index].lines[1].i++;

        vector3 normal = getNormal(triangle[index].lines[0].a, triangle[index].lines[1].a, triangle[index].lines[2].a);
        triangle[index].normal = normal;
        verticesNormal[i + 1][j + RESOLUTIONJ] = midVector(verticesNormal[i+ 1][j + RESOLUTIONJ], normal);
        verticesNormal[i][j + RESOLUTIONJ + 1] = midVector(verticesNormal[i][j + RESOLUTIONJ + 1], normal);
        verticesNormal[i][j + RESOLUTIONJ ] = midVector(verticesNormal[i][j + RESOLUTIONJ], normal);

        index++;

        triangle = (Triangle *)realloc(triangle ,sizeof(Triangle) * ++numberObjects);
        triangle[index].lines = malloc(sizeof(Line) * 3);

        triangle[index].lines[2].a = (Point){ outp[i + 1][j].x, outp[i + 1][j].y, outp[i + 1][j].z };
        triangle[index].lines[0].a = (Point){ outp[i][j + 1].x, outp[i][j + 1].y, outp[i][j + 1].z };
        triangle[index].lines[1].a = (Point){ outp[i + 1][j + 1].x, outp[i + 1][j + 1].y, outp[i + 1][j + 1].z };

        triangle[index].lines[0].i = triangle[index].lines[1].i = triangle[index].lines[2].i = i;
        triangle[index].lines[0].j = triangle[index].lines[1].j = triangle[index].lines[2].j = j + RESOLUTIONJ;
        triangle[index].lines[2].j++;
        triangle[index].lines[0].i++;
        triangle[index].lines[1].j++;
        triangle[index].lines[1].i++;
        triangle[index].even = 1;



        normal = getNormal(triangle[index].lines[0].a, triangle[index].lines[1].a, triangle[index].lines[2].a);
        triangle[index].normal = normal;
        verticesNormal[i + 1][j + RESOLUTIONJ] = midVector(verticesNormal[i + 1][j + RESOLUTIONJ], normal);
        verticesNormal[i][j + RESOLUTIONJ + 1] = midVector(verticesNormal[i][j + RESOLUTIONJ + 1], normal);
        verticesNormal[i + 1][j + RESOLUTIONJ + 1] = midVector(verticesNormal[i + 1][j + RESOLUTIONJ + 1], normal);


        index++;
      }
   }

   //for(i = 0;i < RESOLUTIONJ; i ++){
    //vector3 normal = midVector(verticesNormal[RESOLUTIONI][i + RESOLUTIONJ],verticesNormal[RESOLUTIONI * 2 - 1][i + RESOLUTIONJ]);
    //verticesNormal[RESOLUTIONI][i + RESOLUTIONJ] = verticesNormal[RESOLUTIONI * 2 - 1][i + RESOLUTIONJ]  = normal;
   //}

   float t, b, l, r;

   computeScreenCoordinates(
          filmApertureWidth, filmApertureHeight,
          WIDTH, HEIGTH,
          1,
          nearClippingPLane,
          focalLength,
          &t, &b, &l, &r);

   int noFrame = 0;
   for(noFrame = 0; noFrame < FRAMES; noFrame++) {
       createFrame(triangle,noFrame,numberObjects,3,l, r, t, b, nearClippingPLane,verticesNormal);
   }
}
