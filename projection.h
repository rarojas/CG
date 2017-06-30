static const float inchToMm = 25.4;


#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
      _a < _b ? _a : _b; })


double min3( double a,  double b,  double c) {
   return min(a, min(b, c));
}

double max3( double a,  double b,  double c) {
   return max(a, max(b, c));
}

typedef struct {
   double  x;
 	 double  y;
   double  z;
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
  int  i;
  int  j;
} Line;


typedef struct {
	int red;
	int green;
	int blue;
} Pixel;


typedef union {
  double v[2];
  struct
  {
    double x;
    double y;
  };
} vector2;

typedef union {
  double v[3];
  struct
  {
    double x;
    double y;
    double z;
  };
} vector3;


typedef union {
  double v[4];
  struct
  {
    double x;
    double y;
    double z;
    double w;
  };
} vector4;


typedef struct {
	Point a;
	Point b;
  Point c;
  Line* lines;
  int even;
  vector3 normal;
	Color color;
} Triangle;

typedef struct {
  vector3 m_pos;
  vector3 m_tex;
  vector3 m_normal;
} Vertex;

typedef struct {
  Vertex vertexs[3];
} Polygone;

typedef struct {
  Polygone* polygones;
  vector3 m_color;
} Object;


struct Matrix {
    int rows;
    int cols;
    double** data;
};
typedef struct Matrix Matrix;


void drawLine(Pixel **matrix, Line line,int HEIGTH,int WIDTH){
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
				if(yCord >= HEIGTH || yCord <= 0 || xCord >= WIDTH || xCord <= 0)
			 		return;
				matrix[xCord][yCord].red = line.color.red;
				matrix[xCord][yCord].green = line.color.green;
				matrix[xCord][yCord].blue = line.color.blue;
		}
}


vector3 minus(Point a, Point b) {
    vector3 result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return result;
}


void saveMatrix(FILE *fp, Pixel** matrix,int HEIGTH,int WIDTH){
	int row, col;
	for(col = 0; col < HEIGTH; col++){
		for(row = 0; row < WIDTH;  row++) {
			fprintf(fp,"%i %i %i ",matrix[row][col].red, matrix[row][col].green, matrix[row][col].blue);
	 	}
	 	fputs("\n", fp);
	}
}

FILE* createFile(char*  filename,int HEIGTH,int WIDTH,int MAX_NUMBER,char* MAGIC_NUMBER) {
	FILE *fp;
	fp = fopen(filename, "w+" );
	fprintf(fp, "%s\n# %s\n%i %i\n%i\n",MAGIC_NUMBER, filename, WIDTH,HEIGTH, MAX_NUMBER);
	return fp;
}


void putpixel (int xCord,int yCord,Pixel** matrix, int color,int HEIGTH, int WIDTH){
  if(yCord >= HEIGTH || yCord <= 0 || xCord >= WIDTH || xCord <= 0)
    return;
  matrix[xCord][yCord].red = color;
  matrix[xCord][yCord].green = color;
  matrix[xCord][yCord].blue = color;
}


Pixel** rasterBackground(int HEIGTH, int WIDTH,vector3 background) {
  int i, j;
	Pixel **matrix = (Pixel **) malloc(WIDTH * sizeof(Pixel *));
	for(i = 0;i < WIDTH;i++) {
		matrix[i] = (Pixel *) malloc(HEIGTH * sizeof(Pixel));
	}
	for(i = 0; i < WIDTH;  i++) {
		for(j = 0; j < HEIGTH; j++){
			matrix[i][j].red = background.x;
			matrix[i][j].green = background.y;
			matrix[i][j].blue = background.z;
		}
	}
	return matrix;
}



Pixel** raster(int HEIGTH, int WIDTH) {
  vector3 background  = { 255,255,255};
  return rasterBackground(HEIGTH,WIDTH,background);
}


double dotProduct(vector3 a, vector3 b) {
   return (a.x * b.x) + (a.y * b.y) + (a.z * a.z);
}

vector3 crossProduct(vector3 a,vector3 b) {
  vector3 result;
  result.x = (a.y * b.z) - (a.z * b.y);
  result.y = (a.z * b.x) - (a.x * b.z);
  result.z = (a.x * b.y) - (a.y * b.x);
  return result;
}

float norm(vector3 a) {
   return (a.x * a.x) + (a.y * a.y) + (a.z * a.z);
}

vector3 normalize(vector3 a) {
       float n = norm(a);
       if (n > 0) {
           float factor = 1 / sqrt(n);
           a.x *= factor, a.y *= factor, a.z *= factor;
       }
       return a;
   }

Matrix* make_matrix(int n_rows, int n_cols) {
    struct Matrix matrix;
    matrix.rows = n_rows;
    matrix.cols = n_cols;
    matrix.data = (double **)malloc(sizeof(double *) * n_rows);
    for(int x = 0; x < n_rows; x++){
        matrix.data[x] = (double *)calloc(n_cols, sizeof(double));
    }
    struct Matrix *m;
    m = &matrix;
    return m;
}

double edgeFunction(const Point a, const Point b, const Point c){
  return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
}


double edge(const vector3 a, const vector3 b, const vector3 c){
  return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
}


void multiply (double a[4][4],double b[4], double *resultmatrix){
    double c[4] = { b[0],b[1],b[2],b[3] };
    for(size_t i = 0; i < 4; i++) {
      double result = 0;
      for(size_t j = 0; j < 4; j++){
      	 result += a[i][j] * c[j];
      }
      resultmatrix[i] = result;
    }
}



void multiplyMatrix(Matrix* m,vector4 b, vector4 resultVector)
{   int i,j;
    for(i = 0; i <  m->rows; i++) {
      double result  = 0;
      for(j = 0; j < m->cols; j++)
        result += m->data[i][j] * b.v[j];
      resultVector.v[i] = result;
    }
}




Line scale(Line line, double scale){
		double matrixScale[4][4] = {
			{	scale, 		0,  	0,     0},
			{	0	,		scale, 		0,     0},
			{ 0 , 	0,    scale ,    0},
			{ 0 , 	0,    0    ,   	 1}
		 };
		double aMatrix[4] = { line.a.x , line.a.y, line.a.z ,1};
		multiply(matrixScale,aMatrix,(double *)&aMatrix);
		line.a.x = aMatrix[0];
		line.a.y = aMatrix[1];
		line.a.z = aMatrix[2];

		double bMatrix[4] = { line.b.x , line.b.y, line.b.z, 1};
		multiply(matrixScale,bMatrix, (double *)&bMatrix);
		line.b.x = bMatrix[0];
		line.b.y = bMatrix[1];
		line.b.z = bMatrix[2];
		return line;
}

Line translateLine(Line line, double x,double y,double z){
		double matrixTranlate[4][4] = {
			{	1, 		0,  	0,     x},
			{	0	,		1, 		0,     y},
			{ 0 , 	0,    1,     z},
			{ 0 , 	0,    0,  	 1}
		 };
		double aMatrix[4] = { line.a.x , line.a.y, line.a.z ,1};
		multiply(matrixTranlate,aMatrix,(double *)&aMatrix);
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

Line projection(Line line,int WIDTH,int HEIGTH){
  double f = 300;
  double matrixProjection[4][4] = {
			{	f, 		0,   WIDTH  /  2 , 	0 },
			{	0	,		f ,  HEIGTH /  2 ,  0 },
		 	{ 0 , 	0,   f   ,  0 },
			{ 0 , 	0,   0   ,  f },
	 };
	double oneOverZ = 1 / line.a.z;
  double aMatrix[4] = { line.a.x * oneOverZ, line.a.y  * oneOverZ, 1 , oneOverZ };
  double bMatrix[4] = { line.b.x / line.b.z, line.b.y / line.b.z, 1 , 1 / line.b.z };
  double* resultA = (double *)malloc(sizeof(double) * 4);
  multiply(matrixProjection,aMatrix, resultA);
  double* resultB = (double *)malloc(sizeof(double) * 4);
  multiply(matrixProjection,bMatrix, resultB);
  Line projection;
  projection.a.x = resultA[0];
  projection.a.y = resultA[1];
	projection.a.z = resultA[2];
  projection.b.x = resultB[0];
  projection.b.y = resultB[1];
	projection.b.z = resultA[2];
  free(resultA);
  free(resultB);
  return projection;
}

vector3 projectionFocal(Line line,int WIDTH,int HEIGTH,
     double f,
     float l,
     float r,
     float t,
     float b,
     float near) {
  	float oneOverZ = near / line.a.z;

    vector3 vertexScreen = { line.a.x * oneOverZ, line.a.y  * oneOverZ, 1 / line.a.z };
    vector2 vertexNDC;
    vertexNDC.x = 2 * vertexScreen.x / (r - l) - (r + l) / (r - l);
    vertexNDC.y = 2 * vertexScreen.y / (t - b) - (t + b) / (t - b);
    vector2 vertexRaster;
    vertexRaster.x = (vertexNDC.x + 1) / 2 * WIDTH;
    vertexRaster.y = (1 - vertexNDC.y) / 2 * HEIGTH;

    vector3 projection;
    projection.x =  vertexRaster.x;
    projection.y =  vertexRaster.y;
    projection.z =  -vertexScreen.z;
    return projection;
}


float distance(vector3 a){
  return sqrt(pow(a.x,2) + pow(a.y,2) + pow(a.z,2));
}

vector3 minusVector3(vector3 a, vector3 b){
  vector3 diff;
  diff.x = a.x - b.x;
  diff.y = a.y - b.y;
  diff.z = a.z - b.z;
  return diff;
}


void computeScreenCoordinates(
    float filmApertureWidth,
    float filmApertureHeight,
    int imageWidth,
    int imageHeight,
    int fitFilm,
    float nearClippingPLane,
    float focalLength,
    float *top, float *bottom, float *left, float *right) {

    float filmAspectRatio = filmApertureWidth / filmApertureHeight;
    float deviceAspectRatio = imageWidth / (float)imageHeight;

    *top = ((filmApertureHeight * inchToMm / 2) / focalLength) * nearClippingPLane;
    *right = ((filmApertureWidth * inchToMm / 2) / focalLength) * nearClippingPLane;
    float xscale = 1;
    float yscale = 1;

    switch (fitFilm) {
        default:
        case 0:
            if (filmAspectRatio > deviceAspectRatio) {
                xscale = deviceAspectRatio / filmAspectRatio;
            }
            else {
                yscale = filmAspectRatio / deviceAspectRatio;
            }
            break;
        case 1:
            if (filmAspectRatio > deviceAspectRatio) {
                yscale = filmAspectRatio / deviceAspectRatio;
            }
            else {
                xscale = deviceAspectRatio / filmAspectRatio;
            }
            break;
    }
    *right *= xscale;
    *top *= yscale;
    *bottom = -*top;
    *left = -*right;
}



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

void getPointFromBezierPoints(Point** outp,int out_i,int out_j,Point **inp, int in_i,int in_j){
    int i,j,ki,kj;
    double mui,muj,bi,bj;
    for (i = 0;i < out_i;i++) {
       mui = i / (double)(out_i - 1);
       for (j = 0;j < out_j; j++) {
          muj = j / (double)(out_j - 1);
          outp[i][j].x = 0;
          outp[i][j].y = 0;
          outp[i][j].z = 0;
          for (ki=0; ki <= in_i; ki++) {
             bi = BezierBlend(ki, mui, in_i);
             for (kj = 0; kj <= in_j; kj++) {
                bj = BezierBlend(kj,muj,in_j);
                outp[i][j].x += (inp[ki][kj].x * bi * bj);
                outp[i][j].y += (inp[ki][kj].y * bi * bj);
                outp[i][j].z += (inp[ki][kj].z * bi * bj);
             }
          }
       }
    }
}


double ** createDeepBuffer(int w, int h, float z) {
    int i,j;
    double **depthBuffer = (double **) malloc(w * sizeof(double *));
    for(i = 0;i < w;i++) depthBuffer[i] = (double *) malloc(h * sizeof(double));
    for(i = 0;i < w;i++) for(j = 0;j < h;j++) depthBuffer[i][j] = z;
    return depthBuffer;
}

double** lookAtCamera(vector3 from, vector3 to){

    vector3 tmp = { 0, 1, 0 };
    vector3 diff = { to.x - from.x, to.y - from.y, to.z - from.z };
    vector3 forward = normalize(diff);
    vector3 right = crossProduct(normalize(tmp), forward);
    vector3 up = crossProduct(forward, right);

    double** camToWorld = (double**)malloc(sizeof(double *) * 4);
    for(size_t i = 0; i  < 4;i++) camToWorld[i] = (double*)malloc(sizeof(double) * 4);

    camToWorld[0][0] = right.x;
    camToWorld[0][1] = right.y;
    camToWorld[0][2] = right.z;
    camToWorld[1][0] = up.x;
    camToWorld[1][1] = up.y;
    camToWorld[1][2] = up.z;
    camToWorld[2][0] = forward.x;
    camToWorld[2][1] = forward.y;
    camToWorld[2][2] = forward.z;

    camToWorld[3][0] = from.x;
    camToWorld[3][1] = from.y;
    camToWorld[3][2] = from.z;
    return camToWorld;
}


vector3 getVectorToCamera(vector3 camera,vector3 lookAt, vector3 vector){
    double dy = lookAt.y - camera.y;
    double dz = lookAt.z - camera.z;
    double dx = lookAt.x - camera.x;
    double d = sqrt(pow(dy, 2) + pow(dz, 2));
    double cosine = dz / d;
    double sinus =  dy / d;
    double worldToCamera [4][4] = {
      {   cosine,  0,  sinus,     dx   },
      {   0,       1,    0,       dy  },
      {   -sinus,  0,   cosine,   dz  },
      {   0,       0,    0,        1 },
    };

   double vertexB[4] =  { vector.x, vector.y, vector.z, 1};
   multiply((double(*)[4]) worldToCamera, vertexB, (double *)&vertexB);
   vector3 result = { vertexB[0], vertexB[1], vertexB[2] };
   return result;
}


Line getPointToCamera(vector3 camera,vector3 lookAt, Line vector){
    double dy = lookAt.y - camera.y;
    double dz = lookAt.z - camera.z;
    double dx = lookAt.x - camera.x;

    double d = sqrt(pow(dy, 2) + pow(dz, 2));
    double cosine = dz / d;
    double sinus =  dy / d;
    double worldToCamera [4][4] = {
      {   cosine,  0,  sinus,     dx   },
      {   0,       1,    0,       dy  },
      {   -sinus,  0,   cosine,   dz  },
      {   0,       0,    0,        1 },
    };

   double vertexB[4] =  { vector.a.x, vector.a.y, vector.a.z, 1};
   multiply((double(*)[4]) worldToCamera, vertexB, (double *)&vertexB);
   vector3 result = { vertexB[0], vertexB[1], vertexB[2] };
   vector.a.x = result.x;
   vector.a.y = result.y;
   vector.a.z = result.z;
   return vector;
}

vector3 midVector(vector3 a,vector3 b) {
  if(a.x == 0 && a.y == 0 &&  a.z == 0)
    return b;
  if(b.x == 0 && b.y == 0 &&  b.z == 0)
    return a;
  vector3 mid = { a.x + b.x, a.y + b.y ,a.z + b.z };
  mid.x /= 2;
  mid.y /= 2;
  mid.z /= 2;
  return normalize(mid);
}

vector3 getNormal(Point a, Point b,Point c) {
  vector3 aux0 = minus(b, a);
  vector3 aux1 = minus(c, a);
  return normalize(crossProduct(aux0, aux1));
}

vector3 multiplyVector(vector3 a, double factor) {
    a.x *= factor;
    a.y *= factor;
    a.z *= factor;
    return a;
}

vector3 lerp( vector3 a,vector3 b, float t){
  a = multiplyVector(a, t);
  b = multiplyVector(b, (1.f - t));
  vector3 inter = { a.x + b.x, a.y + b.y, a.z + b.z};
  return normalize(inter);
}


vector3 interpolateVector(vector3 a, vector3 b,vector3 c,float w0,float w1,float w2){
    return lerp(lerp(lerp(a, b, w0), b, w1), c, w2);
}


double getH1(double s) { return (2.0 * pow(s,3)) - (3.0 * s * s) + 1; }
double getH2(double s) { return (-2.0 * pow(s,3)) + (3.0 * s * s); }
double getH3(double s) { return pow(s,3) - (2.0 * s * s) + s; }
double getH4(double s) { return pow(s,3) - (s * s) ; }

Point hermite(Point a,vector3 at,Point b, vector3 bt, double s){
  double h1 = getH1(s);
  double h2 = getH2(s);
  double h3 = getH3(s) * 0.1;
  double h4 = getH4(s) * 0.1;

  double xt = (h1 * a.x) + (h2 * b.x) + (h3 * at.x) + (h4 * bt.x);
  double yt = (h1 * a.y) + (h2 * b.y) + (h3 * at.y) + (h4 * bt.y);
  double zt = (h1 * a.z) + (h2 * b.z) + (h3 * at.z) + (h4 * bt.z);
  Point interpolateVector = { xt, yt, zt };
  return interpolateVector;
}

Vertex interpolatePhong(vector3 a, vector3 b,vector3 c,Point ap, Point bp,Point cp,float w0,float w1,float w2) {

    Point aux0 = hermite(ap,a,bp,b,w0);
    vector3 aux0Normal = lerp(a, b, w0);

    Point aux1 = hermite(aux0,aux0Normal, cp,c,w1);
    vector3 aux1Normal = lerp(aux0Normal, c, w1);

    Point aux2 = hermite(aux1,aux1Normal,ap,a,w2);
    vector3 aux2Normal = lerp(aux1Normal, a, w2);

    vector3 m_pos = { aux2.x, aux2.y, aux2.z};

    Vertex vertex;
    vertex.m_normal = aux2Normal;
    vertex.m_pos = m_pos;
    return vertex;
}


Line rotateLineY(Line line,double angle){
	double cosine = cos( angle );
	double sinus = sin( angle);
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


Line translateTo(Line line, vector3 translateVector){
		line.a.x += translateVector.x;
		line.a.y += translateVector.y;
		line.a.z += translateVector.z;
    return line;
}
