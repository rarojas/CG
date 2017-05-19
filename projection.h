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


typedef struct{
	Point a;
	Point b;
  Point c;
  Line* lines;
	Color color;
} Triangle;


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


Pixel** raster(int HEIGTH, int WIDTH) {
  int i, j;

	Pixel **matrix = (Pixel **) malloc(WIDTH * sizeof(Pixel *));
	for(i = 0;i < WIDTH;i++) {
		matrix[i] = (Pixel *) malloc(HEIGTH * sizeof(Pixel));
	}
	for(i = 0; i < WIDTH;  i++) {
		for(j = 0; j < HEIGTH; j++){
			matrix[i][j].red = 255;
			matrix[i][j].green = 255;
			matrix[i][j].blue = 255;
		}
	}
	return matrix;
}


float dotProduct(vector3 a, vector3 b) {
   return a.x * b.x + a.y * b.y + a.z * a.z;
}

vector3 crossProduct(vector3 a,vector3 b) {
  vector3 result;
  result.x = a.y * b.z - a.z * b.y;
  result.y = a.z * b.x - a.x * b.z;
  result.z = a.x * b.y - a.y * b.x;
  return result;
}

float norm(vector3 a) {
   return a.x * a.x + a.y * a.y + a.z * a.z;
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




void multiply (double a[4][4],double b[4], double *resultmatrix)
{   int i,j;
    for(i = 0; i < 4; i++) {
      double result  = 0;
      for(j = 0; j < 4; j++)
      	result += a[i][j] * b[j];
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
		double scale = 1;
		double matrixTranlate[4][4] = {
			{	scale, 		0,  	0,     x},
			{	0	,		scale, 		0,     y},
			{ 0 , 	0,    scale ,    z},
			{ 0 , 	0,    0    ,   	 1}
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

  double* resultA = malloc(sizeof(double) * 4);
  multiply(matrixProjection,aMatrix, resultA);

  double* resultB = malloc(sizeof(double) * 4);
  multiply(matrixProjection,bMatrix, resultB);

  Line projection;
  projection.a.x = resultA[0];
  projection.a.y = resultA[1];
	projection.a.z = resultA[2];

  projection.b.x = resultB[0];
  projection.b.y = resultB[1];
	projection.a.z = resultA[2];


  free(resultA);
  free(resultB);
  return projection;
}
