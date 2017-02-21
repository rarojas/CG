#include <stdlib.h>
#include <stdio.h>
#define  MAX_NUMBER 255
#define  WIDTH 1920
#define  HEIGTH 1080
#define MAGIC_NUMBER "P3"
#define FILE_NAME "image.ppm"

typedef struct {
	int red;
	int green;
	int blue;
} Pixel;

int main(){
	int i,j;
	FILE *fp;
	fp = fopen(FILE_NAME, "w+" );
	fprintf(fp, "%s\n# %s\n%i %i\n%i\n",MAGIC_NUMBER, FILE_NAME, WIDTH,HEIGTH, MAX_NUMBER);
	for(j = 0; j < HEIGTH; j++) {
	   	for(i = 0 ; i < WIDTH ; i++ ) {
	   		struct Pixel pixel;
	   		pixel.red = rand() % MAX_NUMBER;
	   		pixel.green = rand() % MAX_NUMBER;
	   		pixel.blue =  rand() % MAX_NUMBER;
   			fprintf(fp,"%i %i %i ",pixel.red, pixel.green, pixel.blue);
	   }
	   fputs("\n", fp);
	}

	fclose(fp);
  return(0);
}
