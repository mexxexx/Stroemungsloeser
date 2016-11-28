#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "init.h"

#pragma pack(push, 1)

struct BitMap
{
  int16_t Type;
  int32_t Size;
  int16_t Reserve1;
  int16_t Reserve2;
  int32_t OffBits;
  int32_t biSize;
  int32_t biWidth;
  int32_t biHeight;
  int16_t biPlanes;
  int32_t biBitCount;
  int32_t biCompression;
  int32_t biSizeImage;
  int32_t biXPelsPerMeter;
  int32_t biYPelsPerMeter;
  int32_t biClrUsed;
  int32_t biClrImportant;
} Header;

int allocateVector(double **vector, int size) {
	*vector = (double*)malloc(size * sizeof(double));
	if (*vector == NULL) {
		printf("Konnte keinen Speicher allokieren.\n");
		return 1;
	}

	return 0;
}

void parameterError(char *param) {
	printf("Error reading parameter %s\n", param);
	exit(EXIT_FAILURE);
}

void readParameter(char *filename, char *simulationName, char *heightMap, double *xlength, double *ylength, int *imax, int *jmax, 
	double *delx, double *dely, double *delt, double *del_vec, double *t_end, double *tau, int *itermax,
	double *eps, double *omg, double *alpha, double *Re, double *GX, double *GY, double *UI, double *VI, double *PI) {
	FILE *f = fopen(filename, "r");
	if (f == NULL) {
		printf("Error opening parameter file\n");
		exit(EXIT_FAILURE);
	}
	
	if (fscanf(f, "%*s = %s\n", simulationName) != 1)
		parameterError("simulationName");
	if (fscanf(f, "%*s = %s\n", heightMap) != 1)
		parameterError("heightMap");
	
	if (fscanf(f, "%*s = %lf\n", xlength) != 1)
		parameterError("xlength");
	if (fscanf(f, "%*s = %lf\n", ylength) != 1)
		parameterError("ylength");
	if (fscanf(f, "%*s = %i\n", imax) != 1)
		parameterError("imax");
	if (fscanf(f, "%*s = %i\n", jmax) != 1)
		parameterError("jmax");
	(*delx) = (*xlength) / (*imax);	
	(*dely) = (*ylength) / (*jmax);	
		
	if (fscanf(f, "%*s = %lf\n", delt) != 1)
		parameterError("delt");
	if (fscanf(f, "%*s = %lf\n", t_end) != 1)
		parameterError("t_end");
	if (fscanf(f, "%*s = %lf\n", del_vec) != 1)
		parameterError("del_vec");
	if (fscanf(f, "%*s = %lf\n", tau) != 1)
		parameterError("tau");
	if (fscanf(f, "%*s = %i\n", itermax) != 1)
		parameterError("itermax");
	if (fscanf(f, "%*s = %lf\n", eps) != 1)
		parameterError("eps");
	if (fscanf(f, "%*s = %lf\n", omg) != 1)
		parameterError("omg");
	if (fscanf(f, "%*s = %lf\n", alpha) != 1)
		parameterError("alpha");
	if (fscanf(f, "%*s = %lf\n", Re) != 1)
		parameterError("Re");
	if (fscanf(f, "%*s = %lf\n", GX) != 1)
		parameterError("GX");
	if (fscanf(f, "%*s = %lf\n", GY) != 1)
		parameterError("GY");
	if (fscanf(f, "%*s = %lf\n", UI) != 1)
		parameterError("UI");
	if (fscanf(f, "%*s = %lf\n", VI) != 1)
		parameterError("VI");
	if (fscanf(f, "%*s = %lf\n", PI) != 1)
		parameterError("PI");
		
	fclose(f);
}

void initField(double *field, int imax, int jmax, double value) {
	for (int j = 0; j <= jmax + 1; j++) 
		for (int i = 0; i <= imax + 1; i++)
			field[POS2D(i, j, imax + 2)] = value;
}

unsigned char *loadBitmapFile(char *filename)
{
    FILE *filePtr; //our file pointer
    //BITMAPFILEHEADER bitmapFileHeader; //our bitmap file header
    unsigned char *bitmapImage;  //store image data
    int imageIdx=0;  //image index counter
    unsigned char tempRGB;  //our swap variable

    //open filename in read binary mode
    filePtr = fopen(filename,"rb");
    if (filePtr == NULL)
        return NULL;
	memset(&Header, 0, sizeof(Header));
   
	fread(&Header.Type, 2, 1, filePtr);
	fread(&Header.Size, 4, 1, filePtr);
	fread(&Header.Reserve1, 2, 1, filePtr);
	fread(&Header.Reserve2, 2, 1, filePtr);
	fread(&Header.OffBits, 4, 1, filePtr);
	fread(&Header.biSize, 4, 1, filePtr);
	fread(&Header.biWidth, 4, 1, filePtr);
	fread(&Header.biHeight, 4, 1, filePtr);
	fread(&Header.biPlanes, 2, 1, filePtr);
	fread(&Header.biBitCount, 2, 1, filePtr);
	fread(&Header.biCompression, 4, 1, filePtr);
	fread(&Header.biSizeImage, 4, 1, filePtr);
	fread(&Header.biXPelsPerMeter, 4, 1, filePtr);
	fread(&Header.biYPelsPerMeter, 4, 1, filePtr);
	fread(&Header.biClrUsed, 4, 1, filePtr);
	fread(&Header.biClrImportant, 4, 1, filePtr);
	printf("%i, %i\n", Header.biHeight, Header.biWidth);
	printf("%i\n", Header.biSizeImage);
	
	//verify that this is a bmp file by check bitmap id
    if (Header.Type !=0x4D42)
    {
        fclose(filePtr);
        return NULL;
    }

    //move file point to the begging of bitmap data
    fseek(filePtr, Header.OffBits, SEEK_SET);

    //allocate enough memory for the bitmap image data
    bitmapImage = (unsigned char*)malloc(Header.biSizeImage);

    //verify memory allocation
    if (!bitmapImage)
    {
        free(bitmapImage);
        fclose(filePtr);
        return NULL;
    }

    //read in the bitmap image data
    fread(bitmapImage, sizeof(unsigned char), Header.biSizeImage, filePtr);

    //make sure bitmap image data was read
    if (bitmapImage == NULL)
    {
        fclose(filePtr);
        return NULL;
    }

    //swap the r and b values to get RGB (bitmap is BGR)
    for (imageIdx = 0;imageIdx < Header.biSizeImage;imageIdx+=3) // fixed semicolon
    {
        tempRGB = bitmapImage[imageIdx];
        bitmapImage[imageIdx] = bitmapImage[imageIdx + 2];
        bitmapImage[imageIdx + 2] = tempRGB;
    }
	
	for (int imageIdx = 0; imageIdx < Header.biSizeImage; imageIdx+=3) 
		printf ("%d %d %d\n", bitmapImage[imageIdx], bitmapImage[imageIdx+1], bitmapImage[imageIdx+2]);

    //close file and return bitmap iamge data
    fclose(filePtr);
    return bitmapImage;
}

void initFlag(char *heightMap, char *FLAG, int imax, int jmax) {
	unsigned char *data = loadBitmapFile(heightMap);
	if (data == NULL) {
		printf("Error reading Height Map");
		free(FLAG);
		exit(EXIT_FAILURE);
	}
	
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= imax; j++) {
			FLAG[POS2D(i, j, imax+2)] = (data[POS2D(i * 3, j * 3, imax * 3)] == 0) ? 1 : 0;
		}
	}
	/*for (int j = 1; j <= imax; j++) {
		for (int i = 1; i <= imax; i++) {
			printf ("%d ", data[(i-1) * 3 + (j-1) * imax * 3]);
		}
		printf("\n");
	}*/
	
	free(data);
}
