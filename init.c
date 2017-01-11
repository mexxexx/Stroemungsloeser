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
	double *eps, double *omg, double *alpha, double *Re, double *Pr, double *beta, double *GX, double *GY, double *UI, double *VI, double *PI, 
	double *TI, int *wl, int *wr, int *wt, int *wb, double *posx1, double *posx2, double *posy1, double *posy2,
	int *tl, double *tl_value, int *tr, double *tr_value, int *tt, double *tt_value, int *tb, double *tb_value) {
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
	if (fscanf(f, "%*s = %lf\n", Pr) != 1)
		parameterError("Pr");
	if (fscanf(f, "%*s = %lf\n", beta) != 1)
		parameterError("beta");
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
	if (fscanf(f, "%*s = %lf\n", TI) != 1)
		parameterError("TI");
	if (fscanf(f, "%*s = %i %i %lf\n", wl, tl, tl_value) != 3)
		parameterError("wl");
	if (fscanf(f, "%*s = %i %i %lf\n", wr, tr, tr_value) != 3)
		parameterError("wr");
	if (fscanf(f, "%*s = %i %i %lf\n", wt, tt, tt_value) != 3)
		parameterError("wt");
	if (fscanf(f, "%*s = %i %i %lf\n", wb, tb, tb_value) != 3)
		parameterError("wb");
	if (fscanf(f, "%*s = %lf\n", posx1) != 1)
		parameterError("posx1");
	if (fscanf(f, "%*s = %lf\n", posy1) != 1)
		parameterError("posy1");
	if (fscanf(f, "%*s = %lf\n", posx2) != 1)
		parameterError("posx2");
	if (fscanf(f, "%*s = %lf\n", posy2) != 1)
		parameterError("posy2");
		
	fclose(f);
}

void initField(double *field, int imax, int jmax, double value) {
	for (int j = 0; j <= jmax + 1; j++) 
		for (int i = 0; i <= imax + 1; i++)
			field[POS2D(i, j, imax + 2)] = value;
}

unsigned char *loadBitmapFile(char *filename) {
    FILE *filePtr; 
    unsigned char *bitmapImage; 
    
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
	
    if (Header.Type !=0x4D42)
    {
        fclose(filePtr);
        return NULL;
    }

    fseek(filePtr, Header.OffBits, SEEK_SET);

    int imageSize = Header.biWidth * Header.biHeight * 3;
    bitmapImage = (unsigned char*)malloc(imageSize);

    if (!bitmapImage)
    {
        free(bitmapImage);
        fclose(filePtr);
        return NULL;
    }
    
    int padding = (4-(Header.biWidth * 3)%4)%4;
    
    for (int pixel = 0; pixel < imageSize/3; pixel++) 
    {
        fread(&bitmapImage[3*pixel+1], 1, 1, filePtr);
        fread(&bitmapImage[3*pixel+2], 1, 1, filePtr);
        fread(&bitmapImage[3*pixel], 1, 1, filePtr);
        
		if ((pixel+1)%Header.biWidth==0)
			fseek(filePtr, padding, SEEK_CUR);
    }
    fclose(filePtr);
    return bitmapImage;
}

void initFlag(char *heightMap, char *FLAG, int imax, int jmax, int *numFluidCells) {
	unsigned char *data = loadBitmapFile(heightMap);
	if (data == NULL) {
		printf("Error reading Height Map");
		free(FLAG);
		exit(EXIT_FAILURE);
	}
	
	(*numFluidCells)=0;
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			if (!data[POS2D(i-1, j-1, imax)*3]) 
				FLAG[POS2D(i, j, imax+2)] = 1;
			else {
				FLAG[POS2D(i, j, imax+2)] = 0;
				(*numFluidCells)++;
			}
		}
	}
	free(data);
	
	for (int i = 0; i <= imax+1; i++) {
		FLAG[POS2D(i, 0, imax+2)]=0;
		FLAG[POS2D(i, jmax+1, imax+2)]=0;
	}
	
	for (int j = 0; j <= jmax+1; j++) {
		FLAG[POS2D(0, j, imax+2)]=0;
		FLAG[POS2D(imax+1, j, imax+2)]=0;
	}
	
	for (int i = 1; i <= imax; i++) {  
		for (int j = 1; j <= jmax; j++) {
			int posij = POS2D(i, j, imax+2);
			if (FLAG[posij]) {		
				if (i==1)
					FLAG[posij-1] = 9;
				else if (i==imax)
					FLAG[posij+1] = 17;
				if (j==1)
					FLAG[posij-imax-2] = 5;
				else if (j==jmax)
					FLAG[posij+imax+2] = 3;
					
				if (FLAG[posij+1]) //OSTEN
					FLAG[posij] |= 17;
				if (FLAG[posij-1]) //WESTEN
					FLAG[posij] |= 9;
				if (FLAG[posij+imax+2]) //NORDEN
					FLAG[posij] |= 3;
				if (FLAG[posij-imax-2]) //SÜDEN
					FLAG[posij] |= 5;
			}
		}
	}
	
	for (int i = 1; i <= imax; i++) {
		if (FLAG[POS2D(i+1, 0, imax+2)]) //OSTEN
			FLAG[POS2D(i, 0, imax+2)] |= 17;
		if (FLAG[POS2D(i-1, 0, imax+2)]) //WESTEN
			FLAG[POS2D(i, 0, imax+2)] |= 9;
		if (FLAG[POS2D(i, 1, imax+2)]) //NORDEN
			FLAG[POS2D(i, 0, imax+2)] |= 3;
			
		if (FLAG[POS2D(i+1, jmax+1, imax+2)]) //OSTEN
			FLAG[POS2D(i, jmax+1, imax+2)] |= 17;
		if (FLAG[POS2D(i-1, jmax+1, imax+2)]) //WESTEN
			FLAG[POS2D(i, jmax+1, imax+2)] |= 9;
		if (FLAG[POS2D(i, jmax, imax+2)]) //SÜDEN
			FLAG[POS2D(i, jmax+1, imax+2)] |= 5;
		
	}
	
	for (int j = 1; j <= jmax; j++) {
		if (FLAG[POS2D(1, j, imax+2)]) //OSTEN
			FLAG[POS2D(0, j, imax+2)] |= 17;
		if (FLAG[POS2D(0, j-1, imax+2)]) //SÜDEN
			FLAG[POS2D(0, j, imax+2)] |= 5;
		if (FLAG[POS2D(0, j+1, imax+2)]) //NORDEN
			FLAG[POS2D(0, j, imax+2)] |= 3;
			
		if (FLAG[POS2D(imax, j, imax+2)]) //WESTEN
			FLAG[POS2D(imax+1, j, imax+2)] |= 9;
		if (FLAG[POS2D(imax+1, j-1, imax+2)]) //SÜDEN
			FLAG[POS2D(imax+1, j, imax+2)] |= 5;
		if (FLAG[POS2D(imax+1, j+1, imax+2)]) //NORDEN
			FLAG[POS2D(imax+1, j, imax+2)] |= 3;
	}
}
