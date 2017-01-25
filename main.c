#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "boundary.h"
#include "init.h"
#include "uvp.h"
#include "visual.h"
#include "particles.h"

#define PROGRESS_DISPLAY 0.05
#define PARTICLE_DELTA 0.05

int currentProgressStep = 1;
time_t start;

char simulationName[256], problem[256];
int imax, jmax, itermax, wl, wr, wt, wb, numFluidCells, partCount, anzahl, tl, tr, tt, tb;
double xlength, ylength, delx, dely, delt, delt_min, t_end, del_vec, tau, eps, omg, alpha, Re, Pr, beta, GX, GY, UI, VI, PI, TI, posx1, posx2, posy1, posy2;
double tl_value, tr_value, tt_value, tb_value;

void createSimulationDirectory() {
	struct stat st = {0};
	if (stat(simulationName, &st) == -1) {
		int result;
		#ifdef __MSDOS__
		result = mkdir(simulationName);
		#else
		result = mkdir(simulationName, 0700);
		#endif
		
		if (result) {
			printf("Verzeichnis für die Simulation konnte nicht angelegt werden\n");
			exit(EXIT_FAILURE);
		}
	}
}


void printSnapshot(int currentSnapshot, double *U, double *V, double *P, double *TEMP, 
					double *PSI, double *ZETA, double *HEAT, Particle *particles, int partCount, double t) {

	char nameOfVelocity[256];
	char nameOfPressure[256];
	char nameOfParticles[256];
	char nameOfPSI[256];
	char nameOfZETA[256];
	char nameOfHEAT[256];
	char nameOfTEMP[256];
	sprintf(nameOfVelocity, "%s/%s_velocity_%05d.vtk", simulationName, simulationName, currentSnapshot);
	sprintf(nameOfParticles, "%s/%s_particles_%05d.vtk", simulationName, simulationName, currentSnapshot);
	sprintf(nameOfPressure, "%s/%s_pressure_%05d.vtk", simulationName, simulationName, currentSnapshot++);
	sprintf(nameOfPSI, "%s/%s_PSI_%05d.vtk", simulationName, simulationName, currentSnapshot);
	sprintf(nameOfZETA, "%s/%s_ZETA_%05d.vtk", simulationName, simulationName, currentSnapshot);
	sprintf(nameOfHEAT, "%s/%s_HEAT_%05d.vtk", simulationName, simulationName, currentSnapshot);	
	sprintf(nameOfTEMP, "%s/%s_TEMP_%05d.vtk", simulationName, simulationName, currentSnapshot);	
	printVectorField(U, V, imax, jmax, xlength, ylength, nameOfVelocity); 
	printScalarField(P, imax, jmax, xlength, ylength, nameOfPressure);
	printParticles(particles, partCount, nameOfParticles);
	printScalarField(PSI, imax, jmax, xlength, ylength, nameOfPSI);
	printScalarField(ZETA, imax, jmax, xlength, ylength, nameOfZETA);
	printScalarField(HEAT, imax, jmax, xlength, ylength, nameOfHEAT);
	printScalarField(TEMP, imax, jmax, xlength, ylength, nameOfTEMP);	

	double progress = t/t_end; 
	if (progress >= PROGRESS_DISPLAY * currentProgressStep) {
		printf("%.1f%% completed\n", progress * 100);		
		currentProgressStep++;		
	}
}

void printField(double *field, int imax, int jmax) {
	for (int j = jmax+1; j >= 0; j--) {
		for (int i = 0; i <= imax+1; i++)
			printf("%f ", field[POS2D(i, j, imax+2)]);
		printf("\n");
	}
}

void printCharField(char *field, int imax, int jmax) {
	for (int j = jmax+1; j >= 0; j--) {
		for (int i = 0; i <= imax+1; i++)
			printf("%d ", field[POS2D(i, j, imax+2)]);
		printf("\n");
	}
}

int calculateFluidDynamics(double* U, double* V, double* P, double *TEMP, char* FLAG, double *PSI, double *ZETA, 
							double *HEAT, Particle *particles, int partCount) {
	initField(U, imax, jmax, UI);
	if (strcmp("Stufe", problem) == 0) {
		for (int i = 1; i <= imax; i++) 
			for (int j = 0; j <= jmax/2; j++) 
				U[POS2D(i,j,imax+2)]=0;
	}
	
	initField(V, imax, jmax, VI);
	initField(P, imax, jmax, PI);
	initField(TEMP, imax, jmax, TI);
	
	double *F, *G, *rhs;
	if (allocateVector(&F, (imax+2) * (jmax+2)))
		return 1;
	if (allocateVector(&G, (imax+2) * (jmax+2))) {
		free(F);
		return 1;
	}
	if (allocateVector(&rhs, (imax+2) * (jmax+2))) {
		free(F);
		free(G);
		return 1;
	}		
			
	int currentSnapshot = 1;
	double t = 0;
	double frameDuration = 0;
	double seedTime = 0;
	double umax = (fabs(UI) > eps) ? UI : eps;
	double vmax = (fabs(VI) > eps) ? VI : eps;
	start = time(NULL);
	
	int iter;
	double res;
	
	while (t < t_end) {
		computeDelt(&delt, delt_min, imax, jmax, delx, dely, umax, vmax, Re, Pr, tau);
		setBoundaryCond(U, V, TEMP, FLAG, delx, dely, imax, jmax, wl, wr, wt, wb, tl, tl_value, tr, tr_value, tt, tt_value, tb, tb_value);
		setSpecialBoundaryCond(U, V, TEMP, imax, jmax, problem);
		computeTEMP(U, V, TEMP, FLAG, imax, jmax, delt, delx, dely, alpha, Re, Pr);
		computeFG(U, V, F, G, TEMP, FLAG, imax, jmax, delt, delx, dely, GX, GY, alpha, Re, beta);
		computeRHS(F, G, rhs, FLAG, imax, jmax, delt, delx, dely);
		solvePoisson(P, rhs, FLAG, omg, eps, itermax, delx, dely, imax, jmax, numFluidCells, &iter, &res);
		adapUV(U, V, F, G, P, FLAG, imax, jmax, delt, delx, dely, &umax, &vmax); 
		COMP_PSI_ZETA(U, V, imax, jmax, xlength, ylength, PSI, ZETA, FLAG);
		COMP_HEAT(U, V, TEMP, HEAT, FLAG, Re, Pr, imax, jmax, delx, dely);
		
		if (partCount > 0) {
			if (seedTime >= PARTICLE_DELTA) {
				particleSeed(particles, posx1, posy1, posx2, posy2, partCount, anzahl);
				seedTime -= PARTICLE_DELTA;
			}
			particleVelocity(U, V, delx, dely, imax, jmax, particles, partCount);
			particleTransport(particles, delt, partCount, xlength, ylength);
		}

		if (frameDuration >= del_vec) {
			printf("t: %f, delt: %f, iter: %i, res: %f\n", t, delt, iter, res);
			printSnapshot(currentSnapshot++, U, V, P, TEMP, PSI, ZETA, HEAT, particles, partCount, t);
			frameDuration -= del_vec;
		}
		t += delt;
		frameDuration += delt;
		seedTime += delt;
	}
	
	printSnapshot(currentSnapshot, U, V, P, TEMP, PSI, ZETA, HEAT, particles, partCount, t_end);
	printf("%i frames taken in a time period of %.2f seconds (%.1f FPS)\n", currentSnapshot, t_end, currentSnapshot/t_end);
	printf("Duration: %.1f seconds\n", (double)(time(NULL)-start)); 
	free(F);
	free(G);
	free(rhs);
	return 0;
}

void readParams(int argc, char** argv, char* file) {int fileLoaded = 0;
	int problemSelected = 0;
	partCount = 20000;
	for (int i = 1; i < argc; i++) {
		switch (argv[i][0]) {
			case 'f':
				argv[i]+=2;
				sprintf(file, "%s", argv[i]);
				fileLoaded = 1;
				break;
			case 'p':
				argv[i]+=2;
				sprintf(problem, "%s", argv[i]);
				problemSelected = 1;
				break;
			case 'n':
				argv[i]+=2;
				sscanf(argv[i], "%i", &partCount); 
			default:
				break;
		}
	}
	
	if (!fileLoaded) {
		printf("Bitte eine Datei auswählen: ");
		scanf("%s", file);
	}
	
	if (!problemSelected){
		printf("Bitte ein Problem angeben: ");
		scanf("%s", problem);
	}	
		
	if (partCount < 200) 
		partCount = 0;
}

int main(int argc, char** argv) {
	char file[256];
	readParams(argc, argv, file);
	
	Particle *particles = malloc(sizeof(Particle)*partCount);
	if (particles == NULL) {
		printf("Konnte keinen Speicherplatz fuer particles allokieren");
		return 1;
	}
	
	particleInit(particles, partCount);
	char obstacelsMap[256];
	readParameter(file, simulationName, obstacelsMap, &xlength, &ylength, &imax, &jmax, &delx, &dely, &delt_min, 
							&del_vec, &t_end, &tau, &itermax, &eps, &omg, &alpha, &Re, &Pr, &beta, &GX, &GY, 
							&UI, &VI, &PI, &TI, &wl, &wr, &wt, &wb, &posx1, &posx2, &posy1, &posy2,
							&tl, &tl_value, &tr, &tr_value, &tt, &tt_value, &tb, &tb_value);
	
	createSimulationDirectory();

	anzahl = (int)((sqrt((posx1-posx2)*(posx1-posx2)+(posy1-posy2)*(posy1-posy2)))*10)+1;
	
	char *FLAG = (char*)malloc((imax+2) * (jmax+2) * sizeof(char));
	if (FLAG == NULL) {
		free(particles);
		printf("Konnte keinen Speicherplatz fuer FLAG allokieren");
		return 1;
	}
	initFlag(obstacelsMap, FLAG, imax, jmax, &numFluidCells);
	bmp_verify(FLAG,imax, jmax);

	char obstacleFile[256];
	sprintf(obstacleFile, "%s/obstacles.vtk", simulationName);
	printObstacles(FLAG, imax, jmax, xlength, ylength, obstacleFile);
	
	double *U, *V, *P, *TEMP, *PSI, *ZETA, *HEAT;
	if (allocateVector(&U, (imax+2) * (jmax+2))) {
		free(particles);
		free(FLAG);
		printf("Konnte keinen Speicherplatz fuer U allokieren");
		return 1;
	}
	if (allocateVector(&V, (imax+2) * (jmax+2))) {
		free(particles);
		free(FLAG);
		free(U);
		printf("Konnte keinen Speicherplatz fuer V allokieren");
		return 1;
	}
	if (allocateVector(&P, (imax+2) * (jmax+2))) {
		free(particles);
		free(FLAG);
		free(U);
		free(V);
		printf("Konnte keinen Speicherplatz fuer P allokieren");
		return 1;
	}
	if (allocateVector(&PSI, (imax+2) * (jmax+2))) {
		free(particles);
		free(FLAG);
		free(U);
		free(V);
		free(P);
		printf("Konnte keinen Speicherplatz fuer PSI allokieren");
		return 1;
	}	
	if (allocateVector(&ZETA, (imax+2) * (jmax+2))) {
		free(particles);
		free(FLAG);
		free(U);
		free(V);
		free(P);
		free(PSI);
		printf("Konnte keinen Speicherplatz fuer ZETA allokieren");
		return 1;
	}
	if (allocateVector(&TEMP, (imax+2) * (jmax+2))) {
		free(particles);
		free(FLAG);
		free(U);
		free(V);
		free(P);
		free(PSI);
		free(ZETA);
		printf("Konnte keinen Speicherplatz fuer TEMP allokieren");
		return 1;
	}
	if (allocateVector(&HEAT, (imax+2) * (jmax+2))) {
		free(particles);
		free(FLAG);
		free(U);
		free(V);
		free(P);
		free(PSI);
		free(ZETA);
		free(TEMP);
		printf("Konnte keinen Speicherplatz fuer HEAT allokieren");
		return 1;
	}
	
	calculateFluidDynamics(U, V, P, TEMP, FLAG, PSI, ZETA, HEAT, particles, partCount);
	
	free(particles);
	free(FLAG);
	free(U);
	free(V);
	free(P);
	free(PSI);
	free(ZETA);
	free(TEMP);
	free(HEAT);
	return 0;
}
