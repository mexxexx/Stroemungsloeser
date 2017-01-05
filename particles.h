#ifndef PARTICLES_H
#define PARTICLES_H

typedef struct Particle{
	double x, y;
	double u, v;
	char isActive;
}Particle; 

void particleSeed(Particle *particles, double posx1, double posy1, double posx2, double posy2, int partCount, int anzahl);
void particleVelocity(double *U, double *V, double delx, double dely, int imax, int jmax, Particle *particles, int partCount);
void particleTransport(Particle *particles, double delt, int partCount, double xlength, double ylength);
void particleInit(Particle *particles, int partCount);
void printParticles(Particle *particles, int partCount, char *filename);
void bmp_verify(char* FLAG, int imax, int jmax);

#endif
