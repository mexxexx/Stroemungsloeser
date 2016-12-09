#ifndef PARTICLES_H
#define PARTICLES_H

typedef struct Partikel{

	double x, y;
	double u, v;
	char isActive;


}particle; 


void ParticleSeed(Partikel *particles, double posx1, double posy1, double posx2, double posy2, int partcount, int anzahl);



#endif
