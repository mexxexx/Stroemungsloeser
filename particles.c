#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void ParticleSeed(particle *particles, double posx1, double posy1, double posx2, double posy2, int partcount, int anzahl){


	int count=0;
	double horizontal=(1./anzahl)*fabs(posx1-posx2);
	double vertical=(1./anzahl)*fabs(posy1-posy2);
	for (int i=0;i<partcount;i++) {
		
		if(count>anzahl) break;

		if(particles[i].isActive) continue;

		else{
			particles[i].x=count*horizontal;
			particles[i].y=count*vertical;
			particles[i].isActive=1;	
			count++;	
		}

	}





}
