void particleVelocity(double *U, double *V, double delx, double dely, int imax, int jmax, Particle *particles, int partCount) {
	double oneOverDelxDely = 1.0 / (delx * dely);
	for (int k = 0; k < partCount; k++) {
		double x = particles[k].x;
		double y = particles[k].y;
		
		int i = (int)(x/delx)+1;
		int j = (int)(y+0.5*dely)/dely)+1;
		
		double x1 = (i-1)*delx;
		double x2 = i*delx;
		double y1 = (j-1.5)*delx;
		double y2 = (j-0.5)*delx;
		
		double u1 = U[POS2D(i-1, j-1, imax+2)];
		double u2 = U[POS2D(i, j-1, imax+2)];
		double u3 = U[POS2D(i-1, j, imax+2)];
		double u4 = U[POS2D(i, j, imax+2)];
		
		particles[k].u = oneOverDelxDely * ((x2-x)*(y2-y)*u1+(x-x1)*(y2-y)*u2+(x2-x)*(y-y1)*u3+(x-x1)*(y-y1)*u4);
		
		i=(int)((x+0.5*delx)/delx)+1;
		j=(int)(y/dely)+1;
		
		x1 = (i-1.5)*delx;
		x2 = (i-0.5)*delx;
		y1 = (j-1)*delx;
		y2 = j*delx;
		
		double v1 = V[POS2D(i-1, j-1, imax+2)];
		double v2 = V[POS2D(i, j-1, imax+2)];
		double v3 = V[POS2D(i-1, j, imax+2)];
		double v4 = V[POS2D(i, j, imax+2)];
		
		particles[k].v = oneOverDelxDely * ((x2-x)*(y2-y)*v1+(x-x1)*(y2-y)*v2+(x2-x)*(y-y1)*v3+(x-x1)*(y-y1)*v4);
	}
}

void particleTransport(Particle *particles, double delt, int partCount, double xlength, double ylength) {
	for (int k = 0; k < partCount; k++) {
		Particle *p = &particles[k];
		if (!p->isActive) continue;
		
		p->x += delt * p->u;
		p->y += delt * p->v;
		
		if (p->x < 0 || p->x >= xlength || p->y < 0 || p->y >= ylength)
			p->isActive = 0;
	}
}
