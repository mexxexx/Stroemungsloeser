#ifndef BOUNDARY_H
#define BOUNDARY_H

void applyHomogenousNeumannBC(double *p, int imax, int jmax);

void setBoundaryCond(double *U, double *V, double *TEMP, char *FLAG, double delx, double dely, int imax, int jmax, int wl, int wr, int wt, int wb,
					int tl, double tl_value, int tr, double tr_value, int tt, double tt_value, int tb, double tb_value);

void setSpecialBoundaryCond(double *U, double *V, double *TEMP, int imax, int jmax, char *problem);

void initDrivenCavity(double *U, double *V, int imax, int jmax);

void initKarman(double *U, double *V, int imax, int jmax);

void obstacleBC(double *p, char *FLAG, int imax, int jmax, double deltaX, double deltaY);

#endif
