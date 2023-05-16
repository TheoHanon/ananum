#ifndef DESIGN_H // le header guard
#define DESIGN_H 

void designTuningFork(double r1, double r2, double e, double l, double meshSizeFactor, char * filename);

void designTuningForkSym(double r1, double r2, double e, double l, double meshSizeFactor, char * filename);

void designTuningForkSym2(double r1, double r2, double e, double l, double c, double meshSizeFactor, char * filename);

void designTuningFork2(double r1, double r2, double e, double l, double c, double meshSizeFactor, char * filename);

#endif