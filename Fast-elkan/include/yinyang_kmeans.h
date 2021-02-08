#ifndef YINYANG_KMEANS_H
#define YINYANG_KMEANS_H

/* Authors: Greg Hamerly and Jonathan Drake
 * Feedback: hamerly@cs.baylor.edu
 * See: http://cs.baylor.edu/~hamerly/software/kmeans.php
 * Copyright 2014
 *
 * Elkan's k-means algorithm that uses k lower bounds per point to prune
 * distance calculations.
 */


#include "triangle_inequality_base_kmeans.h"

class yinyang_Kmeans : public TriangleInequalityBaseKmeans {
    public:
        yinyang_Kmeans() : centerCenterDistDiv2(NULL) {}
        //ElkanKmeans_newbound() : oldcenterCenterDistDiv2(NULL) {}
        virtual ~yinyang_Kmeans() { free(); }
        virtual void free();
        virtual void initialize(Dataset const *aX, unsigned short aK, unsigned short *initialAssignment, int aNumThreads);
        virtual std::string getName() const { return "yinyang"; }

    protected:
        virtual int runThread(int threadId, int maxIterations);

        // Update the distances between each pair of centers.
        void update_center_dists(int threadId);

        // Update the upper and lower bounds for the range of points given.
        void update_bounds(int startNdx, int endNdx);

        // Keep track of the distance (divided by 2) between each pair of
        // points.
        double *centerCenterDistDiv2;
        double *oldcenter2newcenterDis;
        double *oldcenters;
        //double *lower2;
        double *lower;
        double *ub_old;
        double *upperc;
        double *grouplower;

        int move_centers2(int *sortindex,bool sorting,bool *nochanged);
        int move_centers_newbound(double *oldcenters,double *oldcenter2newcenterDis);
};

#endif

