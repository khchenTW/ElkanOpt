/* Authors: Greg Hamerly and Jonathan Drake
 * Feedback: hamerly@cs.baylor.edu
 * See: http://cs.baylor.edu/~hamerly/software/kmeans.php
 * Copyright 2014
 */

#include "MO_elkan_kmeans.h"
#include "general_functions.h"
#include <cmath>
#include <chrono>
//using namespace std::chrono;

#define Time 0
#define Countdistance 0
void MO_ElkanKmeans::update_center_dists(int threadId) {
    // find the inter-center distances
    for (int c1 = 0; c1 < k; ++c1) {
        if (c1 % numThreads == threadId) {
            s[c1] = std::numeric_limits<double>::max();

            for (int c2 = 0; c2 < k; ++c2) {
                // we do not need to consider the case when c1 == c2 as centerCenterDistDiv2[c1*k+c1]
                // is equal to zero from initialization, also this distance should not be used for s[c1]
                if (c1 != c2) {
                    // divide by 2 here since we always use the inter-center
                    // distances divided by 2
                    //std::cout <<sqrt(centerCenterDist2(c1, c2))<< "\n";
                    oldcenterCenterDistDiv2[c1 * k + c2] = centerCenterDistDiv2[c1 * k + c2];
                    centerCenterDistDiv2[c1 * k + c2] = sqrt(centerCenterDist2(c1, c2)) / 2.0;

                    if (centerCenterDistDiv2[c1 * k + c2] < s[c1]) {
                        s[c1] = centerCenterDistDiv2[c1 * k + c2];
                    }
                }
            }
        }
    }
}

int MO_ElkanKmeans::runThread(int threadId, int maxIterations) {
    int iterations = 0;

    int startNdx = start(threadId);
    int endNdx = end(threadId);
    //x->print();
    //auto start_time = std::chrono::high_resolution_clock::now();

    ub_old = new double[n];
    std::fill(ub_old, ub_old + n, std::numeric_limits<double>::max());

    oldcenterCenterDistDiv2 = new double[k * k];
    std::fill(oldcenterCenterDistDiv2, oldcenterCenterDistDiv2 + k * k, 0.0);

    oldcenter2newcenterDis = new double[k * k];
    std::fill(oldcenter2newcenterDis, oldcenter2newcenterDis + k * k, 0.0);
    oldcenters = new double[k * d];
    //oldcenters->fill(0.0);
    std::fill(oldcenters, oldcenters + k * d, 0.0);
#if Time
    auto start_time = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> total_elkan_time{};
    std::chrono::duration<double> elapsed_seconds = end-start;
#endif
    while ((iterations < maxIterations) && ! converged) {
#if Time
        start_time = std::chrono::high_resolution_clock::now();
        start = std::chrono::system_clock::now();
#endif
        ++iterations;
        #if Countdistance
        int numberdistances =0;
        #endif
        update_center_dists(threadId);
        synchronizeAllThreads();

        for (int i = startNdx; i < endNdx; ++i) {
            //std::cout << d << "\n";
            unsigned short closest = assignment[i];
            bool r = true;

            if (upper[i] <= s[closest]) {
                continue;
            }

            for (int j = 0; j < k; ++j) {
                if (j == closest) { continue; }
                if (upper[i] <= 2.0*(oldcenterCenterDistDiv2[assignment[i] * k + j]) - ub_old[i]- centerMovement[j]) { continue; }
                if ( upper[i] <= oldcenter2newcenterDis[assignment[i] * k + j] - ub_old[i]) { continue; }  //upper[i] <= lower[i * k + j] ||
                if (upper[i] <= centerCenterDistDiv2[closest * k + j]) { continue; }
                #if Countdistance
                numberdistances ++;
                #endif
                // ELKAN 3(a)
                if (r) {
                    upper[i] = sqrt(pointCenterDist2(i, closest));
                    //lower[i * k + closest] = upper[i];
                    //lower2[i * k + closest] = upper[i];
                    r = false;
                    //if ((upper[i] <= lower[i * k + j]) || (upper[i] <= centerCenterDistDiv2[closest * k + j])) {
                        //continue;
                    //}
                }

                // ELKAN 3(b)
                //lower[i * k + j] = sqrt(pointCenterDist2(i, j));

                if (sqrt(pointCenterDist2(i, j)) < upper[i]) {
                    closest = j;
                    upper[i] = sqrt(pointCenterDist2(i, j));
                }
            }
            if (assignment[i] != closest) {
                changeAssignment(i, closest, threadId);
            }


        }

        #if Countdistance
        std::cout <<numberdistances<< "\n";
        #endif

#if Time
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout <<elapsed_seconds.count()<< "\n";
    total_elkan_time += (std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - start_time));
#endif


        //verifyAssignment(iterations, startNdx, endNdx);

        // ELKAN 4, 5, AND 6
        synchronizeAllThreads();

        if (threadId == 0) {
            int furthestMovingCenter = move_centers_newbound(oldcenters,oldcenter2newcenterDis);
            converged = (0.0 == centerMovement[furthestMovingCenter]);
        }

        synchronizeAllThreads();
        //total_elkan_time += (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start_time));
        if (! converged) {
            update_bounds(startNdx, endNdx);
        }
        else{
            //std::cout << iterations << "\n";
            #if Time
            std::cout <<total_elkan_time.count()<< "\n";
            #endif
        }
        synchronizeAllThreads();

    }

    return iterations;
}

void MO_ElkanKmeans::update_bounds(int startNdx, int endNdx) {
    //lower3
    /*
    for (int i = startNdx; i < endNdx; ++i) {

        for (int j = 0; j < k; ++j) {

            lower3[i * numLowerBounds + j] =  2.0*(centerCenterDistDiv2[assignment[i] * k + j]) - upper[i]- centerMovement[j];

        }
    }
    */
    for (int i = startNdx; i < endNdx; ++i) {
        ub_old[i] = upper[i];
        upper[i] += centerMovement[assignment[i]];

    }
    /*
    for (int i = startNdx; i < endNdx; ++i) {

        for (int j = 0; j < k; ++j) {

            lower2[i * numLowerBounds + j] = oldcenter2newcenterDis[assignment[i] * k + j] - upper[i];

        }
    }
    */
    for (int i = startNdx; i < endNdx; ++i) {
        upper[i] += centerMovement[assignment[i]];
        //for (int j = 0; j < k; ++j) {
            //lower[i * numLowerBounds + j] -= centerMovement[j];

        //}
    }
}

void MO_ElkanKmeans::initialize(Dataset const *aX, unsigned short aK, unsigned short *initialAssignment, int aNumThreads) {
    numLowerBounds = aK;
    TriangleInequalityBaseKmeans::initialize(aX, aK, initialAssignment, aNumThreads);
    centerCenterDistDiv2 = new double[k * k];
    std::fill(centerCenterDistDiv2, centerCenterDistDiv2 + k * k, 0.0);

}

void MO_ElkanKmeans::free() {
    TriangleInequalityBaseKmeans::free();
    delete [] centerCenterDistDiv2;
    centerCenterDistDiv2 = NULL;
    //delete [] oldcenterCenterDistDiv2;
    //oldcenterCenterDistDiv2 = NULL;
    delete centers;
    centers = NULL;
}
int MO_ElkanKmeans::move_centers_newbound(double *oldcenters, double *oldcenter2newcenterDis) {
    int furthestMovingCenter = 0;
    /*
    for (int j = 0; j < k; ++j) {

        //std::cout << oldcenters[ j]<< "\n";
        for (int dim = 0; dim < d; ++dim) {
            std::cout << oldcenters[j] << "\n";
        }
    }
    */
    for (int j = 0; j < k; ++j) {
        centerMovement[j] = 0.0;
        int totalClusterSize = 0;
        double old = 0;
        for (int t = 0; t < numThreads; ++t) {
            totalClusterSize += clusterSize[t][j];
        }
        if (totalClusterSize > 0) {
            for (int dim = 0; dim < d; ++dim) {
                double z = 0.0;
                for (int t = 0; t < numThreads; ++t) {
                    z += (*sumNewCenters[t])(j,dim);
                }
                z /= totalClusterSize;
                //std::cout << z << "\n";
                //std::cout << (*centers)(j, dim) << "\n";
                centerMovement[j] += (z - (*centers)(j, dim)) * (z - (*centers)(j, dim));//calculate distance
                //std::cout << (*centers)(j, dim) << "\n";
                old = (*centers)(j, dim);
                //std::cout << (*oldcenters)(j, dim) << "\n";
                oldcenters[j* d+ dim] = old;
                //std::cout << (*centers)(j, dim) << "\n";
                (*centers)(j, dim) = z; //update new centers
            }
        }
        centerMovement[j] = sqrt(centerMovement[j]);

        if (centerMovement[furthestMovingCenter] < centerMovement[j]) {
            furthestMovingCenter = j;
        }
    }

    for (int c1 = 0; c1 < k; ++c1) {

        for (int c2 = 0; c2 < k; ++c2)
            if (c1 != c2) {
                oldcenter2newcenterDis[c1 * k + c2]=0.0;
                for (int dim = 0; dim < d; ++dim) {
                    oldcenter2newcenterDis[c1 * k + c2] += (oldcenters[c1* d+ dim] - (*centers)(c2, dim)) * (oldcenters[c1* d+ dim] - (*centers)(c2, dim));
                }
                oldcenter2newcenterDis[c1 * k + c2] = sqrt(oldcenter2newcenterDis[c1 * k + c2]);
            }
    }

    #ifdef COUNT_DISTANCES
    numDistances += k;
    #endif

    return furthestMovingCenter;
}
/*
int ElkanKmeans_newbound_opt::move_centers2(int *sortindex,bool sorting, bool *nochanged) {
    int furthestMovingCenter = 0;
    for (int j = 0; j < k; ++j) {
        if (sorting == true){
        //j=sortindex[j];
            auto itr = std::find(sortindex,sortindex + k,j);


            centerMovement[j] = 0.0;
            int totalClusterSize = 0;
            for (int t = 0; t < numThreads; ++t) {
                totalClusterSize += clusterSize[t][j];
            }
            if (totalClusterSize > 0) {
                if (nochanged[j]==false){


                for (int dim = 0; dim < d; ++dim) {
                    double z = 0.0;
                    for (int t = 0; t < numThreads; ++t) {
                        z += (*sumNewCenters[t])(j,dim);
                    }
                    z /= totalClusterSize;
                    //std::cout << z << "\n";

                    centerMovement[j] += (z - (*centers)(std::distance(sortindex,itr), dim)) * (z - (*centers)(std::distance(sortindex,itr), dim));//calculate distance
                    (*centers)(std::distance(sortindex,itr), dim) = z; //update new centers
                }

            centerMovement[j] = sqrt(centerMovement[j]);
                }
            if (centerMovement[furthestMovingCenter] < centerMovement[j]) {
                furthestMovingCenter = j;
            }
            }
        }
        else{
        centerMovement[j] = 0.0;
        int totalClusterSize = 0;
        for (int t = 0; t < numThreads; ++t) {
            totalClusterSize += clusterSize[t][j];
        }
        if (totalClusterSize > 0) {
            if (nochanged[j]==false){
            for (int dim = 0; dim < d; ++dim) {
                double z = 0.0;
                for (int t = 0; t < numThreads; ++t) {
                    z += (*sumNewCenters[t])(j,dim);
                }
                z /= totalClusterSize;
                //std::cout << z << "\n";
                centerMovement[j] += (z - (*centers)(j, dim)) * (z - (*centers)(j, dim));//calculate distance
                (*centers)(j, dim) = z; //update new centers
            }

        centerMovement[j] = sqrt(centerMovement[j]);
            }
        if (centerMovement[furthestMovingCenter] < centerMovement[j]) {
            furthestMovingCenter = j;
        }
        }

        }
    }

    #ifdef COUNT_DISTANCES
    numDistances += k;
    #endif

    return furthestMovingCenter;
}
*/
