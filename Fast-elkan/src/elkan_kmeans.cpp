/* Authors: Greg Hamerly and Jonathan Drake
 * Feedback: hamerly@cs.baylor.edu
 * See: http://cs.baylor.edu/~hamerly/software/kmeans.php
 * Copyright 2014
 */

#include "elkan_kmeans.h"
#include "general_functions.h"
#include <cmath>
#include <chrono>
//using namespace std::chrono;

#define Time 0
#define Countdistance 0
void ElkanKmeans::update_center_dists(int threadId) {
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
                    centerCenterDistDiv2[c1 * k + c2] = sqrt(centerCenterDist2(c1, c2)) / 2.0;

                    if (centerCenterDistDiv2[c1 * k + c2] < s[c1]) {
                        s[c1] = centerCenterDistDiv2[c1 * k + c2];
                    }
                }
            }
        }
    }
}

int ElkanKmeans::runThread(int threadId, int maxIterations) {
    int iterations = 0;

    int startNdx = start(threadId);
    int endNdx = end(threadId);
    lower = new double[n * k];
    std::fill(lower, lower + n * k, 0.0);
    //x->print();
    //auto start_time = std::chrono::high_resolution_clock::now();
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
                if (upper[i] <= lower[i * k + j]) { continue; }
                if (upper[i] <= centerCenterDistDiv2[closest * k + j]) { continue; }
                #if Countdistance
                numberdistances ++;
                #endif
                // ELKAN 3(a)
                if (r) {
                    upper[i] = sqrt(pointCenterDist2(i, closest));
                    lower[i * k + closest] = upper[i];
                    r = false;
                    if ((upper[i] <= lower[i * k + j]) || (upper[i] <= centerCenterDistDiv2[closest * k + j])) {
                        continue;
                    }
                }

                // ELKAN 3(b)
                lower[i * k + j] = sqrt(pointCenterDist2(i, j));
                if (lower[i * k + j] < upper[i]) {
                    closest = j;
                    upper[i] = lower[i * k + j];
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
            int furthestMovingCenter = move_centers();
            converged = (0.0 == centerMovement[furthestMovingCenter]);
        }

        synchronizeAllThreads();
        //total_elkan_time += (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start_time));
        if (! converged) {
            update_bounds(startNdx, endNdx);
        }
        else{
            std::cout << iterations << "\n";
            #if Time
            std::cout <<total_elkan_time.count()<< "\n";
            #endif
        }
        synchronizeAllThreads();

    }

    return iterations;
}

void ElkanKmeans::update_bounds(int startNdx, int endNdx) {
    for (int i = startNdx; i < endNdx; ++i) {
        upper[i] += centerMovement[assignment[i]];
        for (int j = 0; j < k; ++j) {
            lower[i * numLowerBounds + j] -= centerMovement[j];
        }
    }
}

void ElkanKmeans::initialize(Dataset const *aX, unsigned short aK, unsigned short *initialAssignment, int aNumThreads) {
    numLowerBounds = aK;
    TriangleInequalityBaseKmeans::initialize(aX, aK, initialAssignment, aNumThreads);
    centerCenterDistDiv2 = new double[k * k];
    std::fill(centerCenterDistDiv2, centerCenterDistDiv2 + k * k, 0.0);
}

void ElkanKmeans::free() {
    TriangleInequalityBaseKmeans::free();
    delete [] centerCenterDistDiv2;
    centerCenterDistDiv2 = NULL;
    //delete centers;
    //centers = NULL;
}

