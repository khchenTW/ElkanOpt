/* Authors: Greg Hamerly and Jonathan Drake
 * Feedback: hamerly@cs.baylor.edu
 * See: http://cs.baylor.edu/~hamerly/software/kmeans.php
 * Copyright 2014
 */

#include "yinyang_kmeans.h"
#include "general_functions.h"
#include <cmath>
#include <chrono>
#include <vector>
//using namespace std::chrono;

#define Time 0
#define Countdistance 0
void yinyang_Kmeans::update_center_dists(int threadId) {
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
                    centerCenterDistDiv2[c1 * k + c2] = sqrt(centerCenterDist2(c1, c2)) / 2.0;

                    if (centerCenterDistDiv2[c1 * k + c2] < s[c1]) {
                        s[c1] = centerCenterDistDiv2[c1 * k + c2];
                    }
                }
            }
        }
    }
}
int yinyang_Kmeans::runThread(int threadId, int maxIterations) {
    int iterations = 0;

    int startNdx = start(threadId);
    int endNdx = end(threadId);
    ub_old = new double[n];
    std::fill(ub_old, ub_old + n, std::numeric_limits<double>::max());
    lower = new double[n * k];
    std::fill(lower, lower + n * k, 0.0);
    //lower2 = new double[n * k];
    //std::fill(lower2, lower2 + n * k, 0.0);
    oldcenter2newcenterDis = new double[k * k];
    std::fill(oldcenter2newcenterDis, oldcenter2newcenterDis + k * k, 0.0);
    oldcenters = new double[k * d];
    //oldcenters->fill(0.0);
    std::fill(oldcenters, oldcenters + k * d, 0.0);
    upperc = new double[k];
    std::fill(upperc, upperc + k, 0.0);
    std::vector<std::vector<int>> clusterp;
    std::vector<std::vector<int>> oldclusterp;
    std::vector<std::vector<int>> clusterassignment;
    int g = 0;
    if ((k-1)/10 == 0) {
        g = k;
    }
    else{
        g=((k-1)/10)+1;
    }

    std::cout <<g<< "\n";
    grouplower = new double[n * g];
    std::fill(lower, lower + n * g, 0.0);
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
        if (iterations ==1){
            for (int c1 =0; c1 < k; ++c1){
                clusterp.push_back(std::vector<int>());
            }
            for (int c1 =0; c1 < k; ++c1){
                clusterassignment.push_back(std::vector<int>());
            }
            for (int i = startNdx; i < endNdx; ++i) {
                //std::cout << d << "\n";
                unsigned short closest = assignment[i];
                bool r = true;

                if (upper[i] <= s[closest]) {
                    continue;
                }

                for (int j = 0; j < k; ++j) {
                    if (j == closest) { continue; }
                    //if (upper[i] <= lower3[i * k + j]) { continue; }
                    if (upper[i] <= lower[i * k + j] || upper[i] <= oldcenter2newcenterDis[assignment[i] * k + j] - ub_old[i]) { continue; }
                    if (upper[i] <= centerCenterDistDiv2[closest * k + j]) { continue; }
                    #if Countdistance
                    numberdistances ++;
                    #endif
                    // ELKAN 3(a)
                    if (r) {
                        upper[i] = sqrt(pointCenterDist2(i, closest));
                        lower[i * k + closest] = upper[i];
                        //lower2[i * k + closest] = upper[i];
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
                    if ((j+1)%10==0){
                        grouplower[i * g + j/10] = upper[i];
                        //groupminimum=
                    }
                }

                clusterp[closest].push_back(i);
                //clusterp[closest].push_back(oldclusterp[c1][i]);
                if (assignment[i] != closest) {
                    changeAssignment(i, closest, threadId);

                }
            }

            //std::cout <<clusterassignment.size()<< "\n";
            for (int c1 =0; c1 < k; ++c1){
                oldclusterp.push_back(std::vector<int>());
                for (int i= 0; i < clusterp[c1].size();++i){
                    oldclusterp[c1].push_back(clusterp[c1][i]);
                }
            }
        }
        else{
            for (int c1 = 0; c1 < k; ++c1) {
                for (int i = 0; i < oldclusterp[c1].size();++i){
                    if (upperc[c1]<upper[oldclusterp[c1][i]]){
                        upperc[c1] = upper[oldclusterp[c1][i]];
                    }
                }
            }
            for (int c1 = 0; c1 < k; ++c1) {
                clusterassignment[c1].clear();
                for (int c2 = 0; c2 < k; ++c2){
                    if (c1 == c2) { continue; }
                    if (upperc[c1]>=centerCenterDistDiv2[c1 * k + c2]){
                        clusterassignment[c1].push_back(c2);
                    }
                }
            }
            for (int c1 =0; c1 < k; ++c1){
                clusterp[c1].clear();
            }

            for (int c1 = 0; c1 < k; ++c1) {
                //std::cout <<oldclusterp[c1].size()<< "\n";
                if (clusterassignment[c1].size()!=0){
                    for (int i = 0; i < oldclusterp[c1].size();++i){
                        unsigned short closest = assignment[oldclusterp[c1][i]];
                        bool r = true;

                        if (upper[oldclusterp[c1][i]] <= s[closest]) {
                            clusterp[closest].push_back(oldclusterp[c1][i]);
                            continue;
                        }

                        for (int j = 0; j < clusterassignment[c1].size(); ++j) {
                            if (clusterassignment[c1][j] == closest) { continue; }
                            //if (upper[i] <= lower3[i * k + j]) { continue; }
                            if (upper[oldclusterp[c1][i]] <= lower[oldclusterp[c1][i] * k + clusterassignment[c1][j]] || upper[oldclusterp[c1][i]] <= oldcenter2newcenterDis[assignment[oldclusterp[c1][i]] * k + clusterassignment[c1][j]] - ub_old[oldclusterp[c1][i]]) { continue; }
                            if (upper[oldclusterp[c1][i]] <= centerCenterDistDiv2[closest * k + clusterassignment[c1][j]]) { continue; }
                            #if Countdistance
                            numberdistances ++;
                            #endif
                            // ELKAN 3(a)
                            if (r) {
                                upper[oldclusterp[c1][i]] = sqrt(pointCenterDist2(oldclusterp[c1][i], closest));
                                //lower[oldclusterp[c1][i] * k + closest] = upper[oldclusterp[c1][i]];
                                //lower2[i * k + closest] = upper[i];
                                r = false;
                                //if ((upper[oldclusterp[c1][i]] <= lower[oldclusterp[c1][i] * k + clusterassignment[c1][j]]) || (upper[oldclusterp[c1][i]] <= centerCenterDistDiv2[closest * k + clusterassignment[c1][j]])|| upper[oldclusterp[c1][i]] <= oldcenter2newcenterDis[assignment[oldclusterp[c1][i]] * k + clusterassignment[c1][j]] - ub_old[oldclusterp[c1][i]]) {
                                    //continue;
                                //}
                            }

                            // ELKAN 3(b)
                            lower[oldclusterp[c1][i] * k + clusterassignment[c1][j]] = sqrt(pointCenterDist2(oldclusterp[c1][i], clusterassignment[c1][j]));

                            if (lower[oldclusterp[c1][i] * k + clusterassignment[c1][j]] < upper[oldclusterp[c1][i]]) {
                                closest = clusterassignment[c1][j];
                                upper[oldclusterp[c1][i]] = lower[oldclusterp[c1][i] * k + clusterassignment[c1][j]];
                            }
                        }
                        clusterp[closest].push_back(oldclusterp[c1][i]);
                        //clusterp[closest].push_back(oldclusterp[c1][i]);
                        if (assignment[oldclusterp[c1][i]] != closest) {
                            changeAssignment(oldclusterp[c1][i], closest, threadId);

                        }
                    }
                }
                //else{
                    //std::cout <<iterations<< "\n";
                    //std::cout <<c1<< "\n";
                //}

            }

            //std::cout <<clusterassignment.size()<< "\n";
            for (int c1 =0; c1 < k; ++c1){
                oldclusterp[c1].clear();
                //std::cout <<clusterp[c1].size()<< "\n";
                for (int i= 0; i < clusterp[c1].size();++i){
                    oldclusterp[c1].push_back(clusterp[c1][i]);
                }
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
            std::cout << iterations << "\n";
            #if Time
            std::cout <<total_elkan_time.count()<< "\n";
            #endif
        }
        synchronizeAllThreads();

    }

    return iterations;
}

void yinyang_Kmeans::update_bounds(int startNdx, int endNdx) {
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
        for (int j = 0; j < k; ++j) {
            lower[i * numLowerBounds + j] -= centerMovement[j];
            //lower2[i * numLowerBounds + j] = oldcenter2newcenterDis[assignment[i] * k + j] - upper[i];
            //std::cout << lower[i * numLowerBounds + j] << "\n";
        }
    }
}

void yinyang_Kmeans::initialize(Dataset const *aX, unsigned short aK, unsigned short *initialAssignment, int aNumThreads) {
    numLowerBounds = aK;
    TriangleInequalityBaseKmeans::initialize(aX, aK, initialAssignment, aNumThreads);
    centerCenterDistDiv2 = new double[k * k];
    std::fill(centerCenterDistDiv2, centerCenterDistDiv2 + k * k, 0.0);

}

void yinyang_Kmeans::free() {
    TriangleInequalityBaseKmeans::free();
    delete [] centerCenterDistDiv2;
    centerCenterDistDiv2 = NULL;
    //delete [] oldcenterCenterDistDiv2;
    //oldcenterCenterDistDiv2 = NULL;
    delete centers;
    centers = NULL;
}
int yinyang_Kmeans::move_centers_newbound(double *oldcenters, double *oldcenter2newcenterDis) {
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
                //old = (*centers)(j, dim);
                //std::cout << (*oldcenters)(j, dim) << "\n";
                oldcenters[j* d+ dim] = (*centers)(j, dim);
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

