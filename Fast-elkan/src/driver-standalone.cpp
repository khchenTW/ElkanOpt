#include <iostream>
#include <fstream>
#include <string>
#include <cassert>

#include "general_functions.h"
#include "kmeans.h"
#include <chrono>

#include "dataset.h"
#include "elkan_kmeans.h"
#include "general_functions.h"
#include "FB1_elkan_kmeans.h"
#include "MO_elkan_kmeans.h"
#include "FB2_elkan_kmeans.h"
#include "yinyang_kmeans.h"
#define N 145751

#define D 74

#define STATIC 0

#define multialg 0


//3D_spatial_network_tab2.txt
double datasetGlobal[N*D];
double *backupX;

using namespace std::chrono;
//std::chrono::nanoseconds total_time(0);
Dataset *load_dataset(std::string const &filename) {
    std::ifstream input(filename.c_str());

    int n, d;
    //std::cout << n << "\n";
    //std::cout << d << "\n";
    input >> n >> d;

    Dataset *x = new Dataset(n, d);

    double *copyP1 = x->data;
    backupX = copyP1;
    //std::cout << "Dynamic" << std::endl;
    for (int i = 0; i < n * d; ++i) {
        input >> x->data[i];
        //std::cout << copyP1 << "\n";
        copyP1++;
    }

    #if STATIC
    for (int i = 0; i < N*D; ++i) {
        datasetGlobal[i] = x->data[i];
    }


    x->data = datasetGlobal;
    #endif

    double *copyP = datasetGlobal;
    //std::cout << "Global" << std::endl;
    for (int i = 0; i < N*D; ++i) {
        //std::cout << copyP << "\n";
        copyP++;
    }
    //std::cout << &(x->data[i]) << "\n";

    //std::cout << &(x) << "\n";
    return x;
}

Kmeans *get_algorithm(std::string const &name) {
    if (name == "elkan") return new ElkanKmeans();
    if (name == "FB_elkan") return new FB1_ElkanKmeans();
    if (name == "FB1_elkan") return new FB1_ElkanKmeans();
    if (name == "FB2_elkan") return new FB2_ElkanKmeans();
    if (name == "MO_elkan") return new MO_ElkanKmeans();
    if (name == "yinyang") return new yinyang_Kmeans();
    assert(false);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc != 5) {
        std::cout << "usage: " << argv[0] << " algorithm dataset k [centers|assignment]\n";
        return 1;
    }

    std::string algorithm_name(argv[1]);
    std::string filename(argv[2]);
    int k = std::stoi(argv[3]);
    std::string output(argv[4]);

    Dataset *x = load_dataset(filename);

    //std::cout << x << "\n";
    if (algorithm_name=="FB_elkan" && k >=300){
        algorithm_name="FB2_elkan";
    }
    Kmeans *algorithm = get_algorithm(algorithm_name);
#if multialg
    //Kmeans *algorithm2 = get_algorithm("elkan_new2");
    Kmeans *algorithm3 = get_algorithm("elkan_new3");
#endif
    //Dataset *initialCenters = init_centers_kmeanspp_v2(*x, k);
    Dataset *initialCenters = init_centers(*x, k);

    unsigned short *assignment = new unsigned short[x->n];

    assign(*x, *initialCenters, assignment);
    //auto start_time = std::chrono::high_resolution_clock::now();

    algorithm->initialize(x, k, assignment, 1);
    //x->print();
    auto start = std::chrono::system_clock::now();
    algorithm->run(10000);
    //x->print();
    //total_time += (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start_time));
    //std::cout <<total_time.count()<< "\n";
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout <<elapsed_seconds.count()<< "\n";

#if multialg
    //<elkan_new2
    x = load_dataset(filename);
    assignment = new unsigned short[x->n];
    assign(*x, *initialCenters, assignment);


    algorithm->initialize(x, k, assignment, 1);

    start = std::chrono::system_clock::now();
    algorithm->run(10000);
    //total_time += (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start_time));
    //std::cout <<total_time.count()<< "\n";
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout <<elapsed_seconds.count()<< "s\n";
    //elkan_new3
    x = load_dataset(filename);
    assignment = new unsigned short[x->n];
    assign(*x, *initialCenters, assignment);


    algorithm3->initialize(x, k, assignment, 1);

    start = std::chrono::system_clock::now();
    algorithm3->run(10000);
    //total_time += (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start_time));
    //std::cout <<total_time.count()<< "\n";
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout <<elapsed_seconds.count()<< "s\n";
#endif
/*
    Dataset const *finalCenters = algorithm->getCenters();

    if (output == "centers") {

        finalCenters->print();
    } else {
        assign(*x, *finalCenters, assignment);
        for (int i = 0; i < x->n; ++i) {
            std::cout << assignment[i] << "\n";
        }
    }

*/
    //if (algorithm_name != "elkan") {
    x->data = backupX;

    delete x;
    delete initialCenters;
    //}
    delete algorithm;

    delete [] assignment;

    return 0;
}
