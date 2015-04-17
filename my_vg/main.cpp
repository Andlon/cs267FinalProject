#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "points.h"

using namespace std;

int get_num_pairs(char* file);
void read_data(char*, point* data);
double calc_dist(point* pt1, point* pt2);
void register_pairs(int numPoints, point* points, vector<ptPair>* pairs);
double* find_dist_range(vector<ptPair>* pairs);
int get_bin(double dist, double* range, int numBins);
//void fill_variogram(vector<ptPair> pairs, double* gamma);
void print_gamma(double* values, int* counts, double* dists, int numBins, char* filename);


int main(int argc, char* argv[])
{
    //cout << "Hello world!" << endl;
    char *inFile = argv[1];
    //1. Read in measurements (value and coordinates)
    int numPoints = get_num_pairs(inFile);
    //cout << "Number of points: " << numPoints << endl;
    point *data = new point[numPoints];  //use vector?
    read_data(inFile, data);
    //for(int i = 0; i < numPoints; i++)
    //    cout << data[i].value << " is at " << data[i].x << "," << data[i].y << endl;
    //2. Determine all the pairs of points
    vector<ptPair> allPairs;
    register_pairs(numPoints, data, &allPairs);
    ////cout << "First pair is this far apart: " << allPairs[0].dist << "\n";
    //3. For each pair of points
    double *range = find_dist_range(&allPairs);
    cout << "Range: " << range[0] << "," << range[1] << "\n";
    int numBins = 15;
    double *gamma = new double[numBins];
    int *numPairs = new int[numBins];
    double *avgDists = new double[numBins];
    for(int i = 0; i < numBins; i++){
        gamma[i] = 0;
        numPairs[i] = 0;
        avgDists[i] = 0;
    }
       // - Calculate squared difference of values
       // - Calculate distance from coordinates
       // - Add the squared difference to the appropriate distance bin
    for(int i = 0 ; i < allPairs.size(); i++)
    {
         if(i==0) cout << "First pair values: " << (allPairs[i].pt2)->value << " and " <<
            (allPairs[i].pt1)->value << "\n";
         double sqdiff = pow((allPairs[i].pt2)->value - (allPairs[i].pt1)->value, 2);
         if(i==0) cout << "First sqdiff: " << sqdiff << "\n";
         int bin = get_bin(allPairs[i].dist, range, numBins);
         if(i==0) cout << "First bin: " << bin << "\n";
         gamma[bin] += sqdiff;
         avgDists[bin] += allPairs[i].dist;
         numPairs[bin]++;
    }
    //4. Divide the sum of squared differences by the number of pairs in bin
    for(int i = 0; i < numBins; i++)
    {
        gamma[i] /= (2*numPairs[i]);
        avgDists[i] /= numPairs[i];
    }
    //5. Print output
    char *outFile = "output.txt";
    print_gamma(gamma, numPairs, avgDists, numBins, outFile);

    delete[] data;
    delete[] gamma;
    delete[] numPairs;
    delete[] avgDists;
    return 0;
}

int get_num_pairs(char* file)
{
    int numPairs=0;
    ifstream f(file);
    string line;
    for (int i = 0; getline(f, line); ++i)
    numPairs++;
    return numPairs;  //Assumes no header lines
}

void read_data(char* file, point* data)
{
    ifstream f(file);
    string line;
    for (int i = 0; getline(f, line); ++i)
    {
        istringstream is(line);
        is >> data[i].value >> data[i].x >>  data[i].y ;
    }
}

//double calc_dist(double x1, double x2, double y1, double y2)
double calc_dist(point* pt1, point* pt2)
{
    double x_comp = pow(abs(pt2->x - pt1->x), 2);
    double y_comp = pow(abs(pt2->y - pt1->y), 2);
    return sqrt(x_comp+y_comp);
}

void register_pairs(int numPoints, point* points, vector<ptPair>* pairs)
{
    for(int i = 0; i < numPoints; i++)
    {
        for(int j = i+1; j < numPoints; j++)
        {
            ptPair newPair;
            newPair.pt1 = &points[i];
            newPair.pt2 = &points[j];
            newPair.dist = calc_dist(&points[i],&points[j]);
            pairs->push_back(newPair);
        }
    }
}

double* find_dist_range(vector<ptPair>* pairs)
{
    double *minnmax = new double[2];
    double *dists = new double[pairs->size()];
    for(int i = 0 ; i < pairs->size(); i++)
    {
        dists[i] = (*pairs)[i].dist;
    }
    //cout << "First distance: " << dists[0] << "\n";
    minnmax[0] = *min_element(dists,dists+pairs->size());
    minnmax[1] = *max_element(dists,dists+pairs->size());
    return minnmax;
}

int get_bin(double dist, double* range, int numBins)
{
    double interval = (range[1] - range[0])/numBins;
    return floor(dist/interval);
}
//void fill_variogram(vector<ptPair> pairs, double* gamma, )
//{

//}

void print_gamma(double* values, int* counts, double* dists, int numBins, char* filename)
{
    ofstream ofile;
    ofile.open (filename);
    ofile << "# \tfrom \tto  \tn_pairs \tav_dist \tsemivariogram \n";
    for(int i = 0; i < numBins; i++)
    {
        ofile << "?" << "\t" << "?" << "\t" << counts[i] << "\t" <<
            dists[i] << "\t" << values[i] << "\n";
    }
    ofile.close();
}
