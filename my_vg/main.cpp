#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string.h>

using namespace std;

struct point {
  double value;
  double x,y;
} ;

struct ptPair {
  point *pt1, *pt2;
  double dist;
} ;

int get_num_pairs(string file);
void read_data(string file, vector<point>* data);
double calc_dist(point* pt1, point* pt2);
void register_pairs(int numPoints, vector<point>* points, vector<ptPair>* pairs);
double* find_dist_range(vector<ptPair>* pairs);
int get_bin(double dist, double* range, int numBins);
//void fill_variogram(vector<ptPair> pairs, double* gamma);
void print_gamma(vector<double>* values, vector<int>* counts, vector<double>* dists, int numBins, string filename);


int main(int argc, char* argv[])
{
    string inFile;
    inFile = argv[1];
    //1. Read in measurements (value and coordinates)
    int numPoints = get_num_pairs(inFile);
    //point *data = new point[numPoints];  //use vector?
    vector<point> data (numPoints);
    read_data(inFile, &data);
    //2. Determine all the pairs of points
    vector<ptPair> allPairs;
    register_pairs(numPoints, &data, &allPairs);
    //3. For each pair of points
    double *range = find_dist_range(&allPairs);
    int numBins = 15;
    //double *gamma = new double[numBins];
    //int *numPairs = new int[numBins];
    //double *avgDists = new double[numBins];
    vector<double> gamma (numBins,0);
    vector<int> numPairs (numBins,0);
    vector<double> avgDists (numBins,0);
    /*
    for(int i = 0; i < numBins; i++){
        gamma[i] = 0;
        numPairs[i] = 0;
        avgDists[i] = 0;
    }
    */
    // - Calculate squared difference of values
    // - Calculate distance from coordinates
    // - Add the squared difference to the appropriate distance bin
    for(int i = 0 ; i < allPairs.size(); i++)
    {
         double sqdiff = pow((allPairs[i].pt2)->value - (allPairs[i].pt1)->value, 2);
         int bin = get_bin(allPairs[i].dist, range, numBins);
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
    //char *outFile = new char[strlen("output.txt") + 1];
    string outFile = "output.txt";
    //strcpy(outFile, "output.txt");
    print_gamma(&gamma, &numPairs, &avgDists, numBins, outFile);
  //cout << "Made it here!\n";
    /*
    delete[] data;
    delete[] gamma;
    delete[] numPairs;
    delete[] avgDists;
    delete[] outFile; 
    */    
    return 0;
}

int get_num_pairs(string file)
{
    int numPairs=0;
    ifstream f(file.c_str());
    string line;
    for (int i = 0; getline(f, line); ++i)
    numPairs++;
    return numPairs;  //Assumes no header lines
}

void read_data(string file,vector<point>* data)
{
    ifstream f(file.c_str());
    string line;
    for (int i = 0; getline(f, line); ++i)
    {
        istringstream is(line);
	point pt;
        is >> pt.value >> pt.x >>  pt.y ;
	(*data)[i] = pt;
    }
}

double calc_dist(point* pt1, point* pt2)
{
    double x_comp = pow(abs(pt2->x - pt1->x), 2);
    double y_comp = pow(abs(pt2->y - pt1->y), 2);
    return sqrt(x_comp+y_comp);
}

void register_pairs(int numPoints, vector<point>* points, vector<ptPair>* pairs)
{
    for(int i = 0; i < numPoints; i++)
    {
        for(int j = i+1; j < numPoints; j++)
        {
            ptPair newPair;
            newPair.pt1 = &(*points)[i];
            newPair.pt2 = &(*points)[j];
            newPair.dist = calc_dist(&(*points)[i],&(*points)[j]);
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
    minnmax[0] = *min_element(dists,dists+pairs->size());
    minnmax[1] = *max_element(dists,dists+pairs->size());
    delete[] dists;
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

void print_gamma(vector<double>* values, vector<int>* counts, vector<double>* dists, int numBins, string filename)
{
    ofstream ofile;
    ofile.open(filename.c_str());
    ofile << "#from \tto  \tn_pairs \tav_dist \tsemivariogram \n";
    for(int i = 0; i < numBins; i++)
    {
        ofile << "?" << "\t" << "?" << "\t" << (*counts)[i] << "\t" <<
            (*dists)[i] << "\t" << (*values)[i] << "\n";
    }
    cout << "About to close file\n";
    ofile.close();
    cout << "Closed file\n";
}
