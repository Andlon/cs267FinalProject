# cs267FinalProject
Parallelizing Gstat for CS267 Final Project

Pseudo code for serial implementation:
1. Read in measurements (value and coordinates)
2. Determine all the pairs of points
3. For each pair of points
	- Calculate squared difference of values
	- Calculate distance from coordinates
	- Add the squared differnce to the appropriate distance bin
4. Divide the sum of squared differences by the number of pairs in bin
5. Print output