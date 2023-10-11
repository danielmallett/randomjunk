#include <iostream>
#include <cstdio>
using namespace std;

int main(int argc, char** argv){
	long long x, y, r, tmp, ans;
	bool verbose;

	verbose = false;

	if (argc != 4){
		cout << "usage: euclid a b y|n (verbose)" << endl;
		return 1;
	}

	if (sscanf(argv[1], "%lld", &x) != 1){
		cout << "usage: euclid a b y|n (verbose)" << endl;
		return 1;
	}

	if (sscanf(argv[2], "%lld", &y) != 1){
		cout << "usage: euclid a b y|n (verbose)" << endl;
		return 1;
	}

	if (argv[3][0] == 'y'){
		verbose = true;
	}


	if (y < x){
		tmp = x;
		x = y;
		y = tmp;
	}	

	r = y % x;

	if (verbose){
		printf("x = %-10lld y = %-10lld r = %-10lld \n", x, y, r);
	}

	while (r != 0){
		y = x;
		x = r;
		r = y % x;
		if (verbose){
			printf("x = %-10lld y = %-10lld r = %-10lld \n", x, y, r);
		}
	}

	cout << "GCD(" << x << ", " << y << ") = " << x << endl;
	return 0;
}
