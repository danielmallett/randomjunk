/* Takes an integer, n and a base, b on the command line. It outputs the Base b expansion of n */
#include <iostream>
#include <vector>
#include <cstdio>
using namespace std;

int main(int argc, char** argv){
	int n, b;
	int i;
	vector<int> revdigits;

	if (argc != 3){
		cerr << "./expansion n b" << endl;
		return 1;
	}

	sscanf(argv[1], "%d", &n);
	sscanf(argv[2], "%d", &b);
	
	while (n > 0){
		revdigits.push_back(n % b);
		n /= b;
	}

	for (i = revdigits.size()-1; i >= 0; i--){
		cout << revdigits[i];
	}
	cout << endl;

	return 0;
}
