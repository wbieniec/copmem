// MemoryFill.cpp : Defines the entry point for the console application.
//

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>

void displayHelp(char* name) {
	std::cout << "Memory fill.\n Use: " << name << " <bytes to fill\n";
	std::cout << name << " 1024\n";
	std::cout << name << " 53k\n";
	std::cout << name << " 53M\n";
	std::cout << name << " 53G\n";
}

int main(int ac, char* av[])
{
	std::cout.setf(std::ios_base::unitbuf);
	std::cout << "Memory fill to disable disk caches (--> 'cold start')\n";
	if (ac == 1) {
		displayHelp(av[0]);
		return 0;
	}

	int l = strlen(av[1]);
	size_t mult = 1;
	switch (av[1][l - 1]) {
	case 'k':
		mult = 1 << 10;
		break;
	case 'M':
		mult = 1 << 20;
		break;
	case 'G':
		mult = 1 << 30;
		break;
	}
	size_t N = mult * std::atoll(av[1]);
	if (N == 0) {
		displayHelp(av[0]);
		return 0;
	}
	std::cout << N << " bytes\n";
	const size_t num = 100;
	const size_t fillSize = N / num;
	char** arr = new char*[num];
	for (size_t i = 0; i < num; i++) {
		std::cout << ".";
		arr[i] = new char[fillSize];
		if (arr[i] == 0) {
			std::cerr << "out of memory";
			return -1;
		}
		memset(arr[i], 170, fillSize);
	}
	std::cout << " Done!\n";
	return 0;
}
