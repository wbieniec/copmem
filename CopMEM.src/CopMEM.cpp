/* ================================================================= *
*  CopMEM.cpp : Main file                                           *
*                                                                   *
*  copMEM is a program for efficient computation of MEMs            *
*  (Maximal Exact Matches) in a pair of genomes.                    *
*  Its main algorithmic idea requires that two internal parameters  *
*  (k1 and k2) are coprime, hence the name.                         *
*                                                                   *
*                                                                   *
*  Copyright (c) 2018, Szymon Grabowski and Wojciech Bieniecki      *
*  All rights reserved                                              *
*                                                                   *
*  This program is free software: you can redistribute it and/or    *
*  modify it under the terms of the GNU General Public License as   *
*  published by the Free Software Foundation, either version 3 of   *
*  the License, or (at your option) any later version.              *
*                                                                   *
*  This program is distributed in the hope that it will be useful,  *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of   *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
*  GNU General Public License for more details.                     *
*                                                                   *
*  You should have received a copy of the GNU General Public        *
*  License along with this program.                                 *
*                                                                   *
*  This file is subject to the terms and conditions defined in the  *
*  file 'license', which is part of this source code package.       *
* ================================================================= */

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>
#include <cmath>
#include <cstring>
#include <cassert>

#if (defined(linux) || defined(__linux) || defined(__linux__))
#define _prefetch(x,y) __builtin_prefetch(x,1,(4-y))
#else
#include <xmmintrin.h>
#define _prefetch(x,y) _mm_prefetch(x,y)
#endif

#include "StopWatch.h"

//typedef std::uint32_t myuint;
typedef std::pair<std::string, size_t> SequenceItem;
typedef std::tuple<std::string, size_t, size_t> SequenceItem2;
typedef std::tuple<std::string, char*, size_t> SequenceItem3;
typedef std::vector<SequenceItem> Sequences;
typedef std::vector<SequenceItem2> Sequences2;
typedef std::tuple<size_t, char*, Sequences> GenomeData; //<size of buffer, memory pointer, starting pointer, seqnence list>
template<class MyUINT1, class MyUINT2>
 using HashBuffer =  std::pair<MyUINT1*, MyUINT2* >;


/////////////////// FUNCTIONS ///////////////////////////
SequenceItem3 readSequence(std::ifstream& f, SequenceItem2& seq, const char paddingChar, bool removeNs);
GenomeData readMultiFasta(std::string fn, const char paddingChar, bool removeNs);
void displayHelp(const char* progName);
void displayParams();
void processCmd(int argc, char* argv[]);
size_t getFileSize(std::string& fn);


//////////////////// GLOBALS ////////////////////////////
const int K = 44;
const int K_OVER_4 = K / 4;
const std::uint32_t HASH_SIZE = 1U << 29;
const std::uint32_t HASH_SIZE_MINUS_ONE = HASH_SIZE - 1;

int L = 100;
int k1 = 8;
int k2 = 7;

std::string matchesFN;
std::ofstream f_matches;
std::string R_FN;
std::string Q_FN;

enum verbosity{v0, v1, v2};
enum reverse{no, yes, both};


verbosity isVerbose = v1;
reverse isRC = no;
bool isTight = false;
//////////////////// GLOBALS ////////////////////////////


void displayHelp(const char* progName) {
	std::cout << "copMEM 0.1.1, by Szymon Grabowski and Wojciech Bieniecki, June 2018." << std::endl;
	std::cout << "Usage: " << progName << " [-l n] [-v|-q] [-t [-b]|[-r]] <-o MEMs_file> <Ref_genome> <Query_genome>\n";
	std::cout << "Attention: -o is a required parameter. l is optional (default: 100).\n";
	std::cout << "-o MEMs_file - REQUIRED parameter. Output files with matches\n";
	std::cout << "-v - verbose mode. Display more details\n";
	std::cout << "-q - quiet mode. No screen output\n";
	std::cout << "-t - tight mode. Saves memory. May slow down I/O.\n";
	std::cout << "-b - compute forward and reverse complement matches. Available only with -t.\n";
	std::cout << "-r - compute only reverse complement matches. Available only with -t.\n";
	std::cout << "-l n - minimal length of matches. Default value is 100.\n";
}


void displayParams() {
	std::cout << "PARAMETERS" << std::endl;
	std::cout << "Reference filename: " << R_FN << std::endl;
	std::cout << "Query filename: " << Q_FN << std::endl;
	std::cout << "l = " << L << std::endl;
	std::cout << "K = " << K << std::endl;
	std::cout << "HASH_SIZE = " << HASH_SIZE << std::endl;
	std::cout << "k1 = " << k1 << std::endl;
	std::cout << "k2 = " << k2 << std::endl;
	std::cout << "Memory save = " << (isTight? "Yes":"No") << std::endl;
	std::cout << "Reverse Complement = " << ((isRC==both)? "Both" : (isRC == yes) ? "Yes" : "No") << std::endl;
}


void processCmd(int argc, char* argv[]) {
	assert(K % 4 == 0);
	bool isOset = false;
	const char* incompleteCmd = " option requires one integer argument.\n";
	for (int i = 1; i < argc - 2; ++i) {
		std::string arg = argv[i];

		if (arg == "-o") {
			if (i + 1 < argc) {
				matchesFN = argv[++i];
				isOset = true;
			}
			else {
				std::cerr << "-o requires file name.\n";
				exit(1);
			}
		}
		if (arg == "-v") {
			if (isVerbose == v0) {
				std::cerr << "-v and -q parameters are exclusive.\n";
				exit(1);
			}
			isVerbose = v2;
		}
		if (arg == "-q") {
			if (isVerbose == v2) {
				std::cerr << "-v and -q parameters are exclusive.\n";
				exit(1);
			}
			isVerbose = v0;
		}
		if (arg == "-t") {
			isTight = true;
		}
		if (arg == "-b") {
			if (isRC == yes) {
				std::cerr << "-b and -r parameters are exclusive.\n";
				exit(1);
			}
			isRC = both;
		}
		if (arg == "-r") {
			if (isRC == both) {
				std::cerr << "-b and -r parameters are exclusive.\n";
				exit(1);
			}
			isRC = yes;
		}

		if (arg == "-h") {
			displayHelp("");
			exit(0);
		}
		if (arg == "-l") {
			if (i + 1 < argc) {
				L = atoi(argv[++i]);
				if (L < 50) {
					std::cerr << "Incorrect L value (must be >= 50).\n";
					exit(1);
				}
			}
			else {
				std::cerr << "L" << incompleteCmd;
				exit(1);
			}
		}
	}

	if (isOset == false) {
		std::cerr << "-o not given or specified correctly.\n";
		exit(1);
	}
	if (!isTight && (isRC!=no)) {
		std::cerr << "-r or -b parameters go with -t only.\n";
		exit(1);
	}

	//touching files:
	R_FN = argv[argc - 2];
	Q_FN = argv[argc - 1];
	std::ifstream f;

	f.open(R_FN);
	if (f.fail()) {
		std::cerr << "\nReference file '" << R_FN << "' does not exist. Quit.";
		exit(1);
	}
	f.close();
	f.open(Q_FN);
	if (f.fail()) {
		std::cerr << "\nQuery file '" << Q_FN << "' does not exist. Quit.";
		exit(1);
	}
	f.close();
	f_matches.open(matchesFN);
	if (f_matches.fail()) {
		std::cerr << "\nFAILED. The -o parameter specifies a file that cannot be created.\n";
		exit(1);
	}
	f_matches.close();

	/* setting k1 and k2 */
	int tempVar = L - K + 1;
	k1 = (int)(pow(tempVar, 0.5)) + 1;
	k2 = k1 - 1;
	if (k1 * k2 > tempVar) {
		--k1;
		--k2;
	}
}


bool arrayComp(size_t* p1, size_t* p2) {
	if (p1[1] < p2[1]) return true;
	if (p1[1] > p2[1]) return false;
	return p1[0] < p2[0];
}


// based on http://www.amsoftware.narod.ru/algo2.html
inline std::uint32_t maRushPrime1HashSimplified(const char *str) {
	std::uint64_t hash = K;
	for (std::uint32_t j = 0; j < K_OVER_4; ) {
		std::uint32_t k;
		memcpy(&k, str, 4);
		k += j++;
		hash ^= k;
		hash *= 171717;
		str += 4;
	}
	return (std::uint32_t)(hash & HASH_SIZE_MINUS_ONE);
}

template <class MyUINT>
void genCumm(GenomeData& genome, MyUINT* cumm, std::uint32_t(*hashFunc)(const char*)) {
	const size_t MULTI1 = 128;
	const size_t k1MULTI1 = k1 * MULTI1;

	size_t N = std::get<0>(genome);
	char* gen = std::get<1>(genome);

	std::fill(cumm, cumm + HASH_SIZE + 2, 0);

	uint32_t hashPositions[MULTI1];
	size_t i;

	//for (i = 0; i < N - K - k1MULTI1; i += k1MULTI1) {
	for (i = 0; i + K + k1MULTI1 < N ; i += k1MULTI1) {
		char* tempPointer = gen + i;
		for (size_t temp = 0; temp < MULTI1; ++temp) {
			hashPositions[temp] = hashFunc(tempPointer) + 2;
			tempPointer += k1;
			_prefetch((char*)(cumm + hashPositions[temp]), 1);
		}

		for (size_t temp = 0; temp < MULTI1; ++temp) {
			++cumm[hashPositions[temp]];
		}
	}

	//////////////////// processing the end part of R  //////////////////////
	for (; i < N - K + 1; i += k1) {
		uint32_t h = hashFunc(gen + i) + 2;
		++cumm[h];
	}
	//////////////////// processing the end part of R //////////////////////

	std::partial_sum(cumm, cumm + HASH_SIZE + 1, cumm);
}


void dumpMEM(SequenceItem& item1, SequenceItem& item2, size_t* match, std::string &s) {
	size_t baseindex1 = match[0] - item1.second;
	size_t baseindex2 = match[1] -item2.second;
	s.append(" ");
	s.append(item1.first);
	s.append("\t");
	s.append(std::to_string(baseindex1));
	s.append("\t");
	s.append(std::to_string(baseindex2));
	s.append("\t");
	s.append(std::to_string(match[2]));
	s.append("\n");
	f_matches << s;
	s.clear();
}

void dumpMEMTight(SequenceItem& item1, size_t* match, std::string &s, size_t counter) {
	size_t baseindex1 = match[0] - item1.second;
	size_t baseindex2 = match[1] - counter;
	s.append(" ");
	s.append(item1.first);
	s.append("\t");
	s.append(std::to_string(baseindex1));
	s.append("\t");
	s.append(std::to_string(baseindex2));
	s.append("\t");
	s.append(std::to_string(match[2]));
	s.append("\n");
	f_matches << s;
	s.clear();
}



bool lowerBoundComp(SequenceItem lhs, SequenceItem rhs) {
	return lhs.second < rhs.second;
}


SequenceItem findSeqDesc(size_t index, Sequences& seq) {
	SequenceItem dummySequenceItem = { "", index };
	SequenceItem item = seq[0];
	auto lower = std::lower_bound(seq.begin(), seq.end(), dummySequenceItem, lowerBoundComp);
	size_t diff = lower - seq.begin();
	return seq[diff - 1];
}

void displayMatchInfo(std::string& name, size_t count) {
	if (isVerbose == v2) {
		switch (count) {
		case 0:  std::cout << name << ": no matches.\n"; break;
		case 1:  std::cout << name << ": match.\n"; break;
		default: std::cout << name << ": matches.\n"; break;
		}
	}
}


void postProcess(std::vector<size_t*> &matches, GenomeData& t1, GenomeData& t2, size_t index2) {
	Sequences& seq2 = std::get<2>(t2);
	SequenceItem seq2item = seq2[index2];
	f_matches << "> " << seq2item.first << std::endl;

	if (matches.size() == 0) {
		displayMatchInfo(seq2item.first, 0);
		return;
	}
	char* gen1 = std::get<1>(t1);
	Sequences& seq1 = std::get<2>(t1);
	if (matches.size() > 1)
		std::sort(matches.begin(), matches.end(), arrayComp);

	SequenceItem seq1item;
	std::string dumpString;
	dumpString.reserve(256);

	size_t* prev_match = matches[0];
	seq1item = findSeqDesc(prev_match[0], seq1);

	auto foundPos = std::find(gen1 + prev_match[0], gen1 + prev_match[0] + prev_match[2], 'N');
	if (foundPos == gen1 + prev_match[0] + prev_match[2])
		dumpMEM(seq1item, seq2item, prev_match, dumpString);

	std::uint64_t count = 1ULL;

	for (auto match: matches) {
		if (memcmp(prev_match, match, 3 * sizeof(size_t)) != 0) {
			auto foundPos = std::find(gen1 + match[0], gen1 + match[0] + match[2], 'N');
			if (foundPos != gen1 + match[0] + match[2])
				continue;
			if (prev_match[0] != match[0])
				seq1item = findSeqDesc(match[0], seq1);
			dumpMEM(seq1item, seq2item, match, dumpString);
			++count;
		}
		prev_match = match;
	}
	displayMatchInfo(seq2item.first, count);
	for (auto match : matches)
		delete[] match;
	matches.clear();
}

void postProcessTight(std::vector<size_t*> &matches, GenomeData& t1, std::string seqName, size_t counter) {
	f_matches << "> " << seqName << std::endl;

	if (matches.size() == 0) {
		displayMatchInfo(seqName, 0);
		return;
	}
	char* gen1 = std::get<1>(t1);
	Sequences& seq1 = std::get<2>(t1);
	if (matches.size() > 1)
		std::sort(matches.begin(), matches.end(), arrayComp);

	SequenceItem seq1item;
	std::string dumpString;
	dumpString.reserve(256);

	size_t* prev_match = matches[0];
	seq1item = findSeqDesc(prev_match[0], seq1);

	auto foundPos = std::find(gen1 + prev_match[0], gen1 + prev_match[0] + prev_match[2], 'N');
	if (foundPos == gen1 + prev_match[0] + prev_match[2])
		dumpMEMTight(seq1item,  prev_match, dumpString, counter);

	std::uint64_t count = 1ULL;

	for (auto match : matches) {
		if (memcmp(prev_match, match, 3 * sizeof(size_t)) != 0) 
		{
			auto foundPos = std::find(gen1 + match[0], gen1 + match[0] + match[2], 'N');
			if (foundPos != gen1 + match[0] + match[2])
				continue;
			if (prev_match[0] != match[0])
				seq1item = findSeqDesc(match[0], seq1);
			dumpMEMTight(seq1item,  match, dumpString, counter);
			++count;
		}
		prev_match = match;
	}
	displayMatchInfo(seqName, count);
	for (auto match : matches)
		delete[] match;
	matches.clear();
}


template<class MyUINT1, class MyUINT2>
HashBuffer<MyUINT1, MyUINT2> processRef(GenomeData& rGenome, std::uint32_t(*hashFunc)(const char*)) {
	CStopWatch stopwatch;
	stopwatch.start();
	size_t N = std::get<0>(rGenome);
	char* start = std::get<1>(rGenome);

	const size_t hashCount = (N - K + 1) / k1;
	const unsigned int MULTI2 = 128;
	const unsigned int k1MULTI2 = k1 * MULTI2;

	if (isVerbose>v0)
		std::cout << "Hash count = " << hashCount << std::endl;

	MyUINT1* sampledPositions = new MyUINT1[hashCount];
	MyUINT2* cumm = new MyUINT2[HASH_SIZE + 2];
	if (isVerbose>v1) {
		std::cout << "\tprocessRef: init = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
	}
	genCumm(rGenome, cumm, hashFunc);

	if (isVerbose>v1) {
		std::cout << "\tprocessRef: cumm(1) = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
	}
	uint32_t hashPositions[MULTI2];
	MyUINT1 i1;

	for (i1 = 0; i1 + K + k1MULTI2 < N ; i1 += k1MULTI2) {
		char* tempPointer = start + i1;
		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
			hashPositions[temp] = hashFunc(tempPointer) + 1;
			tempPointer += k1;
			_prefetch((char*)(cumm + hashPositions[temp]), 1);
		}

		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
			_prefetch((char*)(sampledPositions + *(cumm + hashPositions[temp])), 1);
		}

		MyUINT1 i2 = i1;
		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
			sampledPositions[cumm[hashPositions[temp]]] = i2;
			++cumm[hashPositions[temp]];
			i2 += k1;
		}
	}

	//////////////////// processing the end part of R
	for (; i1 < N - K + 1; i1 += k1) {
		char* tempPointer = start + i1;
		uint32_t h = hashFunc(start + i1) + 1;
		sampledPositions[cumm[h]] = i1;
		++cumm[h];
	}
	if (isVerbose>v1)
		std::cout << "\tprocessRef: cumm(2) = " << stopwatch.stop() << std::endl;

	return { sampledPositions, cumm };
}


template<class MyUINT1, class MyUINT2>
void processQuery(HashBuffer<MyUINT1, MyUINT2> buffer, GenomeData& rGenome, GenomeData& qGenome, std::uint32_t(*hashFunc)(const char*)) {
	const int L_PLUS_ONE = L + 1;
	const int L2 = L / 2;
	const int LK2 = (L - K) / 2;
	const int LK2_MINUS_2 = LK2 - 2;
	const int LK2_MINUS_4 = LK2 - 4;
	const int K_PLUS_LK22 = K + LK2_MINUS_2;
	const int K_PLUS_LK24 = K + LK2_MINUS_4;

	size_t N1 = std::get<0>(rGenome);
	size_t N2 = std::get<0>(qGenome);
	char* start1 = std::get<1>(rGenome);
	char* start2 = std::get<1>(qGenome);
	Sequences seq1 = std::get<2>(rGenome);
	Sequences seq2 = std::get<2>(qGenome);

	const size_t hashCount = (N1 - K + 1) / k1;
	const std::uint32_t MULTI = 256;
	const std::uint32_t k2MULTI = k2 * MULTI;

	char* curr1 = start1;
	char* curr2 = start2;

	MyUINT1* sampledPositions = buffer.first;
	MyUINT2* cumm = buffer.second;

	std::vector<size_t*> matches;

	std::uint32_t hArray[MULTI];
	MyUINT2 posArray[MULTI * 2];

	size_t currSeqIndex = 0;
	size_t nextSeqBeginPos = (seq2.size() > 1 ? seq2[currSeqIndex + 1].second : SIZE_MAX);

	f_matches.open(matchesFN);

	std::uint32_t l1, l2, r1, r2;

	size_t i1;
	for (i1 = 0; i1 + K + k2MULTI < N2 + 1 ; i1 += k2MULTI) {
		size_t i1temp = i1;
		char* prev = start2 + i1;
		curr2 = start2 + i1;
		size_t tempCount = 0;
		for (size_t i2 = 0; i2 < MULTI; ++i2) {
			if (*curr2 != 'N') {
				hArray[tempCount++] = hashFunc(curr2);
			}
			curr2 += k2;
		}
		for (size_t i2 = 0; i2 < tempCount; ++i2) {
			memcpy(posArray + i2 * 2, cumm + hArray[i2], sizeof(MyUINT2) * 2);
		}

		curr2 = start2 + i1;	// !!!
		for (size_t i2 = 0; i2 < tempCount; ++i2) {
			if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
				curr2 += k2;
				continue;
			}

			memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
			memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

			for (MyUINT1 j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
				curr1 = start1 + sampledPositions[j];

				memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
				memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

				if (r1 == r2 || l1 == l2) {
					char* p1 = curr1 + K;
					char* p2 = curr2 + K;

					while (*p1 == *p2) {
						++p1;
						++p2;
					}
					char* right = p1;

					p1 = curr1 - 1;
					p2 = curr2 - 1;
					while (*p1 == *p2) {
						--p1;
						--p2;
					}

					if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
						size_t* tempMatch = new size_t[3];
						tempMatch[0] = p1 + 1 - start1;
						tempMatch[1] = p2 + 1 - start2;
						tempMatch[2] = right - p1 - 1;
						while ((size_t)(p2 + 1 - start2) >= nextSeqBeginPos) {
							if (currSeqIndex + 1 < seq2.size()) {
								postProcess(matches, rGenome, qGenome, currSeqIndex);
								++currSeqIndex;
								nextSeqBeginPos = seq2[currSeqIndex + 1].second;
							}
							else
								break;
						}
						matches.push_back(tempMatch);

					}
				}
			}
			curr2 += k2;
		}
	}

	//////////////////// processing the end part of Q  //////////////////////
	for (; i1 + K < N2 + 1; i1 += k2) {
		char* prev = start2 + i1;
		curr2 = start2 + i1;
		size_t tempCount = 0;

		if (*curr2 != 'N') {
			memcpy(posArray + tempCount * 2, cumm + hashFunc(curr2), sizeof(MyUINT2) * 2);
			++tempCount;
		}

		curr2 = start2 + i1;
		for (size_t i2 = 0; i2 < tempCount; ++i2) {
			if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
				curr2 += k2;
				continue;
			}

			memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
			memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

			for (MyUINT1 j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
				curr1 = start1 + sampledPositions[j];

				memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
				memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

				if (r1 == r2 || l1 == l2) {
					char* p1 = curr1 + K;
					char* p2 = curr2 + K;

					while (*p1 == *p2) {
						++p1;
						++p2;
					}
					char* right = p1;

					p1 = curr1 - 1;
					p2 = curr2 - 1;
					while (*p1 == *p2) {
						--p1;
						--p2;
					}

					if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
						size_t* tempMatch = new size_t[3];
						tempMatch[0] = p1 + 1 - start1;
						tempMatch[1] = p2 + 1 - start2;
						tempMatch[2] = right - p1 - 1;

						while ((size_t)(p2 + 1 - start2) > nextSeqBeginPos) {
							if (currSeqIndex + 1 < seq2.size()) {
								postProcess(matches, rGenome, qGenome, currSeqIndex);
								++currSeqIndex;
								nextSeqBeginPos = seq2[currSeqIndex + 1].second;
							}
							else
								break;
						}
						matches.push_back(tempMatch);
					}
				}
			}
			curr2 += k2;
		}
	}
	//////////////////// processing the end part of Q  //////////////////////
	while (currSeqIndex + 2 <= seq2.size()) {
		postProcess(matches, rGenome, qGenome, currSeqIndex);
		++currSeqIndex;
		nextSeqBeginPos = seq2[currSeqIndex + 1].second;
	}
	postProcess(matches, rGenome, qGenome, currSeqIndex);
	f_matches.close();
}

void reverseComplement(char* start, const char* lut, const std::size_t N) {
	char* left = start + 1; //sequence starts from paddingChar
	char* right = start + N - 1;
	while (right > left) {
		char tmp = lut[*left];
		*left = lut[*right];
		*right = tmp;
		++left;
		--right;
	}
}

void reverseComplement_new(char* start, const char* lut, const std::size_t N) {
	char* left = start + 1; //sequence starts from paddingChar
	char* right = start + N - 1;

	std::reverse(left, right + 1);
	for (size_t i = 1; i < N; ++i)
		*(start + i) = lut[*(start + i)];
}

template<class MyUINT1, class MyUINT2>
void processQueryTight(HashBuffer<MyUINT1, MyUINT2> buffer, GenomeData& rGenome, Sequences2& seq2, std::uint32_t(*hashFunc)(const char*)) {
	const int L_PLUS_ONE = L + 1;
	const int L2 = L / 2;
	const int LK2 = (L - K) / 2;
	const int LK2_MINUS_2 = LK2 - 2;
	const int LK2_MINUS_4 = LK2 - 4;
	const int K_PLUS_LK22 = K + LK2_MINUS_2;
	const int K_PLUS_LK24 = K + LK2_MINUS_4;
	const int paddingSize = L_PLUS_ONE;
	const char paddingChar = 125;
	const bool removeNs = true;

	size_t N1 = std::get<0>(rGenome);

	char* start1 = std::get<1>(rGenome);
	Sequences seq1 = std::get<2>(rGenome);
	const size_t hashCount = (N1 - K + 1) / k1;
	const unsigned int MULTI = 256;
	const unsigned int k2MULTI = k2 * MULTI;

	char* curr1 = start1;
	char complement[256];
	memset(complement, paddingChar, 256);
	complement['A'] = 'T';
	complement['C'] = 'G';
	complement['G'] = 'C';
	complement['T'] = 'A';
	complement['N'] = 'N';


	MyUINT1* sampledPositions = buffer.first;
	MyUINT2* cumm = buffer.second;

	std::vector<size_t*> matches;

	uint32_t hArray[MULTI];
	MyUINT2 posArray[MULTI * 2];

	f_matches.open(matchesFN);
	std::ifstream qf(Q_FN, std::ios::binary);

	std::uint32_t l1, l2, r1, r2;

	size_t i1;


	size_t counter = 0ULL;

	for (auto se2i : seq2) {
	//Read sequence to memory from file
		SequenceItem3 si3 = readSequence(qf, se2i, paddingChar, removeNs);
		size_t N2 = std::get<2>(si3);
		char* start2 = std::get<1>(si3);
		if (isRC!=yes) {
			char* curr2 = start2;
			i1 = 0;
			for (i1 = 0; i1 + K + k2MULTI < N2 + 1; i1 += k2MULTI) {
				size_t i1temp = i1;
				char* prev = start2 + i1;
				curr2 = start2 + i1;
				size_t tempCount = 0;
				for (size_t i2 = 0; i2 < MULTI; ++i2) {
					hArray[tempCount++] = hashFunc(curr2);
					curr2 += k2;
				}
				for (size_t i2 = 0; i2 < tempCount; ++i2) {
					memcpy(posArray + i2 * 2, cumm + hArray[i2], sizeof(MyUINT2) * 2);
				}

				curr2 = start2 + i1;	// !!!
				for (size_t i2 = 0; i2 < tempCount; ++i2) {
					if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
						curr2 += k2;
						continue;
					}

					memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
					memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

					for (MyUINT1 j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
						curr1 = start1 + sampledPositions[j];

						memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
						memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

						if (r1 == r2 || l1 == l2) {
							char* p1 = curr1 + K;
							char* p2 = curr2 + K;

							while (*p1 == *p2) {
								++p1;
								++p2;
							}
							char* right = p1;

							p1 = curr1 - 1;
							p2 = curr2 - 1;
							while (*p1 == *p2) {
								--p1;
								--p2;
							}

							if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
								size_t* tempMatch = new size_t[3];
								tempMatch[0] = p1 + 1 - start1;
								tempMatch[1] = p2 + 1 - start2 + counter;
								tempMatch[2] = right - p1 - 1;
								matches.push_back(tempMatch);
							}
						}
					}
					curr2 += k2;
				}
			}
			//////////////////// processing the end part of Q  //////////////////////
			for (; i1 + K < N2 + 1; i1 += k2) {
				char* prev = start2 + i1;
				curr2 = start2 + i1;
				size_t tempCount = 0;

				memcpy(posArray + tempCount * 2, cumm + hashFunc(curr2), sizeof(MyUINT2) * 2);
				++tempCount;

				curr2 = start2 + i1;
				for (size_t i2 = 0; i2 < tempCount; ++i2) {
					if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
						curr2 += k2;
						continue;
					}

					memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
					memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

					for (MyUINT1 j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
						curr1 = start1 + sampledPositions[j];

						memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
						memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

						if (r1 == r2 || l1 == l2) {
							char* p1 = curr1 + K;
							char* p2 = curr2 + K;

							while (*p1 == *p2) {
								++p1;
								++p2;
							}
							char* right = p1;

							p1 = curr1 - 1;
							p2 = curr2 - 1;
							while (*p1 == *p2) {
								--p1;
								--p2;
							}

							if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
								size_t* tempMatch = new size_t[3];
								tempMatch[0] = p1 + 1 - start1;
								tempMatch[1] = p2 + 1 - start2 + counter;
								tempMatch[2] = right - p1 - 1;
								//tempMatch[2] &= (~0x80000000ULL);
								matches.push_back(tempMatch);
							}
						}
					}
					curr2 += k2;
				}
			}
			//////////////////// processing the end part of Q  //////////////////////
			postProcessTight(matches, rGenome, std::get<0>(se2i), counter);
		}
		if (isRC!=no) {
			reverseComplement(start2, complement, N2);
			char* curr2 = start2;
			i1 = 0;
			for (i1 = 0; i1 + K + k2MULTI < N2 + 1; i1 += k2MULTI) {
				size_t i1temp = i1;
				char* prev = start2 + i1;
				curr2 = start2 + i1;
				size_t tempCount = 0;
				for (size_t i2 = 0; i2 < MULTI; ++i2) {
					hArray[tempCount++] = hashFunc(curr2);
					curr2 += k2;
				}
				for (size_t i2 = 0; i2 < tempCount; ++i2) {
					memcpy(posArray + i2 * 2, cumm + hArray[i2], sizeof(MyUINT2) * 2);
				}

				curr2 = start2 + i1;	// !!!
				for (size_t i2 = 0; i2 < tempCount; ++i2) {
					if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
						curr2 += k2;
						continue;
					}

					memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
					memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

					for (MyUINT1 j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
						curr1 = start1 + sampledPositions[j];

						memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
						memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

						if (r1 == r2 || l1 == l2) {
							char* p1 = curr1 + K;
							char* p2 = curr2 + K;

							while (*p1 == *p2) {
								++p1;
								++p2;
							}
							char* right = p1;

							p1 = curr1 - 1;
							p2 = curr2 - 1;
							while (*p1 == *p2) {
								--p1;
								--p2;
							}

							if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
								size_t* tempMatch = new size_t[3];
								tempMatch[0] = p1 + 1 - start1;
								tempMatch[1] = (p2 + 1 - start2) + counter;
								tempMatch[2] = right - p1 - 1;
								matches.push_back(tempMatch);
							}
						}
					}
					curr2 += k2;
				}
			}
			//////////////////// processing the end part of Q  //////////////////////
			for (; i1 + K < N2 + 1; i1 += k2) {
				char* prev = start2 + i1;
				curr2 = start2 + i1;
				size_t tempCount = 0;

				memcpy(posArray + tempCount * 2, cumm + hashFunc(curr2), sizeof(MyUINT2) * 2);
				++tempCount;

				curr2 = start2 + i1;
				for (size_t i2 = 0; i2 < tempCount; ++i2) {
					if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
						curr2 += k2;
						continue;
					}

					memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
					memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

					for (MyUINT1 j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
						curr1 = start1 + sampledPositions[j];

						memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
						memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

						if (r1 == r2 || l1 == l2) {
							char* p1 = curr1 + K;
							char* p2 = curr2 + K;

							while (*p1 == *p2) {
								++p1;
								++p2;
							}
							char* right = p1;

							p1 = curr1 - 1;
							p2 = curr2 - 1;
							while (*p1 == *p2) {
								--p1;
								--p2;
							}

							if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
								size_t* tempMatch = new size_t[3];
								tempMatch[0] = p1 + 1 - start1;
								tempMatch[1] = (p2 + 1 - start2) + counter;
								tempMatch[2] = right - p1 - 1;
								matches.push_back(tempMatch);
							}
						}
					}
					curr2 += k2;
				}
			}
			//////////////////// processing the end part of Q  //////////////////////
			postProcessTight(matches, rGenome, std::get<0>(se2i) + " Reverse", counter);
		}

		delete[](start2 - paddingSize);

		counter += N2;
	}
	f_matches.close();
}


// invoked only in readWholeFAFile function
void replaceBadSymbol(char* gen, char* dst, char symbol, char paddingChar) {
	char* movingPtr = gen;
	while (1) {
		char* tempPtr = std::find(movingPtr, dst, symbol);
		while (*tempPtr == symbol) {
			*tempPtr = paddingChar;
			++tempPtr;
		}
		if (tempPtr == dst)
			break;
		movingPtr = tempPtr + 1;
	}
}


GenomeData readMultiFasta(std::string fn, const char paddingChar, bool removeNs) {
	CStopWatch stopWatch;
	stopWatch.start();
	const char beginChar = '>';
	const char terminatorChar = 0;
	const char eolChar1 = 10;
	const char eolChar2 = 13;
	const char spaceChar = ' ';
	const char filterChar = (char)0xDF;
	const int paddingSize = L + 1;

	//create a buffer for the whole file + padding at left and right
	std::ifstream f(fn, std::ios::ate | std::ios::binary);
	if (f.fail()) {
		std::cerr << "\nFile '" << fn << "' does not exist. Quit.";
		exit(1);
	}
	size_t N = f.tellg();
	f.seekg(0, std::ios::beg);
	char* buf1 = new char[N + 2 * paddingSize];
	if (buf1 == nullptr) {
		std::cerr << "\nFile '" << fn << "' is too large. Quit.";
		exit(1);
	}
	memset(buf1, paddingChar, paddingSize);
	f.read(buf1 + paddingSize, N);
	if (isVerbose>v1) {
		std::cout << "\treadMultiFasta: Reading file from disk " << stopWatch.stop() << "\n";
		stopWatch.resume();
	}
	buf1[paddingSize + N] = terminatorChar; // null-terminate the string
	memset(buf1 + paddingSize + N + 1, paddingChar, paddingSize - 1);
	memset(buf1 + paddingSize + N, terminatorChar, 10);  // >= sizeof(uint_64) is enough
	f.close();

	char* gen = buf1 + paddingSize;
	Sequences seq;

	char* dst = gen;
	char* src = gen;

	char tempLine[512];  // FASTA lines shouldn't be that long (even 128 should be ok)

	while (1) {
		if (*src == beginChar) {
			size_t idx = 0;
			while (*src != eolChar1 && *src != eolChar2 && *src != ' ' && *src != terminatorChar) {
				tempLine[idx] = *src++;
				++idx;
			}
			tempLine[idx] = 0;
			seq.push_back({ tempLine + 1, (dst - gen) });  // + 1, as we omit the starting '>'
														   //search for EOL
			while (*src != eolChar1 && *src != eolChar2)
				src++;

			*dst++ = paddingChar;
		}
		else {
			while (*src == eolChar1 || *src == eolChar2) {
				++src;
			}
			if (*src == beginChar)
				continue;
			if (*src == terminatorChar)
				break;

			uint64_t temp2;
			memcpy(&temp2, src, 8);
			while (((~temp2) & 0x4040404040404040) == 0) {
				temp2 = temp2 & 0xDFDFDFDFDFDFDFDF;
				memcpy(dst, &temp2, 8);
				dst += 8;
				src += 8;
				memcpy(&temp2, src, 8);
			}
			while (((~(*src)) & (char)0x40) == 0) {
				*dst++ = (*src++) & (char)0xDF;
			}
		}
	}

	memset(dst, paddingChar, N + paddingSize - (dst - gen));

	if (removeNs == true) {
		replaceBadSymbol(gen, dst, 'N', paddingChar);
		replaceBadSymbol(gen, dst, 'n', paddingChar);
	}
	if (isVerbose>v1) {
		std::cout << "\treadMultiFasta: Analysing file " << stopWatch.stop() << "\n";
		std::cout << "\treadMultiFasta: Genome data size " << (dst - gen) << std::endl;
		std::cout << "\treadMultiFasta: Sequences " << seq.size() << std::endl;
	}
	return { (dst - gen), gen, seq };
}



SequenceItem3 readSequence(std::ifstream& f, SequenceItem2& seq, const char paddingChar, bool removeNs) {
	const char terminatorChar = 0;
	const char eolChar1 = 10;
	const char eolChar2 = 13;
	const int paddingSize = L + 1;

	size_t N = std::get<2>(seq);
	f.seekg(std::get<1>(seq), std::ios::beg);
	char* buf1 = new char[N + 2 * paddingSize + 1];
	if (buf1 == nullptr) {
		std::cerr << "\nSequence '" << std::get<0>(seq) << "' is too large. Quit.";
		exit(1);
	}
	memset(buf1, paddingChar, paddingSize); // reading starts from padding char
	char* gen = buf1 + paddingSize;
	gen[0] = paddingChar;
	f.read(gen + 1, N);
	memset(gen + 1 + N, terminatorChar, 10);  // >= sizeof(uint_64) is enough
	memset(gen + 1 + N + 10, paddingChar, paddingSize - 10);
	

	
	char* dst = gen;
	char* src = gen;

	while (1) {
		while (*src == eolChar1 || *src == eolChar2) {
			++src;
		}
		if (*src == terminatorChar)
			break;

		uint64_t temp2;
		memcpy(&temp2, src, 8);
		while (((~temp2) & 0x4040404040404040) == 0) {
			temp2 = temp2 & 0xDFDFDFDFDFDFDFDF;
			memcpy(dst, &temp2, 8);
			dst += 8;
			src += 8;
			memcpy(&temp2, src, 8);
		}
		while (((~(*src)) & (char)0x40) == 0) {
			*dst++ = (*src++) & (char)0xDF;
		}
	}
	memset(dst, paddingChar, N + paddingSize - (dst - gen));

	if (removeNs == true) {
		replaceBadSymbol(gen, dst, 'N', paddingChar);
		replaceBadSymbol(gen, dst, 'n', paddingChar);
	}
	return { std::get<0>(seq),gen, (dst - gen) };
}


std::tuple<Sequences2, std::size_t,std::size_t> scanMultiFasta(std::string fn){
/* returns: <vector of sequence positions, file size, maximal length of the sequence>
*/
	CStopWatch stopWatch;
	stopWatch.start();
	const char beginChar = '>';
	const char terminatorChar = 0;
	const char eolChar1 = 10;
	const char eolChar2 = 13;
	const char spaceChar = ' ';
	const char filterChar = (char)0xDF;

	//create a buffer for the whole file + padding at left and right
	std::ifstream f(fn, std::ios::ate | std::ios::binary);
	if (f.fail()) {
		std::cerr << "\nFile '" << fn << "' does not exist. Quit.";
		exit(1);
	}
	size_t N = f.tellg();
	f.seekg(0, std::ios::beg);
	char* buf1 = new char[N];
	if (buf1 == nullptr) {
		std::cerr << "\nFile '" << fn << "' is too large. Quit.";
		exit(1);
	}
	f.read(buf1, N);
	if (isVerbose>v1) {
		std::cout << "\tscanMultiFasta: Reading file from disk " << stopWatch.stop() << "\n";
		stopWatch.resume();
	}
	f.close();
	char* gen = buf1;
	Sequences2 seq;
	char* src = gen;

	char tempLine[512];  // FASTA lines shouldn't be that long (even 128 should be ok)

	while (1) {
		if (*src == beginChar) {
			size_t idx = 0;
			size_t seqBeg = src - gen;  // offset of the first symbol of the seq name, just after the '>'
			while (*src != eolChar1 && *src != eolChar2 && *src != ' ' && *src != terminatorChar) {
				tempLine[idx] = *src++;
				++idx;
			}
			tempLine[idx] = 0;
			
			while (*src != eolChar1 && *src != eolChar2)
				src++;
			while (*src == eolChar1 || *src == eolChar2)
				src++;
			seq.push_back({ tempLine + 1, src - gen, seqBeg});
		}
		else {
			while (*src == eolChar1 || *src == eolChar2) {
				++src;
			}
			if (*src == beginChar)
				continue;
			if (*src == terminatorChar)
				break;

			uint64_t temp2;
			memcpy(&temp2, src, 8);
			while (((~temp2) & 0x4040404040404040) == 0) {
				temp2 = temp2 & 0xDFDFDFDFDFDFDFDF;
				src += 8;
				memcpy(&temp2, src, 8);
			}
			while (((~(*src)) & (char)0x40) == 0) {
				++src;
			}
		}
	}
	delete[] buf1;
	seq.push_back({ "dummy_end_sequence", 0, N });

	size_t seqSize = seq.size();
	//rearranging sequnece content:  Name, header_begin, sequence_begin  --> Name, sequence_begin, sequence_len
	size_t maxSeqLen = 0ULL;
	for (size_t i = 0; i < seqSize - 1; i++)
	{
		SequenceItem2 se2i = seq[i];
		SequenceItem2 se2i1 = seq[i + 1];
		size_t le = std::get<2>(se2i1) - std::get<1>(se2i);
		std::get<2>(seq[i]) = le;
		if (maxSeqLen < le)
			maxSeqLen = le;
	}
	seq.pop_back(); //remove dummy sequnece item
	if (isVerbose>v1) {
		std::cout << "\tscanMultiFasta: Analysing file " << stopWatch.stop() << "\n";
		std::cout << "\tscanMultiFasta: Sequences " << seq.size() << std::endl;
	}
	return { seq, N, maxSeqLen };
}

size_t getFileSize(std::string& fn) {
	std::ifstream f(fn, std::ios::ate | std::ios::binary);
	if (f.fail()) {
		std::cerr << "\nFile '" << fn << "' does not exist. Quit.";
		exit(1);
	}
	size_t N = f.tellg();
	f.close();
	return N;
}
template <class MyUINT1, class MyUINT2>
void deleteHashBuffer(HashBuffer<MyUINT1, MyUINT2> & buf) {
	delete[]buf.first;
	delete[]buf.second;
}

void deleteReading(GenomeData & r) {
	delete[](std::get<1>(r)-(L+1));
}


int main(int argc, char* argv[]) {
	std::cout.setf(std::ios_base::unitbuf);
	if (argc < 3) {
		displayHelp(argv[0]);
		return 1;
	}
	processCmd(argc, argv);
	if(isVerbose>v0)
		displayParams();

	CStopWatch stopwatch;
	stopwatch.start();
	Sequences2 se2;
	std::size_t maxSeqSize = 0;
	std::size_t N1 = 0;
	std::size_t N2 = 0;
	GenomeData  rGenome, qGenome;
	if (isTight) {
		if (isVerbose>v0)
			std::cout << "scanning Query genome ...\n";
		std::tie(se2, N2, maxSeqSize) = scanMultiFasta(Q_FN);
	}
	else {
		if (isVerbose>v0)
			std::cout << "Reading Query genome ...\n";
		qGenome = readMultiFasta(Q_FN, 125, true);
	}
	if (isVerbose>v0)
		std::cout << "Reading Reference genome ...\n";
	rGenome = readMultiFasta(R_FN, 123, true);
	
	if (isVerbose>v0)
		std::cout << "Time of I/O = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
	
	int bigRef = 0;
	if ((std::get<0>(rGenome)) >= (1ULL << 32)) {
		bigRef = 1;
		if ((std::get<0>(rGenome)) / k1 >= (1ULL << 32))
			bigRef = 2;
	}

	if (bigRef == 2) {
		if (isVerbose>v0)
			std::cout << "WARNING - LARGE reference file (SIZE / k1 > 4GB), 64-bit arrays\n";
		std::pair<std::uint64_t*, std::uint64_t*> buffer = processRef<std::uint64_t, std::uint64_t>(rGenome, maRushPrime1HashSimplified);
		if (isVerbose>v0)
			std::cout << "Time of processRef = " << stopwatch.stop() << std::endl;
		if (isTight)
			processQueryTight<std::uint64_t, std::uint64_t>(buffer, rGenome, se2, maRushPrime1HashSimplified);
		else
			processQuery<std::uint64_t, std::uint64_t>(buffer, rGenome, qGenome, maRushPrime1HashSimplified);
		if (isVerbose>v0)
			std::cout << "Time of processQuery = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		deleteHashBuffer(buffer);
	}else
	if (bigRef == 1) {
		if (isVerbose>v0)
			std::cout << "WARNING - BIG reference file (>4GB), 64-bit arrays\n";
		std::pair<std::uint64_t*, std::uint32_t*> buffer = processRef<std::uint64_t, std::uint32_t>(rGenome, maRushPrime1HashSimplified);
		if (isVerbose>v0)
			std::cout << "Time of processRef = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		if (isTight)
			processQueryTight<std::uint64_t, std::uint32_t>(buffer, rGenome, se2, maRushPrime1HashSimplified);
		else
			processQuery<std::uint64_t, std::uint32_t>(buffer, rGenome, qGenome, maRushPrime1HashSimplified);
		if (isVerbose>v0)
			std::cout << "Time of processQuery = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		deleteHashBuffer(buffer);
	}
	else {
		std::pair<std::uint32_t*, std::uint32_t*> buffer = processRef<uint32_t, uint32_t>(rGenome, maRushPrime1HashSimplified);
		if (isVerbose>v0)
			std::cout << "Time of processRef = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		if (isTight) 
			processQueryTight<std::uint32_t, std::uint32_t>(buffer, rGenome, se2, maRushPrime1HashSimplified);
		else
			processQuery<std::uint32_t, std::uint32_t>(buffer, rGenome, qGenome, maRushPrime1HashSimplified);
		if (isVerbose>v0)
			std::cout << "Time of processQuery = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		deleteHashBuffer(buffer);
	}
	deleteReading(rGenome);
	if (!isTight)
		deleteReading(qGenome);

	if (isVerbose > v0) {
		std::cout << "Time of deleting = " << stopwatch.stop() << "\nTotal time = " << stopwatch.totalTime() << "\nFINISHED\n";
	}
	return 0;
	
}
