// CopMEM.cpp : Defines the entry point for the console application.


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


typedef std::pair<std::string, size_t> SequenceItem;
typedef std::vector<SequenceItem> Sequences;
typedef std::tuple<size_t, char*, Sequences> Reading;

typedef std::uint32_t myuint;

typedef std::uint64_t myuint2;


Reading readMultiFasta(std::string fn, const size_t paddingSize, const char paddingChar, bool removeNs);


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
bool isVerbose = false;
//////////////////// GLOBALS ////////////////////////////


void displayHelp(const char* progName) {
	std::cout << "Usage: " << progName << " [-l n] [-v] <-o MEMs_file> <Ref_genome> <Query_genome>\n";
	std::cout << "Attention: -o is a required parameter. l is optional (default: 100).\n";
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
			isVerbose = true;
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


template <class MyUINT> void genCumm(Reading& r, MyUINT* cumm, std::uint32_t(*hashFunc)(const char*)) {
	const size_t MULTI1 = 128;
	const size_t k1MULTI1 = k1 * MULTI1;

	size_t N = std::get<0>(r);
	char* gen = std::get<1>(r);

	std::fill(cumm, cumm + HASH_SIZE + 2, 0);

	uint32_t hashPositions[MULTI1];
	size_t i;

	for (i = 0; i < N - K - k1MULTI1; i += k1MULTI1) {
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
	for ( ; i < N - K + 1; i += k1) {
		uint32_t h = hashFunc(gen + i) + 2;
		++cumm[h];
	}
	//////////////////// processing the end part of R //////////////////////

	std::partial_sum(cumm, cumm + HASH_SIZE + 1, cumm);
}


void dumpMEM(SequenceItem& item1, SequenceItem& item2, size_t* match, std::string &s) {
	size_t baseindex1 = match[0] - item1.second;
	size_t baseindex2 = match[1] - item2.second;
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


void postProcess(std::vector<size_t*> &matches, Reading& t1, Reading& t2, size_t index2) {
	Sequences& seq2 = std::get<2>(t2);
	SequenceItem seq2item = seq2[index2];
	f_matches << "> " << seq2item.first << std::endl;

	if (matches.size() == 0)
		return;

	char* gen1 = std::get<1>(t1);
	Sequences& seq1 = std::get<2>(t1);
	if(matches.size() > 1)
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
			if(prev_match[0] != match[0])
				seq1item = findSeqDesc(match[0], seq1);
			dumpMEM(seq1item, seq2item, match, dumpString);
			++count;
		}
		prev_match = match;
	}
	if (isVerbose)
		if (count == 1)
			std::cout << seq2item.first << ": " << count << " match.\n";
		else
			std::cout << seq2item.first << ": " << count << " matches.\n";
	for (auto match: matches)
		delete[] match;
	matches.clear();
}


template<class MyUINT>
MyUINT* processRef(Reading& t1, std::uint32_t(*hashFunc)(const char*)) {
	CStopWatch stopwatch;
	stopwatch.start();
	size_t N1 = std::get<0>(t1);
	const size_t hashCount = (N1 - K + 1) / k1;
	const myuint MULTI2 = 128;
	const myuint k1MULTI2 = k1 * MULTI2;

	std::cout << "Hash count = " << hashCount << std::endl;

	char* gen1 = std::get<1>(t1);

	MyUINT* buf = new MyUINT[hashCount + HASH_SIZE + 2];
	MyUINT* sampledPositions = buf;
	MyUINT* cumm = buf + hashCount;

	std::cout << "\tprocessRef: init = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
	genCumm(t1, cumm, hashFunc);
	std::cout << "\tprocessRef: cumm(1) = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
	uint32_t hashPositions[MULTI2];
	MyUINT i1;

	for (i1 = 0; i1 < N1 - K - k1MULTI2; i1 += k1MULTI2) {
		char* tempPointer = gen1 + i1;
		for (myuint temp = 0; temp < MULTI2; ++temp) {
			hashPositions[temp] = hashFunc(tempPointer) + 1;
			tempPointer += k1;
			_prefetch((char*)(cumm + hashPositions[temp]), 1);
		}

		for (myuint temp = 0; temp < MULTI2; ++temp) {
			_prefetch((char*)(sampledPositions + *(cumm + hashPositions[temp])), 1);
		}

		MyUINT i2 = i1;
		for (myuint temp = 0; temp < MULTI2; ++temp) {
			sampledPositions[cumm[hashPositions[temp]]] = i2;
			++cumm[hashPositions[temp]];
			i2 += k1;
		}
	}

	//////////////////// processing the end part of R
	for ( ; i1 < N1 - K + 1; i1 += k1) {
		char* tempPointer = gen1 + i1;
		uint32_t h = hashFunc(gen1 + i1) + 1;
		sampledPositions[cumm[h]] = i1;
		++cumm[h];
	}
	std::cout << "\tprocessRef: cumm(2) = " << stopwatch.stop() << std::endl;
	return buf;
}


template<class MyUINT>
void processQuery(MyUINT* buffer, Reading& t1, Reading& t2, std::uint32_t(*hashFunc)(const char*)) {
	const int L_PLUS_ONE = L + 1;
	const int L2 = L / 2;
	const int LK2 = (L - K) / 2;
	const int LK2_MINUS_2 = LK2 - 2;
	const int LK2_MINUS_4 = LK2 - 4;
	const int K_PLUS_LK22 = K + LK2_MINUS_2;
	const int K_PLUS_LK24 = K + LK2_MINUS_4;

	size_t N1 = std::get<0>(t1);

	char* gen1 = std::get<1>(t1);
	Sequences seq1 = std::get<2>(t1);
	const size_t hashCount = (N1 - K + 1) / k1;
	const myuint MULTI = 256;
	const myuint k2MULTI = k2 * MULTI;

	char* curr1 = gen1;

	size_t N2 = std::get<0>(t2);
	char* gen2 = std::get<1>(t2);
	Sequences seq2 = std::get<2>(t2);

	char* curr2 = gen2;

	MyUINT* sampledPositions = buffer;
	MyUINT* cumm = buffer + hashCount;

	std::vector<size_t*> matches;

	myuint left2 = 0, right2 = 0;


	uint32_t hArray[MULTI];
	MyUINT posArray[MULTI * 2];

	size_t currSeqIndex = 0;
	size_t nextSeqBeginPos = (seq2.size() > 1 ? seq2[currSeqIndex + 1].second : SIZE_MAX);

	f_matches.open(matchesFN);

	myuint l1, l2, r1, r2;

	size_t i1;
	for (i1 = 0; i1 < N2 - K + 1 - k2MULTI; i1 += k2MULTI) {
		size_t i1temp = i1;
		char* prev = gen2 + i1;
		curr2 = gen2 + i1;
		size_t tempCount = 0;
		for (size_t i2 = 0; i2 < MULTI; ++i2) {
			if (*curr2 != 'N') {
				hArray[tempCount++] = hashFunc(curr2);
			}
			curr2 += k2;
		}
		for (size_t i2 = 0; i2 < tempCount; ++i2) {
			memcpy(posArray + i2 * 2, cumm + hArray[i2], sizeof(MyUINT) * 2);
		}

		curr2 = gen2 + i1;	// !!!
		for (size_t i2 = 0; i2 < tempCount; ++i2) {
			if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
				curr2 += k2;
				continue;
			}

			memcpy(&l2, curr2 - LK2, sizeof(myuint));
			memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(myuint));

			for (MyUINT j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
				curr1 = gen1 + sampledPositions[j];

				memcpy(&l1, curr1 - LK2, sizeof(myuint));
				memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(myuint));

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
						tempMatch[0] = p1 + 1 - gen1;
						tempMatch[1] = p2 + 1 - gen2;
						tempMatch[2] = right - p1 - 1;
						while ((size_t)(p2 + 1 - gen2) >= nextSeqBeginPos) {
							if (currSeqIndex + 1 < seq2.size()) {
								postProcess(matches, t1, t2, currSeqIndex);
								++currSeqIndex;
								nextSeqBeginPos = seq2[currSeqIndex+1].second;
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
	for ( ; i1 < N2 - K + 1; i1 += k2) {
		char* prev = gen2 + i1;
		curr2 = gen2 + i1;
		size_t tempCount = 0;

		if (*curr2 != 'N') {
			memcpy(posArray + tempCount * 2, cumm + hashFunc(curr2), sizeof(MyUINT) * 2);
			++tempCount;
		}

		curr2 = gen2 + i1;
		for (size_t i2 = 0; i2 < tempCount; ++i2) {
			if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
				curr2 += k2;
				continue;
			}

			memcpy(&l2, curr2 - LK2, sizeof(myuint));
			memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(myuint));

			for (MyUINT j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
				curr1 = gen1 + sampledPositions[j];

				memcpy(&l1, curr1 - LK2, sizeof(myuint));
				memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(myuint));

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
						tempMatch[0] = p1 + 1 - gen1;
						tempMatch[1] = p2 + 1 - gen2;
						tempMatch[2] = right - p1 - 1;

						while ((size_t)(p2 + 1 - gen2) > nextSeqBeginPos) {
							if (currSeqIndex + 1 < seq2.size()) {
								postProcess(matches, t1, t2, currSeqIndex);
								++currSeqIndex;
								nextSeqBeginPos = seq2[currSeqIndex+1].second;
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
			postProcess(matches, t1, t2, currSeqIndex);
			++currSeqIndex;
			nextSeqBeginPos = seq2[currSeqIndex + 1].second;
	}
	postProcess(matches, t1, t2, currSeqIndex);
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


Reading readMultiFasta(std::string fn, const size_t paddingSize, const char paddingChar, bool removeNs) {
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
	char* buf1 = new char[N + 2 * paddingSize];
	if (buf1 == nullptr) {
		std::cerr << "\nFile '" << fn << "' is too large. Quit.";
		exit(1);
	}
	memset(buf1, paddingChar, paddingSize);
	f.read(buf1 + paddingSize, N);
	std::cout << "\treadMultiFasta: Reading file from disk " << stopWatch.stop() << "\n";
	stopWatch.resume();
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

	std::cout << "\treadMultiFasta: Analysing file " << stopWatch.stop() << "\n";
	std::cout << "\treadMultiFasta: Genome data size " << (dst - gen) << std::endl;
	std::cout << "\treadMultiFasta: Sequences " << seq.size() << std::endl;

	return { (dst - gen), gen, seq };
}


int main(int argc, char* argv[]) {
	std::cout.setf(std::ios_base::unitbuf);
	if (argc < 3) {
		displayHelp(argv[0]);
		return 1;
	}
	processCmd(argc, argv);
	displayParams();

	CStopWatch stopwatch;
	stopwatch.start();
	std::cout << "Reading Reference genome ...\n";
	Reading gen1Data = readMultiFasta(R_FN, L + 1, 123, true);
	std::cout << "Reading Query genome ...\n";
	Reading gen2Data = readMultiFasta(Q_FN, L + 1, 125, true);
	std::cout << "Time of I/O = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
	bool isFileBig = ((std::get<0>(gen1Data)) >= (1ULL << 32)) || ((std::get<0>(gen2Data)) >= (1ULL << 32));
	if (isFileBig) {
		std::cout << "WARNING - LARGE files (>4GB), 64-bit arrays\n";
		std::uint64_t* buffer = processRef<uint64_t>(gen1Data, maRushPrime1HashSimplified);
		std::cout << "Time of processRef = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		processQuery(buffer, gen1Data, gen2Data, maRushPrime1HashSimplified);
		std::cout << "Time of processQuery = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		delete[] (std::get<1>(gen1Data) - L - 1);
		delete[] (std::get<1>(gen2Data) - L - 1);
		delete[] buffer;
	}
	else {
		std::uint32_t* buffer = processRef<uint32_t>(gen1Data, maRushPrime1HashSimplified);
		std::cout << "Time of processRef = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		processQuery(buffer, gen1Data, gen2Data, maRushPrime1HashSimplified);
		std::cout << "Time of processQuery = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		delete[] (std::get<1>(gen1Data) - L - 1);
		delete[] (std::get<1>(gen2Data) - L - 1);
		delete[] buffer;
	}
	std::cout << "Time of deleting = " << stopwatch.stop() << std::endl;
	std::cout << "Total time = " << stopwatch.totalTime() << std::endl;
	std::cout << "FINISHED\n";
	return 0;
}
