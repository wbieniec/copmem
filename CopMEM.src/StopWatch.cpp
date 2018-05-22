#include "StopWatch.h"


CStopWatch::CStopWatch(){
	running = 0;
}


CStopWatch::~CStopWatch(){
}


void CStopWatch::start(){
	running = 1;
	elapsed.clear();
	t1 = std::chrono::steady_clock::now();
}


double CStopWatch::stop(){
	t2  = std::chrono::steady_clock::now();
	running = 0;
	double el = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0;
	elapsed.push_back(el);
	return el;
}


void CStopWatch::resume(){
	running = 1;
	t1 = std::chrono::steady_clock::now();
}


double CStopWatch::totalTime(){
	double sum = 0;
	for (std::vector<double>::iterator it = elapsed.begin(); it != elapsed.end(); ++it) {
		sum += *it;
	}
	return sum;
}
