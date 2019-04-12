/*
Symbols:
Note:
 The algorithm finds the same points several times. How this can be avoided?
 -Why to disthinguis beetween std::vector and boost::vector how about just use boost vector
  for everything

  TODO: Put string name to the system so apper in the text 
  TODO: Homogenize the arguments in functions like, initial condition, x
  TODO: The algorithm can be made more eficently if at each domain partition, we dont explore the corners
  TODO: The algorithm seem to give slightly diferents graphs if we cahnge the values of M and N at start


*/
#include "al21managerErandi.h"
#include "util.h"
#include <chrono>
#include <ctime>
#include <iostream>

int main() {
	auto start = std::chrono::system_clock::now();
	runAl21NoDelay();
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	LogAndStdout lcout("points_log.txt");
	lcout << "finished computation at " << std::ctime(&end_time)
		<< "elapsed time: " << elapsed_seconds.count() << "s\n";
	return 0;
}