#pragma once

#include "CGnuplot.h"
#include "RCore.h"
#include "IMPES.h"
#include <future>

class IMPES_Parallel : public IMPES {
protected:
	std::vector<std::shared_ptr<std::thread>> threads;
	std::vector<double>s_n1;
	std::vector<double>pc0;
	std::vector<int> wkload;
	int num_of_threads;
	bool isPR1;

public:
	IMPES_Parallel(std::shared_ptr<RCore>& ncore, IMPESArgs ia, bool _isPR1);
	void SolverManifold() override;
	//void parallelPressureSolver(size_t start, size_t end);
	void parallelPressureSolver(std::vector<int> wkload);
	//void SolverPressure(size_t n) override;
	void SolverPressure(int start, int end);
	double SolverSaturation1(double n) override;
	double SolverSaturation2(double n) override;
	void SatIncrementParallel(int& start, int& end, double& n);
	static std::vector<int> Workload(size_t tasks, int num_of_workers = std::thread::hardware_concurrency());
};