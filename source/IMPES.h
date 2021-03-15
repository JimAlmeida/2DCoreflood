#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <assert.h>
#include <map>
#include <exception>
#include <ctime>
#include <thread>

#include "CGnuplot.h"
#include "RCore.h"

struct IMPESArgs {
	double Patm;
	double vazao_inj;
	double uo;
	double uw;
	double DsMax;
	double vpi1;
	double vpi2;
	double si;
	int nbl;
};

class IMPES {
protected:
	bool drenagem; 
	bool embebicao;
	int nbl;
	double uw; 
	double uo; 
	double Pl;
	double Qinj;
	double DSMAX;
	double h;
	double vpi1; 
	double vpi2;
	double t_drenagem; 
	double t_embebicao;
	double tempo_sim;
	double si;

	std::shared_ptr<RCore> core;

	std::vector<std::vector<double>> pc;
	std::vector<std::vector<double>> p;
	std::vector<std::vector<double>> s;
	std::vector<double> g;
	
	std::vector<std::string> tempo_plot;
	std::vector<std::string> tempo_plot_emb;
	std::vector<std::vector<double>> Plot_s;
	std::vector<std::vector<double>> Plot_p;
	std::vector<std::vector<double>> Plot_s_emb;
	std::vector<std::vector<double>> Plot_p_emb;
	std::vector<double> Plot_h;

public:
	IMPES(std::shared_ptr<RCore>& ncore, IMPESArgs ia);
	void Grid();
	virtual void SolverManifold();
	virtual void SolverPressure(size_t n);
	virtual double SolverSaturation1(double n);
	virtual double SolverSaturation2(double n);
	double SatEq1(size_t n, size_t i);
	double SatEq2(size_t n, size_t i);
	std::vector<double> PreEq(size_t n, size_t i);
	std::vector<double> ThomasSolver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, std::vector<double>& d);
	double upstream_weighting(const size_t& n, const size_t& i, const int& orientation);
	double u();
	double dt();
	std::pair<double, double> tempo_inj();
	void PlotSat(bool multiplot = false);
	void PlotPre();

};