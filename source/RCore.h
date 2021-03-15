#pragma once
#include <string>
#include <vector>
#include <exception>
#include <iostream>
#include <fstream>
#include <math.h>
#include "CGnuplot.h"

struct coreArgs {
	double porosidade;
	double diametro;
	double comprimento;
	double permeabilidade;
	double exp_krw;
	double exp_kro;
	double exp_pc;
	double krwMax;
	double kroMax;
	double pcMax;
	double sor;
	double swi;
	double pct;
};

class RCore {
protected:
	double a; 
	double b; 
	double c;
	double porosity;
	double diameter;
	double length;
	double permeability;
	double krwMAX;
	double kroMAX;
	double pcMAX;
	double swi;
	double sor;
	double Pct;
	static const double pi;
	std::string wettability;

public:
	RCore(coreArgs cr, std::string _w= "W");
	inline void setC(double _c) { c = _c; }
	double volPoroso();
	double Area();
	double krw(double s);
	double kro(double s);
	double lambda_w(double s, double uw);
	double lambda_o(double s, double uo);
	double lambda_t(double s, double uw, double uo);
	double pc(double s);
	double fw(double s, double uw, double uo);

	void PlotKr();
	void PlotPc();

	friend class IMPES;
	friend class IMPES_Parallel;
	friend class IMPES_ParallelCV;
};

