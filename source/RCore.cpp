#include "RCore.h"

const double RCore::pi = 3.14159265358979323846;

RCore::RCore(coreArgs cr, std::string _w) : porosity(cr.porosidade),
	diameter(cr.diametro), length(cr.comprimento), permeability(cr.permeabilidade), wettability(_w), swi(cr.swi), sor(cr.sor), a(cr.exp_krw), b(cr.exp_kro),
	c(cr.exp_pc), krwMAX(cr.krwMax), kroMAX(cr.kroMax), pcMAX(cr.pcMax), Pct(cr.pct) {
}
double RCore::volPoroso() { return 0.25 * RCore::pi * length * diameter * diameter * porosity; }
double RCore::Area() { return 0.25 * RCore::pi * pow(diameter, 2); }
double RCore::krw(double s) {
	double r;
	if (s >= 1 - sor) { r = krwMAX; }
	else if (s <= swi) { r = 0; }
	else {
		r = krwMAX * pow(((s - swi) / (1 - sor - swi)), a);
	}
	return r;
}
double RCore::kro(double s) {
	double r;
	if (s >= 1 - sor) { r = 0; }
	else if (s <= swi) { r = kroMAX; }
	else {
		r = kroMAX * pow(((1 - sor - s) / (1 - sor - swi)), a);
	}
	return r;
}
double RCore::lambda_w(double s, double uw) { return krw(s) / uw; }
double RCore::lambda_o(double s, double uo) { return kro(s) / uo; }
double RCore::lambda_t(double s, double uw, double uo) { return lambda_w(s, uw) + lambda_o(s, uo); }
double RCore::pc(double s) {
	double r;
	if (s >= 1 - sor) { r = Pct; }
	else {
		r = pcMAX * pow(((1 - sor - s) / (1 - sor - swi)), c) + Pct;
	}
	return r;
}
double RCore::fw(double s, double uw, double uo) { return lambda_w(s, uw) / lambda_t(s, uw, uo); }
void RCore::PlotKr() {
	try {
		CGnuplot plot;
		std::ofstream tmpfile("kr.txt");
		std::vector<double> s;
		for (double k = swi; k < 1 - sor; k += 0.01) {
			s.push_back(k);
		}
		for (int i = 0; i < s.size(); i++) {
			tmpfile << s[i] << "	" << krw(s[i]) << "	" << kro(s[i]) << '\n';
		}
		tmpfile.close();
		plot << "set title \"Permeabilidade Relativa \"";
		plot << "set xlabel \"Sw\"";
		plot << "set ylabel \"Kr\"";
		//plot << "plot 'krw.txt' title \"krw\" with lines";
		plot << "plot 'kr.txt' using 1:2 title \"krw\" with lines, 'kr.txt' using 1:3 title \"kro\" with lines";
		//plot << "set legend \"Series1\"";
		system("Pause");
		system("erase kr.txt");
	}
	catch (std::exception& e) {
		std::cerr << e.what() << '\n';
	}
}
void RCore::PlotPc() {
	try {
		CGnuplot plot;
		std::ofstream tmpfile("pc.txt");
		std::vector<double> s;
		for (double k = swi; k < 1; k += 0.01) {
			s.push_back(k);
		}
		for (int i = 0; i < s.size(); i++) {
			tmpfile << s[i] << "	" << pc(s[i]) << '\n';
		}
		tmpfile.close();
		plot << "set title \"Pressao Capilar \"";
		plot << "set xlabel \"Sw\"";
		plot << "set ylabel \"Pc (atm)\"";
		plot << "set yrange [0:6]";
		//plot << "plot 'krw.txt' title \"krw\" with lines";
		plot << "plot 'pc.txt' using 1:2 title \"pc\" with lines";
		//plot << "set legend \"Series1\"";
		system("Pause");
		system("erase pc.txt");
	}
	catch (std::exception& e) {
		std::cerr << e.what() << '\n';
	}
}