#include "IMPES.h"

IMPES::IMPES(std::shared_ptr<RCore>& ncore, IMPESArgs ia)
	: core(ncore), uw(ia.uw), uo(ia.uo), Pl(ia.Patm), Qinj(ia.vazao_inj), DSMAX(ia.DsMax), nbl(ia.nbl), vpi1(ia.vpi1), vpi2(ia.vpi2), si(ia.si) {
	std::pair<double, double> tempos = tempo_inj();
	t_drenagem = tempos.first;
	t_embebicao = tempos.second;
}
void IMPES::Grid() {
	h = core->length / nbl;
	for (double k = h; k < nbl; k += h) {
		Plot_h.push_back(k - h / 2);
	}
	std::vector<double> aux(nbl, 0);
	g = aux;
}
void IMPES::SolverManifold() {
	clock_t time_drenagem = clock();
	try {
		//Calculando a saturação S0 inicial de todos os blocos;
		std::vector<double> s0(nbl, si);
		std::vector<double> pc0(nbl, core->Pct);
		s.push_back(s0);
		Plot_s.push_back(s0);
		tempo_plot.push_back(std::to_string(0));
		pc.push_back(pc0);
		//PROBLEMA AQUI 
		SolverPressure(0);
		size_t n = 1;
		drenagem = true;
		double t; double t_fim = 0;
		for (double t_sim = 0; t_sim < t_drenagem; t_sim += t) {
			t = SolverSaturation1(n);
			for (int j = 0; j < nbl; j++) {
				pc0[j] = core->pc(s[n][j]);
			}
			pc.push_back(pc0);
			SolverPressure(n);
			n += 1;
			if (n < 3000 && n % 500 == 0) {
				Plot_p.push_back(p[n - 1]);
				Plot_s.push_back(s[n - 1]);
				tempo_plot.push_back(std::to_string(t_fim));
				std::cout << "Terminou-se o passo de tempo " << n << " com Dt = " << t << '\n';
				std::cout << "Tempo da simulação total: " << t_sim << " s\n";
			}
			else if (n > 3000 && n % 2000 == 0) {
				Plot_p.push_back(p[n - 1]);
				Plot_s.push_back(s[n - 1]);
				tempo_plot.push_back(std::to_string(t_fim));
				std::cout << "Terminou-se o passo de tempo " << n << " com Dt = " << t << '\n';
				std::cout << "Tempo da simulação total: " << t_sim << " s\n";
			}
			t_fim = t_sim;
		}

		drenagem = false;
		std::cout << "TEMPO FINAL DRENAGEM: " << (clock() - time_drenagem) / CLOCKS_PER_SEC << '\n';
		std::vector<double> semb(nbl, core->swi);
		s.push_back(semb);
		for (int j = 0; j < nbl; j++) {
			pc0[j] = core->pcMAX;
		}
		pc.push_back(pc0);
		SolverPressure(n);
		n += 1;
		Plot_s.push_back(s[n - 1]);
		Plot_p.push_back(p[n - 1]);
		Plot_s_emb.push_back(s[n - 1]);
		Plot_p_emb.push_back(p[n - 1]);
		tempo_plot_emb.push_back(std::to_string(t_drenagem));
		tempo_plot.push_back(std::to_string(t_drenagem));
		std::cout << '\n';
		std::cout << "Terminou-se o passo de tempo " << n << " com Dt = " << t << '\n';
		std::cout << "Tempo da simulação total: " << t_fim << " s\n";

		embebicao = true;
		if (embebicao) {
			for (double t_sim = t_drenagem; t_sim < t_embebicao; t_sim += t) {
				t = SolverSaturation1(n);
				for (int j = 0; j < nbl; j++) {
					pc0[j] = core->pc(s[n][j]);
				}
				pc.push_back(pc0);
				SolverPressure(n);
				n += 1;
				if (n % 100 == 0) {
					Plot_s_emb.push_back(s[n - 1]);
					Plot_p_emb.push_back(p[n - 1]);
					tempo_plot_emb.push_back(std::to_string(t_sim));
					std::cout << '\n';
					std::cout << "Terminou-se o passo de tempo " << n << " com Dt = " << t << '\n';
					std::cout << "Tempo da simulação total: " << t_sim << " s\n";
				}
				t_fim = t_sim;
			}
		}
		embebicao = false;
		std::cout << "TEMPO FINAL EMBEBICAO" << '\n';

		Plot_s_emb.push_back(s[n - 1]);
		Plot_p_emb.push_back(p[n - 1]);
		tempo_plot_emb.push_back(std::to_string(t_embebicao));
	}
	catch (std::exception& e) {
		std::cerr << "Erro no método IMPES::SolverManifold() - " << e.what();
	}
}
void IMPES::SolverPressure(size_t n) {
	std::vector<double> a(nbl, 0.0);
	std::vector<double> b(nbl, 0.0);
	std::vector<double> c(nbl, 0.0);
	std::vector<double> d(nbl, 0.0);
	std::vector<double> coefficients;
	for (int i = 0; i < nbl; i++) {
		coefficients = PreEq(n, i);
		if (i > 0) {
			a[i] = coefficients[0];
		}
		if (i < nbl - 1) {
			c[i] = coefficients[2];
		}
		b[i] = coefficients[1];
		d[i] = coefficients[3];
	}
	
	p.push_back(ThomasSolver(a, b, c, d));

	
}
double IMPES::SolverSaturation1(double n) {
	std::vector<double>s_n1(nbl);
	g.clear(); g.resize(nbl);
	for (int i = 0; i < nbl; i++) {
		SatEq2(n, i);
	}
	double t = dt();
	try {
		for (int i = 0; i < nbl; ++i) {
			s_n1[i] = g[i] * t + s[n - 1][i];
		}
		s.push_back(s_n1);
	}
	catch (std::exception& e) {
		std::cerr << "Erro no método IMPES::SolverSaturation" << e.what() << '\n';
	}
	return t;
}
double IMPES::SolverSaturation2(double n) {
	std::vector<double>s_n1(nbl);
	double t = dt();
	g.clear(); g.resize(nbl);
	double gi;
	for (int i = 0; i < nbl; i++) {
		gi = SatEq1(n, i);
		s_n1[i] = gi * t + s[n - 1][i];
	}
	s.push_back(s_n1);
	return dt();
}
double IMPES::SatEq1(size_t n, size_t i) {
	double gi;
	if (i == 0) {
		double Sp = upstream_weighting(n, i, 1);
		double Ao_POS = core->lambda_o(Sp, uo);
		if (drenagem) {
			gi = (-core->permeability / (core->porosity * pow(h, 2))) * (Ao_POS * (p[n - 1][i + 1] - p[n - 1][i]) + (u() * h / core->permeability));
		}
		else if (embebicao) {
			gi = (-core->permeability / (core->porosity * pow(h, 2))) * (Ao_POS * (p[n - 1][i + 1] - p[n - 1][i]));
		}
	}
	else if (i == nbl - 1) {
		double Sp = upstream_weighting(n, i, 1);
		double Sn = upstream_weighting(n, i, -1);
		double Ao_POS = core->lambda_o(Sp, uo);
		double Ao_NEG = core->lambda_o(Sn, uo);
		gi = (-4 * core->permeability / (3 * core->porosity * h * h)) * (2 * Ao_POS * (Pl - p[n - 1][i]) - Ao_NEG * (p[n - 1][i] - p[n - 1][i - 1]));

	}
	else {
		double Sp = upstream_weighting(n, i, 1);
		double Sn = upstream_weighting(n, i, -1);
		double Ao_POS = core->lambda_o(Sp, uo);
		double Ao_NEG = core->lambda_o(Sn, uo);
		gi = (-core->permeability / (core->porosity * pow(h, 2))) * (Ao_POS * (p[n - 1][i + 1] - p[n - 1][i]) - Ao_NEG * (p[n - 1][i] - p[n - 1][i - 1]));
	}
	g[i] = gi;
	return gi;
}
double IMPES::SatEq2(size_t n, size_t i) {
	double gi = 0;
	if (i == 0) {
		double ts = upstream_weighting(n, i, 1);
		double fwPOS = core->fw(ts, uw, uo);
		double lamb_oPOS = core->lambda_o(ts, uo);
		double T = fwPOS * lamb_oPOS;
		if (drenagem) {
			gi = -1 / core->porosity * (core->permeability / pow(h, 2) * (T * (pc[n - 1][i + 1] - pc[n - 1][i])) + u() / h * fwPOS);
		}
		else if (embebicao) {
			gi = -1 / core->porosity * (core->permeability / pow(h, 2) * (T * (pc[n - 1][i + 1] - pc[n - 1][i])) + u() / h * (fwPOS - 1));
		}
	}
	else if (i == nbl - 1) {
		double tsp = upstream_weighting(n, i, 1);
		double tsn = upstream_weighting(n, i, -1);
		double fwPOS = core->fw(tsp, uw, uo);
		double fwNEG = core->fw(tsn, uw, uo);
		double lamb_oPOS = core->lambda_o(tsp, uo);
		double lamb_oNEG = core->lambda_o(tsn, uo);
		double T = fwNEG * lamb_oNEG;
		gi = -1 / core->porosity * ((2 * u() / h) * (fwPOS - fwNEG) + ((4 * core->permeability) / (3 * pow(h, 2))) * ((2 * fwPOS * lamb_oPOS * (Pl - pc[n - 1][i])) - (T * (pc[n - 1][i] - pc[n - 1][i - 1]))));

	}
	else {
		double tsp = upstream_weighting(n, i, 1);
		double tsn = upstream_weighting(n, i, -1);
		double fwPOS = core->fw(tsp, uw, uo);
		double fwNEG = core->fw(tsn, uw, uo);
		double lamb_oPOS = core->lambda_o(tsp, uo);
		double lamb_oNEG = core->lambda_o(tsn, uo);
		gi = -1 / core->porosity * (core->permeability / pow(h, 2) * (lamb_oPOS * fwPOS * pc[n - 1][i + 1] - (lamb_oNEG * fwNEG + lamb_oPOS * fwPOS) * pc[n - 1][i] + lamb_oNEG * fwNEG * pc[n - 1][i - 1]) + u() / h * (fwPOS - fwNEG));
	}
	g[i] = gi;
	return gi;
}
std::vector<double> IMPES::PreEq(size_t n, size_t i) {
	try {
		std::vector<double> c(4, 0);
		//c[0] p_i-1
		//c[1] p_i
		//c[2] p_i+1
		//c[3] lado direito da eq;
		if (i == 0) {
			double tsp = upstream_weighting(n, i, 1);
			double lambda_wPOS = core->lambda_w(tsp, uw);
			double lambda_tPOS = core->lambda_t(tsp, uw, uo);
			c[0] = 0;
			c[1] = -lambda_tPOS;
			c[2] = lambda_tPOS;
			c[3] = ((-u() * h / core->permeability) + lambda_wPOS * (pc[n][i + 1] - pc[n][i]));
		}
		else if (i == nbl - 1) {
			double tsp = upstream_weighting(n, i, 1);
			double tsn = upstream_weighting(n, i, -1);
			double lambda_wPOS = core->lambda_w(tsp, uw);
			double lambda_wNEG = core->lambda_w(tsn, uw);
			double lambda_tPOS = core->lambda_t(tsp, uw, uo);
			double lambda_tNEG = core->lambda_t(tsn, uw, uo);

			c[0] = lambda_tNEG;
			c[1] = -(2 * lambda_tPOS + lambda_tNEG);
			c[2] = 0;
			c[3] = -2 * lambda_tPOS * Pl - (2 * lambda_wPOS + lambda_wNEG) * pc[n][i] + lambda_wNEG * pc[n][i - 1];
		}
		else {
			double tsp = upstream_weighting(n, i, 1);
			double tsn = upstream_weighting(n, i, -1);
			double lambda_wPOS = core->lambda_w(tsp, uw);
			double lambda_wNEG = core->lambda_w(tsn, uw);
			double lambda_tPOS = core->lambda_t(tsp, uw, uo);
			double lambda_tNEG = core->lambda_t(tsn, uw, uo);
			c[0] = lambda_tNEG;
			c[1] = -(lambda_tNEG + lambda_tPOS);
			c[2] = lambda_tPOS;
			c[3] = (lambda_wPOS * pc[n][i + 1] - (lambda_wPOS + lambda_wNEG) * pc[n][i] + lambda_wNEG * pc[n][i - 1]);
		}
		return c;
	}
	catch (std::exception& e) {
		std::cerr << "Erro no método IMPES::PreEq() - " << e.what();
	}
}
std::vector<double> IMPES::ThomasSolver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, std::vector<double>& d) {
	size_t N = d.size();

	c[0] = c[0] / b[0];
	d[0] = d[0] / b[0];

	double m;
	for (int i = 1; i < N; i++) {
		m = 1.0 / (b[i] - a[i] * c[i - 1]);
		c[i] = c[i] * m;
		d[i] = (d[i] - a[i] * d[i - 1]) * m;
	}
	for (int i = N - 1; i-- > 0;) {
		d[i] = d[i] - (c[i] * d[i + 1]);
	}
	return d;
}
double IMPES::upstream_weighting(const size_t& n, const size_t& i, const int& orientation) {
	if (n == 0) {
		return s[n][i];
	}
	if (orientation == 1) {
		return s[n - 1][i];
	}
	else if (orientation == -1) {
		return s[n - 1][i - 1];
	}
}
double IMPES::u() {
	return Qinj / (core->Area());
};
double IMPES::dt() {
	assert(g.size() != 0);
	std::vector<double> absG;
	for (double& e : g) {
		absG.push_back(abs(e));
	}
	double t = DSMAX / *std::max_element(absG.begin(), absG.end());
	if (isinf(t)) {
		return DSMAX;
	}
	if (drenagem) {
		if (tempo_sim + t > t_drenagem) {
			return t_drenagem - tempo_sim;
		}
		else if (t > 20) { return t / 100; }
		else return t;
	}
	else if (embebicao) {
		if (tempo_sim + t > t_embebicao) {
			return t_embebicao - tempo_sim;
		}
		else return t;
	}
}
std::pair<double, double> IMPES::tempo_inj() {
	double param1 = core->volPoroso() / Qinj;
	return std::make_pair<double, double>(vpi1 * param1, vpi2 * param1);
}
void IMPES::PlotSat(bool multiplot) {
	//PLOT DRENAGEM
	try {
		CGnuplot plot; std::string title;
		int NUM_PLOTS = Plot_s.size();
		std::cout << "Número total de plots: " << NUM_PLOTS << '\n';

		std::ofstream tmpfile("sat.txt");
		for (int i = 0; i < nbl; ++i) {
			tmpfile << Plot_h[i] << "\t";
			for (int j = 0; j < NUM_PLOTS; ++j) {
				tmpfile << Plot_s[j][i] << "\t";
			}
			tmpfile << '\n';
		}
		tmpfile.close();

		plot << "set xlabel \"X\"";
		plot << "set ylabel \"Sw\"";

		if (multiplot) {
			int plots = Plot_s.size() / 6;
			plots += 1;
			size_t j = 2;
			for (int p = 1; p < plots; p++) {
				plot << "set multiplot layout 3,2 title \"Distribuicao de Saturacao no Meio Poroso - DRENAGEM\"";
				while (j < 2 + (p * 6)) {
					title = tempo_plot[j - 2];
					//plot_command += ",'sat.txt' using 1:" + std::to_string(j) + "title \"Tempo"+std::to_string(j)+"\" with lines";
					plot << "set xlabel \"X(t = " + title + "s)\"";
					plot << "plot 'sat.txt' using 1:" + std::to_string(j) + " w filledcu above y1=0 lc rgb \"blue\" notitle, 'sat.txt' using 1:" + std::to_string(j) + " w filledcu below y1=1 lc rgb \"black\" notitle";
					j += 1;
				}
				plot << "unset multiplot";
				std::cout << "J - WHILE: " << j << '\n';
				system("Pause");
			}
			std::cout << "J: " << j << '\n';
			plot << "set multiplot layout 3,2 title \"Distribuicao de Saturacao no Meio Poroso - DRENAGEM\"";
			for (int k = j; k <= Plot_s.size() + 1; ++k) {
				//plot_command += ",'sat.txt' using 1:" + std::to_string(j) + "title \"Tempo"+std::to_string(j)+"\" with lines";
				plot << "set xlabel \"X(t = " + title + "s)\"";
				plot << "plot 'sat.txt' using 1:" + std::to_string(k) + " w filledcu above y1=0 lc rgb \"blue\" notitle, 'sat.txt' using 1:" + std::to_string(j) + " w filledcu below y1=1 lc rgb \"black\" notitle";
			}
			plot << "unset multiplot";
			system("Pause");
		}
		else {
			plot << "set title \"Distribuicao de Saturacao no Meio Poroso - DRENAGEM\"";
			plot << "set yrange [0:1]";
			plot << "plot for[i=2:" + std::to_string(Plot_s.size()) + "] 'sat.txt' using 1:i notitle with lines";
			system("Pause");
		}
		//system("erase sat.txt");
	}
	catch (std::exception& e) {
		std::cerr << e.what() << '\n';
	}
	//PLOT EMBEBICAO
	try {
		CGnuplot emb_plot; std::string title;
		int NUM_PLOTS = Plot_s_emb.size();
		std::cout << "Número total de plots: " << NUM_PLOTS << '\n';
		std::ofstream tmpfile("sat_emb.txt");
		for (int i = 0; i < nbl; ++i) {
			tmpfile << Plot_h[i] << "\t";
			for (int j = 0; j < NUM_PLOTS; ++j) {
				tmpfile << Plot_s_emb[j][i] << "\t";
			}
			tmpfile << '\n';
		}
		tmpfile.close();
		emb_plot << "set xlabel \"X\"";
		emb_plot << "set ylabel \"So\"";
		if (multiplot) {
			int plots = Plot_s_emb.size() / 6;
			plots += 1;
			int j = 2;
			for (int p = 1; p < plots; p++) {
				emb_plot << "set multiplot layout 3,2 title \"Distribuicao de Saturacao no Meio Poroso - EMBEBICAO\"";
				while (j < 2 + (p * 6)) {
					//plot_command += ",'sat.txt' using 1:" + std::to_string(j) + "title \"Tempo"+std::to_string(j)+"\" with lines";
					title = tempo_plot_emb[j - 2];
					emb_plot << "set xlabel \"X(t = " + title + "s)\"";
					emb_plot << "plot 'sat_emb.txt' using 1:" + std::to_string(j) + " w filledcu above y1=0 lc rgb \"blue\" notitle, 'sat_emb.txt' using 1:" + std::to_string(j) + " w filledcu below y1=1 lc rgb \"black\" notitle";
					j += 1;
				}
				emb_plot << "unset multiplot";
				std::cout << "J - WHILE: " << j << '\n';
				system("Pause");
			}
			std::cout << "J: " << j << '\n';
			emb_plot << "set multiplot layout 3,2 title \"Distribuicao de Saturacao no Meio Poroso - DRENAGEM\"";
			for (int k = j; k <= Plot_s_emb.size() + 1; ++k) {
				//plot_command += ",'sat.txt' using 1:" + std::to_string(j) + "title \"Tempo"+std::to_string(j)+"\" with lines";
				title = tempo_plot_emb[k - 2];
				emb_plot << "set xlabel \"X(t = " + title + "s)\"";
				emb_plot << "plot 'sat_emb.txt' using 1:" + std::to_string(k) + " w filledcu above y1=0 lc rgb \"blue\" notitle, 'sat_emb.txt' using 1:" + std::to_string(j) + " w filledcu below y1=1 lc rgb \"black\" notitle";
			}
			emb_plot << "unset multiplot";
			system("Pause");
		}
		else {
			emb_plot << "set title \"Distribuicao de Saturacao no Meio Poroso - EMBEBICAO\"";
			emb_plot << "set yrange [0:1]";
			emb_plot << "plot for[i=2:" + std::to_string(Plot_s_emb.size()) + "] 'sat_emb.txt' using 1:i notitle with lines";
			system("Pause");
		}
		//system("erase sat.txt");
	}
	catch (std::exception& e) {
		std::cerr << e.what() << '\n';
	}

}
void IMPES::PlotPre() {
	//PLOT DRENAGEM
	try {
		CGnuplot plot;
		int NUM_PLOTS = Plot_p.size();
		std::cout << "Número total de plots: " << NUM_PLOTS << '\n';
		std::ofstream tmpfile("pre_dre.txt");
		for (int i = 0; i < nbl; ++i) {
			tmpfile << Plot_h[i] << "\t";
			for (int j = 0; j < NUM_PLOTS; ++j) {
				tmpfile << Plot_p[j][i] << "\t";
			}
			tmpfile << '\n';
		}
		tmpfile.close();
		plot << "set xlabel \"X\"";
		plot << "set ylabel \"P(atm)\"";
		plot << "set title \"Distribuicao de Pressao no Meio Poroso - DRENAGEM\"";
		plot << "plot for[i=2:" + std::to_string(Plot_p.size() + 1) + "] 'pre_dre.txt' using 1:i title 'T'.(i-1) with lines";
		system("Pause");
		//system("erase sat.txt");
	}
	catch (std::exception& e) {
		std::cerr << e.what() << '\n';
	}
	//PLOT EMBEBICAO
	try {
		CGnuplot emb_plot;
		int NUM_PLOTS = Plot_p_emb.size();
		std::ofstream tmpfile("pre_emb.txt");
		for (int i = 0; i < nbl; ++i) {
			tmpfile << Plot_h[i] << "\t";
			for (int j = 0; j < NUM_PLOTS; ++j) {
				tmpfile << Plot_p_emb[j][i] << "\t";
			}
			tmpfile << '\n';
		}
		tmpfile.close();
		emb_plot << "set xlabel \"X\"";
		emb_plot << "set ylabel \"P(atm)\"";
		emb_plot << "set title \"Distribuicao de Pressao no Meio Poroso - EMBEBICAO\"";
		emb_plot << "plot for[i=2:" + std::to_string(Plot_p_emb.size()) + "] 'pre_emb.txt' using 1:i title 'T'.(i-1) with lines";
		system("Pause");
		//system("erase sat.txt");
	}
	catch (std::exception& e) {
		std::cerr << e.what() << '\n';
	}

}