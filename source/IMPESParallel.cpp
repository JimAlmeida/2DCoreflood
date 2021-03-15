#include "IMPESParallel.h"

IMPES_Parallel::IMPES_Parallel(std::shared_ptr<RCore>& ncore, IMPESArgs ia, bool _isPR1) : IMPES(ncore, ia), threads(std::thread::hardware_concurrency(), nullptr), s_n1(nbl), pc0(nbl) {
	num_of_threads = std::thread::hardware_concurrency();
	wkload = Workload(nbl, num_of_threads);
	isPR1 = _isPR1;
};
void IMPES_Parallel::SolverManifold() {
	//Calculando a saturação S0 inicial de todos os blocos;
	std::vector<int> dre;
	std::vector<int> emb;
	std::vector<double> s0(nbl, si);
	std::vector<double> pc0(nbl, core->Pct);
	s.push_back(s0);
	Plot_s.push_back(s0);
	tempo_plot.push_back(std::to_string(0));
	pc.push_back(pc0); 
	size_t n = 1;
	drenagem = true;

	double t; double t_fim = 0;
	for (double t_sim = 0; t_sim < t_drenagem; t_sim += t) {
		if (isPR1) {
			t = SolverSaturation1(n);
		}
		else {
			t = SolverSaturation2(n);
		}
		n += 1;
		if (n % 2000 == 0) {
			dre.push_back(n-1);
			std::cout << "Terminou-se o passo de tempo " << n << " com Dt = " << t << '\n';
			std::cout << "Tempo da simulação total: " << t_sim << " s\n";
		}
		t_fim = t_sim;
	}
	dre.push_back(n-1);
	drenagem = false;
	embebicao = true;
	if (embebicao) {
		for (double t_sim = t_drenagem; t_sim < t_embebicao; t_sim += t) {
			if (isPR1) {
				t = SolverSaturation1(n);
			}
			else {
				t = SolverSaturation2(n);
			}
			n += 1;
			if (n % 100 == 0) {
				emb.push_back(n-1);
				std::cout << '\n';
				std::cout << "Terminou-se o passo de tempo " << n << " com Dt = " << t << '\n';
				std::cout << "Tempo da simulação total: " << t_sim << " s\n";
			}
			t_fim = t_sim;
		}
	}
	emb.push_back(n-1);
	embebicao = false;
	
	std::cout << "Resolvendo o campo da pressao...\n";
	parallelPressureSolver(Workload(s.size(), std::thread::hardware_concurrency()));
	std::cout << "Simulacao completa!\n";
	
	for(auto& pos:dre){
		Plot_s.push_back(s[pos]);
		Plot_p.push_back(p[pos]);
	}
	
	for(auto& pos:emb){
		Plot_s_emb.push_back(s[pos]);
		Plot_p_emb.push_back(p[pos]);
	}
}
void IMPES_Parallel::parallelPressureSolver(std::vector<int> wkload) {
	p.resize(s.size());
	std::vector<std::shared_ptr<std::thread>> threads;
	for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
		threads.push_back(std::make_shared<std::thread>(&IMPES_Parallel::SolverPressure, this, wkload[i], wkload[i+1]));
	}
	for (auto& t : threads) {
		t->join();
	}
}
void IMPES_Parallel::SolverPressure(int start, int end) {
	for (int n = start; n < end; n++) {
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
		p[n] = ThomasSolver(a, b, c, d);
	}
}
double IMPES_Parallel::SolverSaturation1(double n) {
	std::vector<double>s_n1(nbl);
	std::vector<double>pc0(nbl);
	g.clear(); g.resize(nbl);
	for (int i = 0; i < nbl; i++) {
		SatEq2(n, i);
	}
	double t = dt();
	
	for (int i = 0; i < nbl; ++i) {
		s_n1[i] = g[i] * t + s[n - 1][i];
		pc0[i] = core->pc(s_n1[i]);
	}
	s.push_back(s_n1);
	pc.push_back(pc0);
	return t;
}
double IMPES_Parallel::SolverSaturation2(double n) {
	auto s1 = std::chrono::high_resolution_clock::now();
	g.clear(); g.resize(nbl);
	auto s2 = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < num_of_threads; i++) {
		threads[i] = std::make_shared<std::thread>(std::thread(&IMPES_Parallel::SatIncrementParallel, this, std::ref(wkload[i]), std::ref(wkload[i + 1]), std::ref(n)));
	}
	for (auto& t : threads) {
		t->join();
	}
	auto s3 = std::chrono::high_resolution_clock::now();
	double time_step = dt();
	auto s4 = std::chrono::high_resolution_clock::now();
	for (int j = 0; j < nbl; j++) {
		s_n1[j] = g[j] * time_step + s[n - 1][j];
		pc0[j] = core->pc(s_n1[j]);
	}
	auto s5 = std::chrono::high_resolution_clock::now();
	s.push_back(s_n1);
	pc.push_back(pc0);
	auto s6 = std::chrono::high_resolution_clock::now();
	
	std::cout << "###################################################\n";
	std::cout << "G Resizing" << std::chrono::duration_cast<std::chrono::nanoseconds>(s2-s1).count() << '\n';
	std::cout << "Threaded G calc" << std::chrono::duration_cast<std::chrono::nanoseconds>(s3-s2).count() << '\n';
	std::cout << "Dt" << std::chrono::duration_cast<std::chrono::nanoseconds>(s4-s3).count() << '\n';
	std::cout << "Time increment" << std::chrono::duration_cast<std::chrono::nanoseconds>(s5-s4).count() << '\n';
	std::cout << "Pushing back" << std::chrono::duration_cast<std::chrono::nanoseconds>(s6-s5).count() << '\n';
	std::cout << "###################################################\n";
	return time_step;
}
void IMPES_Parallel::SatIncrementParallel(int& start, int& end, double& n) {
	for (int i = start; i < end; i++) {
		SatEq2(n, i);
	}
}
std::vector<int> IMPES_Parallel::Workload(size_t tasks, int num_of_workers)  {
	int remainder = tasks % num_of_workers;
	int equal_tasks = (tasks - remainder) / num_of_workers;

	std::vector<int> tasks_for_each_thread(num_of_workers);
	int last = 0;
	for (int j = 0; j < tasks_for_each_thread.size(); j++) {
		tasks_for_each_thread[j] = last + equal_tasks;
		last += equal_tasks;
	}
	for (int i = 0; i < remainder; i++) {
		tasks_for_each_thread[i] += 1;
		for (int k = i + 1; k < tasks_for_each_thread.size(); k++) {
			tasks_for_each_thread[k] += 1;
		}
	}
	tasks_for_each_thread.emplace(tasks_for_each_thread.begin(), 0);
	return tasks_for_each_thread;
}