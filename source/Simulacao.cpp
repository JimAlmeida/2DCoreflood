#include "IMPESParallel.h"
#include <variant>
 

typedef std::chrono::steady_clock cronometro;

void EntradaManual(coreArgs& a, IMPESArgs& b) {
	std::cout << "Entrada manual selecionada!\n";
	std::cout << "Entre com o valor da porosidade: "; std::cin >> a.porosidade; std::cin.get();
	std::cout << "Entre com o valor do diametro: "; std::cin >> a.diametro; std::cin.get();
	std::cout << "Entre com o valor do comprimento: "; std::cin >> a.comprimento; std::cin.get();
	std::cout << "Entre com o valor da permeabilidade: "; std::cin >> a.permeabilidade; std::cin.get();
	std::cout << "Entre com o valor do expoente para a lei de perm. relativa a agua: "; std::cin >> a.exp_krw; std::cin.get();
	std::cout << "Entre com o valor do expoente para a lei de perm. relativa ao oleo: "; std::cin >> a.exp_kro; std::cin.get();
	std::cout << "Entre com o valor do expoente para a lei de pressao capilar: "; std::cin >> a.exp_pc; std::cin.get();
	std::cout << "Entre com o valor máximo de permeabilidade relativa a agua: "; std::cin >> a.krwMax; std::cin.get();
	std::cout << "Entre com o valor máximo de permeabilidade relativa ao oleo: "; std::cin >> a.kroMax; std::cin.get();
	std::cout << "Entre com o valor máximo de pressão capilar: "; std::cin >> a.pcMax; std::cin.get();
	std::cout << "Entre com o valor da saturacao residual de oleo: "; std::cin >> a.sor; std::cin.get();
	std::cout << "Entre com o valor da saturacao irredutível de agua: "; std::cin >> a.swi; std::cin.get();
	std::cout << "Entre com o valor da pressao capilar de threshold: "; std::cin >> a.pct; std::cin.get();

	std::cout << "Entre com o valor da pressao atmosferica: "; std::cin >> b.Patm; std::cin.get();
	std::cout << "Entre com o valor da vazao de injeçao do experimento: "; std::cin >> b.vazao_inj; std::cin.get();
	std::cout << "Entre com o valor da viscosidade do oleo: "; std::cin >> b.uo; std::cin.get();
	std::cout << "Entre com o valor da viscosidade da agua: "; std::cin >> b.uw; std::cin.get();
	std::cout << "Entre com o valor da variacao maxima de saturacao em 1 passo de tempo: "; std::cin >> b.DsMax; std::cin.get();
	std::cout << "Entre com o valor dos volumes porosos de oleo a ser injetado na drenagem: "; std::cin >> b.vpi1; std::cin.get();
	std::cout << "Entre com o valor dos volumes porosos de agua a ser injetado na embebicao: "; std::cin >> b.vpi2; std::cin.get();
	std::cout << "Entre com o valor da saturacao inicial de agua no meio poroso: "; std::cin >> b.si; std::cin.get();
}

void EntradaDisco(coreArgs& a, IMPESArgs& b) {
	std::cout << "Iniciando leitura\n";
	std::ifstream fin("argsDarcy.txt");
	if (fin.good()) {
		for (int i = 0; i < 3; i++) {
			fin.ignore(5000, '\n');
		}
		fin.ignore(17); fin >> a.porosidade; fin.ignore(5000, '\n');
		fin.ignore(17); fin >> a.diametro; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> a.comprimento; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> a.permeabilidade; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> a.exp_krw; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> a.exp_kro; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> a.exp_pc; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> a.krwMax; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> a.kroMax; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> a.pcMax; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> a.sor; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> a.swi; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> a.pct; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> b.Patm; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> b.vazao_inj; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> b.uo; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> b.uw; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> b.DsMax; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> b.vpi1; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> b.vpi2; fin.ignore(5000, '\n');
		fin.ignore(17, '\n'); fin >> b.si; fin.ignore(5000, '\n');
		fin.close();
	}
	else {
		std::cout << "Arquivo de argumentos nao encontrado!\n";
		exit(5);
	}

}

void Teste(coreArgs argsCORE, IMPESArgs argsIMPES) {
	std::cout << "Press enter to initialize performance tests: "; std::cin.get(); std::cout << '\n';
	std::cout << "###################### INITIALIZING PERFORMANCE TESTS #############################\n";

	std::shared_ptr<RCore> tcore1 = std::make_shared<RCore>(argsCORE);
	std::shared_ptr<RCore> tcore2 = std::make_shared<RCore>(argsCORE);
	std::shared_ptr<RCore> tcore3 = std::make_shared<RCore>(argsCORE);

	std::cout << "##################################################################\n";
	std::cout << "Numero de blocos: " << argsIMPES.nbl << '\n';
	std::cout << "##################################################################\n";

	IMPES SR1(tcore1, argsIMPES);
	IMPES_Parallel PR1(tcore2, argsIMPES, true);
	IMPES_Parallel PR2(tcore3, argsIMPES, false);

	SR1.Grid();
	PR1.Grid();
	PR2.Grid();

	auto start = cronometro::now();
	//SR1.SolverManifold();
	auto SR1end = cronometro::now();
	auto PR1start = cronometro::now();
	//PR1.SolverManifold();
	auto PR1end = cronometro::now();
	auto PR2start = cronometro::now();
	PR2.SolverManifold();
	auto PR2end = cronometro::now();

	std::cout << "##################################################################";
	std::cout << "Tempo simulacao serial SR1: " << std::chrono::duration_cast<std::chrono::milliseconds>(SR1end - start).count() << '\n';
	std::cout << "Tempo simulacao paralela PR1: " << std::chrono::duration_cast<std::chrono::milliseconds>(PR1end - PR1start).count() << '\n';
	std::cout << "Tempo simulacao paralela PR2: " << std::chrono::duration_cast<std::chrono::milliseconds>(PR2end - PR2start).count() << '\n';
	std::cout << "##################################################################";
}

int main()
{
	std::cout << "###########################################\n";
	std::cout << "Simulador IMPES" << '\n';
	std::cout << "Simulacao de Reservatorios - LENEP/UENF\n";
	std::cout << "###########################################\n";
	std::cout << "Selecione o modo de entrada de dados:\n";
	std::cout << "1 - Manual\n";
	std::cout << "2 - Por arquivo de disco\n";
	std::cout << "Opcao: "; int op; std::cin >> op; std::cin.get();
	std::cout << "###########################################\n";

	IMPESArgs argsIMPES;
	coreArgs argsCORE;

	switch (op) {
	case 1: EntradaManual(argsCORE, argsIMPES); break;
	case 2: EntradaDisco(argsCORE, argsIMPES); break;
	default: EntradaDisco(argsCORE, argsIMPES); break;
	}
	
	std::cout << "Entre com o número de blocos: "; std::cin >> argsIMPES.nbl; std::cin.get(); std::cout << '\n';

	std::cout << "Selecione o modo de simulacao do experimento:\n";
	std::cout << "1 - Serial\n";
	std::cout << "2 - Paralela para poucos blocos\n";
	std::cout << "3 - Paralela para muitos blocos\n";
	std::cout << "4 - Teste\n";
	std::cout << "Opcao: "; std::cin >> op; std::cin.get();
	std::cout << "###########################################\n";

	std::shared_ptr<RCore> tcore = std::make_shared<RCore>(argsCORE);
	IMPES* sim = nullptr;
	switch (op) {
		case 1: {
			sim = new IMPES(tcore, argsIMPES);
			break;
		}
		case 2: {
			sim = new IMPES_Parallel(tcore, argsIMPES, true);
			break;
		};
		case 3: {
			sim = new IMPES_Parallel(tcore, argsIMPES, false);
			break;
		}
		case 4: {
			Teste(argsCORE, argsIMPES);
			break;
		}
		default: {
			sim = new IMPES_Parallel(tcore, argsIMPES, true);
			break;
		}
	}
	if (sim != nullptr) {
		sim->Grid();
		sim->SolverManifold();
		sim->PlotSat();
		sim->PlotPre();
	}
	delete sim;
}