#include "benchmark.h"
#include <iostream>
#include <vector>
#include <string>

/**
 * Programa principal para ejecutar benchmarks del MSA Aligner
 */
int main(int argc, char* argv[]) {
    std::cout << "============================================================" << std::endl;
    std::cout << "MSA ALIGNER - SISTEMA DE BENCHMARKS v1.0" << std::endl;
    std::cout << "============================================================" << std::endl;
    
    if (argc < 2) {
        std::cout << std::endl;
        std::cout << "Uso: " << argv[0] << " <comando> [opciones]" << std::endl;
        std::cout << std::endl;
        std::cout << "Comandos disponibles:" << std::endl;
        std::cout << "  single <dataset.fasta> [output.fasta]  - Ejecutar benchmark individual" << std::endl;
        std::cout << "  multiple <dataset1> <dataset2> ...     - Ejecutar múltiples benchmarks" << std::endl;
        std::cout << "  scalability <dataset.fasta> [max] [step] - Test de escalabilidad" << std::endl;
        std::cout << "  synthetic <num_seq> <length> <mut_rate> <output.fasta> - Crear dataset sintético" << std::endl;
        std::cout << std::endl;
        std::cout << "Ejemplos:" << std::endl;
        std::cout << "  " << argv[0] << " single benchmarks/datasets/small/dna_sample.fasta" << std::endl;
        std::cout << "  " << argv[0] << " scalability entrada.fasta 50 10" << std::endl;
        std::cout << "  " << argv[0] << " synthetic 20 100 0.1 synthetic_test.fasta" << std::endl;
        std::cout << std::endl;
        return 1;
    }
    
    std::string command = argv[1];
    Benchmark benchmark;
    
    try {
        if (command == "single") {
            if (argc < 3) {
                std::cerr << "Error: Falta especificar el dataset" << std::endl;
                return 1;
            }
            
            std::string dataset = argv[2];
            std::string output = (argc > 3) ? argv[3] : "";
            
            std::cout << "Ejecutando benchmark individual..." << std::endl;
            BenchmarkResult result = benchmark.runSingleBenchmark(dataset, output);
            
            std::vector<BenchmarkResult> results = {result};
            benchmark.generateReport(results);
            
        } else if (command == "multiple") {
            if (argc < 3) {
                std::cerr << "Error: Faltan especificar los datasets" << std::endl;
                return 1;
            }
            
            std::vector<std::string> datasets;
            for (int i = 2; i < argc; ++i) {
                datasets.push_back(argv[i]);
            }
            
            std::cout << "Ejecutando benchmarks múltiples..." << std::endl;
            std::vector<BenchmarkResult> results = benchmark.runMultipleBenchmarks(datasets);
            
            benchmark.generateReport(results, "benchmarks/results/multiple_benchmark_report.txt");
            benchmark.exportToCSV(results, "benchmarks/results/multiple_benchmark_results.csv");
            
        } else if (command == "scalability") {
            if (argc < 3) {
                std::cerr << "Error: Falta especificar el dataset base" << std::endl;
                return 1;
            }
            
            std::string dataset = argv[2];
            int max_sequences = (argc > 3) ? std::stoi(argv[3]) : 50;
            int step = (argc > 4) ? std::stoi(argv[4]) : 10;
            
            // Leer secuencias base
            std::vector<Sequence> base_sequences = FastaIO::readFasta(dataset);
            
            if (base_sequences.empty()) {
                std::cerr << "Error: No se pudieron leer las secuencias del dataset base" << std::endl;
                return 1;
            }
            
            std::cout << "Ejecutando test de escalabilidad..." << std::endl;
            std::vector<BenchmarkResult> results = benchmark.runScalabilityBenchmark(base_sequences, max_sequences, step);
            
            benchmark.generateReport(results, "benchmarks/results/scalability_report.txt");
            benchmark.exportToCSV(results, "benchmarks/results/scalability_results.csv");
            
        } else if (command == "synthetic") {
            if (argc < 6) {
                std::cerr << "Error: Parámetros insuficientes para dataset sintético" << std::endl;
                std::cerr << "Uso: " << argv[0] << " synthetic <num_seq> <length> <mut_rate> <output.fasta>" << std::endl;
                return 1;
            }
            
            int num_sequences = std::stoi(argv[2]);
            int base_length = std::stoi(argv[3]);
            double mutation_rate = std::stod(argv[4]);
            std::string output_path = argv[5];
            
            std::cout << "Creando dataset sintético..." << std::endl;
            benchmark.createSyntheticDataset(num_sequences, base_length, mutation_rate, output_path);
            
        } else {
            std::cerr << "Error: Comando desconocido '" << command << "'" << std::endl;
            std::cerr << "Comandos válidos: single, multiple, scalability, synthetic" << std::endl;
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << std::endl << "Benchmark completado exitosamente!" << std::endl;
    return 0;
}