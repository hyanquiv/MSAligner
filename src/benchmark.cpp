#include "benchmark.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <random>
#include <algorithm>
#include <iomanip>
#include <cstdlib>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#elif __linux__
#include <fstream>
#include <unistd.h>
#elif __APPLE__
#include <mach/mach.h>
#endif

Benchmark::Benchmark() {
    // Constructor vacío, el alineador se inicializa por defecto
}

BenchmarkResult Benchmark::runSingleBenchmark(const std::string& dataset_path,
                                             const std::string& output_path) {
    BenchmarkResult result;
    result.dataset_name = dataset_path;
    result.timestamp = getCurrentTimestamp();
    
    try {
        // Leer secuencias del dataset
        std::vector<Sequence> sequences = FastaIO::readFasta(dataset_path);
        result.num_sequences = sequences.size();
        
        if (sequences.empty()) {
            std::cerr << "Error: No se pudieron leer secuencias de " << dataset_path << std::endl;
            return result;
        }
        
        // Calcular estadísticas originales
        auto stats = calculateSequenceStats(sequences);
        result.original_avg_length = static_cast<int>(stats["avg_length"]);
        
        // Medir memoria inicial
        size_t initial_memory = getCurrentMemoryUsage();
        
        // Medir tiempo de ejecución
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Ejecutar alineamiento
        std::vector<Sequence> aligned_sequences = aligner.alignSequences(sequences);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        
        // Calcular tiempo transcurrido
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        result.execution_time_ms = duration.count() / 1000.0;
        
        // Medir memoria final
        size_t final_memory = getCurrentMemoryUsage();
        result.memory_usage_mb = final_memory > initial_memory ? final_memory - initial_memory : final_memory;
        
        // Obtener estadísticas del alineamiento
        auto alignment_stats = aligner.getAlignmentStats();
        result.final_length = alignment_stats["final_length"];
        result.total_gaps = alignment_stats["total_gaps"];
        
        // Calcular porcentaje de gaps
        if (result.final_length > 0 && result.num_sequences > 0) {
            int total_positions = result.final_length * result.num_sequences;
            result.gap_percentage = (static_cast<double>(result.total_gaps) / total_positions) * 100.0;
        }
        
        // Guardar resultado si se especifica ruta
        if (!output_path.empty()) {
            FastaIO::writeFasta(aligned_sequences, output_path);
        }
        
        std::cout << "Benchmark completado para " << dataset_path << std::endl;
        std::cout << "  Tiempo: " << result.execution_time_ms << " ms" << std::endl;
        std::cout << "  Memoria: " << result.memory_usage_mb << " MB" << std::endl;
        std::cout << "  Secuencias: " << result.num_sequences << std::endl;
        std::cout << "  Gaps: " << result.gap_percentage << "%" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error en benchmark: " << e.what() << std::endl;
    }
    
    return result;
}

std::vector<BenchmarkResult> Benchmark::runMultipleBenchmarks(const std::vector<std::string>& dataset_paths) {
    std::vector<BenchmarkResult> results;
    
    std::cout << "Ejecutando " << dataset_paths.size() << " benchmarks..." << std::endl;
    
    for (size_t i = 0; i < dataset_paths.size(); ++i) {
        std::cout << "\\nBenchmark " << (i + 1) << "/" << dataset_paths.size() << ": " << dataset_paths[i] << std::endl;
        
        BenchmarkResult result = runSingleBenchmark(dataset_paths[i]);
        results.push_back(result);
    }
    
    return results;
}

double Benchmark::compareWithReference(const std::string& alignment_path,
                                     const std::string& reference_path) {
    try {
        std::vector<Sequence> alignment = FastaIO::readFasta(alignment_path);
        std::vector<Sequence> reference = FastaIO::readFasta(reference_path);
        
        return calculateAlignmentAccuracy(alignment, reference);
    } catch (const std::exception& e) {
        std::cerr << "Error comparando con referencia: " << e.what() << std::endl;
        return 0.0;
    }
}

void Benchmark::generateReport(const std::vector<BenchmarkResult>& results,
                              const std::string& output_file) {
    std::ostream* out = &std::cout;
    std::ofstream file;
    
    if (!output_file.empty()) {
        file.open(output_file);
        out = &file;
    }
    
    *out << "\\n" << std::string(80, '=') << std::endl;
    *out << "REPORTE DE BENCHMARKS - MSA Aligner" << std::endl;
    *out << std::string(80, '=') << std::endl;
    
    if (results.empty()) {
        *out << "No hay resultados para mostrar." << std::endl;
        return;
    }
    
    // Estadísticas generales
    double total_time = 0.0;
    size_t total_memory = 0;
    int total_sequences = 0;
    
    for (const auto& result : results) {
        total_time += result.execution_time_ms;
        total_memory += result.memory_usage_mb;
        total_sequences += result.num_sequences;
    }
    
    *out << "\\nRESUMEN GENERAL:" << std::endl;
    *out << "  Total de benchmarks: " << results.size() << std::endl;
    *out << "  Tiempo total: " << total_time << " ms" << std::endl;
    *out << "  Tiempo promedio: " << (total_time / results.size()) << " ms" << std::endl;
    *out << "  Memoria promedio: " << (total_memory / results.size()) << " MB" << std::endl;
    *out << "  Total secuencias procesadas: " << total_sequences << std::endl;
    
    // Resultados detallados
    *out << "\\nRESULTADOS DETALLADOS:" << std::endl;
    *out << std::string(80, '-') << std::endl;
    
    for (const auto& result : results) {
        *out << "Dataset: " << result.dataset_name << std::endl;
        *out << "  Timestamp: " << result.timestamp << std::endl;
        *out << "  Secuencias: " << result.num_sequences << std::endl;
        *out << "  Longitud original promedio: " << result.original_avg_length << std::endl;
        *out << "  Longitud final: " << result.final_length << std::endl;
        *out << "  Tiempo de ejecución: " << result.execution_time_ms << " ms" << std::endl;
        *out << "  Uso de memoria: " << result.memory_usage_mb << " MB" << std::endl;
        *out << "  Total de gaps: " << result.total_gaps << std::endl;
        *out << "  Porcentaje de gaps: " << std::fixed << std::setprecision(2) << result.gap_percentage << "%" << std::endl;
        
        if (result.has_reference) {
            *out << "  Precisión vs referencia: " << std::fixed << std::setprecision(3) << result.accuracy_score << std::endl;
        }
        
        *out << std::string(40, '-') << std::endl;
    }
    
    if (file.is_open()) {
        std::cout << "Reporte guardado en: " << output_file << std::endl;
    }
}

std::vector<BenchmarkResult> Benchmark::runScalabilityBenchmark(const std::vector<Sequence>& base_sequences,
                                                               int max_sequences,
                                                               int step) {
    std::vector<BenchmarkResult> results;
    
    std::cout << "Ejecutando benchmark de escalabilidad..." << std::endl;
    std::cout << "Desde " << step << " hasta " << max_sequences << " secuencias (step: " << step << ")" << std::endl;
    
    for (int n = step; n <= max_sequences && n <= static_cast<int>(base_sequences.size()); n += step) {
        std::cout << "\\nProbando con " << n << " secuencias..." << std::endl;
        
        // Crear subset de secuencias
        std::vector<Sequence> subset(base_sequences.begin(), base_sequences.begin() + n);
        
        // Crear archivo temporal
        std::string temp_file = "temp_scalability_" + std::to_string(n) + ".fasta";
        FastaIO::writeFasta(subset, temp_file);
        
        // Ejecutar benchmark
        BenchmarkResult result = runSingleBenchmark(temp_file);
        result.dataset_name = "Scalability_" + std::to_string(n) + "_sequences";
        
        results.push_back(result);
        
        // Limpiar archivo temporal
        std::remove(temp_file.c_str());
    }
    
    return results;
}

void Benchmark::createSyntheticDataset(int num_sequences, int base_length,
                                      double mutation_rate, const std::string& output_path) {
    std::vector<Sequence> synthetic_sequences;
    
    // Generar secuencia base aleatoria
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 3);
    
    std::string bases = "ATCG";
    std::string base_sequence;
    
    for (int i = 0; i < base_length; ++i) {
        base_sequence += bases[dis(gen)];
    }
    
    // Crear secuencias mutadas
    for (int i = 0; i < num_sequences; ++i) {
        Sequence seq;
        seq.header = "Synthetic_Seq_" + std::to_string(i + 1) + " | Generated with mutation rate " + std::to_string(mutation_rate);
        seq.sequence = mutateSequence(base_sequence, mutation_rate);
        
        synthetic_sequences.push_back(seq);
    }
    
    // Guardar dataset
    FastaIO::writeFasta(synthetic_sequences, output_path);
    
    std::cout << "Dataset sintético creado: " << output_path << std::endl;
    std::cout << "  Secuencias: " << num_sequences << std::endl;
    std::cout << "  Longitud base: " << base_length << std::endl;
    std::cout << "  Tasa de mutación: " << mutation_rate << std::endl;
}

void Benchmark::exportToCSV(const std::vector<BenchmarkResult>& results,
                           const std::string& csv_file) {
    std::ofstream file(csv_file);
    
    if (!file.is_open()) {
        std::cerr << "Error: No se pudo crear el archivo CSV " << csv_file << std::endl;
        return;
    }
    
    // Header CSV
    file << "Dataset,Timestamp,NumSequences,OriginalAvgLength,FinalLength,";
    file << "ExecutionTime_ms,MemoryUsage_MB,TotalGaps,GapPercentage,";
    file << "AccuracyScore,HasReference\\n";
    
    // Datos
    for (const auto& result : results) {
        file << result.dataset_name << ",";
        file << result.timestamp << ",";
        file << result.num_sequences << ",";
        file << result.original_avg_length << ",";
        file << result.final_length << ",";
        file << result.execution_time_ms << ",";
        file << result.memory_usage_mb << ",";
        file << result.total_gaps << ",";
        file << result.gap_percentage << ",";
        file << result.accuracy_score << ",";
        file << (result.has_reference ? "true" : "false") << "\\n";
    }
    
    file.close();
    std::cout << "Resultados exportados a CSV: " << csv_file << std::endl;
}

// Métodos privados

size_t Benchmark::getCurrentMemoryUsage() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    return pmc.WorkingSetSize / (1024 * 1024); // MB
#elif __linux__
    std::ifstream status_file("/proc/self/status");
    std::string line;
    while (std::getline(status_file, line)) {
        if (line.substr(0, 6) == "VmRSS:") {
            std::istringstream iss(line);
            std::string key, value, unit;
            iss >> key >> value >> unit;
            return std::stoul(value) / 1024; // Convert KB to MB
        }
    }
    return 0;
#elif __APPLE__
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) == KERN_SUCCESS) {
        return info.resident_size / (1024 * 1024); // MB
    }
    return 0;
#else
    return 0; // Fallback para otros sistemas
#endif
}

std::map<std::string, double> Benchmark::calculateSequenceStats(const std::vector<Sequence>& sequences) {
    std::map<std::string, double> stats;
    
    if (sequences.empty()) {
        return stats;
    }
    
    double total_length = 0.0;
    int min_length = sequences[0].sequence.length();
    int max_length = sequences[0].sequence.length();
    
    for (const auto& seq : sequences) {
        int length = seq.sequence.length();
        total_length += length;
        min_length = std::min(min_length, length);
        max_length = std::max(max_length, length);
    }
    
    stats["avg_length"] = total_length / sequences.size();
    stats["min_length"] = min_length;
    stats["max_length"] = max_length;
    
    return stats;
}

std::string Benchmark::mutateSequence(const std::string& base_sequence, double mutation_rate) {
    std::string mutated = base_sequence;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> mutation_prob(0.0, 1.0);
    std::uniform_int_distribution<> base_dis(0, 3);
    
    std::string bases = "ATCG";
    
    for (size_t i = 0; i < mutated.length(); ++i) {
        if (mutation_prob(gen) < mutation_rate) {
            char new_base;
            do {
                new_base = bases[base_dis(gen)];
            } while (new_base == mutated[i]); // Asegurar que sea una mutación real
            
            mutated[i] = new_base;
        }
    }
    
    return mutated;
}

double Benchmark::calculateAlignmentAccuracy(const std::vector<Sequence>& alignment1,
                                            const std::vector<Sequence>& alignment2) {
    if (alignment1.size() != alignment2.size()) {
        return 0.0;
    }
    
    int total_positions = 0;
    int matching_positions = 0;
    
    for (size_t i = 0; i < alignment1.size(); ++i) {
        const std::string& seq1 = alignment1[i].sequence;
        const std::string& seq2 = alignment2[i].sequence;
        
        int min_length = std::min(seq1.length(), seq2.length());
        total_positions += min_length;
        
        for (int j = 0; j < min_length; ++j) {
            if (seq1[j] == seq2[j]) {
                matching_positions++;
            }
        }
    }
    
    return total_positions > 0 ? static_cast<double>(matching_positions) / total_positions : 0.0;
}

std::string Benchmark::getCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t time = std::chrono::system_clock::to_time_t(now);
    
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}