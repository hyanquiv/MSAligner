#include <iostream>
#include <chrono>
#include <iomanip>
#include "io.h"
#include "alignment.h"

void printUsage(const char* program_name) {
    std::cout << "\n🧬 Alineador Múltiple de Secuencias (MSA)\n" << std::endl;
    std::cout << "Uso: " << program_name << " <archivo_entrada.fasta> <archivo_salida.fasta>" << std::endl;
    std::cout << "\nDescripción:" << std::endl;
    std::cout << "  Este programa realiza alineamiento múltiple de secuencias usando:" << std::endl;
    std::cout << "  1. Matriz de distancias basada en identidad porcentual" << std::endl;
    std::cout << "  2. Construcción de árbol guía con algoritmo UPGMA" << std::endl;
    std::cout << "  3. Alineamiento progresivo con programación dinámica" << std::endl;
    std::cout << "\nEjemplo:" << std::endl;
    std::cout << "  " << program_name << " sequences.fasta aligned_sequences.fasta" << std::endl;
    std::cout << "\nFormato de entrada:" << std::endl;
    std::cout << "  - Archivo FASTA estándar con múltiples secuencias" << std::endl;
    std::cout << "  - Mínimo 2 secuencias requeridas" << std::endl;
    std::cout << "  - Soporta secuencias de ADN y proteínas" << std::endl;
    std::cout << std::endl;
}

void printHeader() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "🧬 ALINEADOR MÚLTIPLE DE SECUENCIAS (MSA) v1.0" << std::endl;
    std::cout << "   Implementación en C++ con algoritmo progresivo" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
}

void printSummary(const std::chrono::duration<double>& duration, 
                 const std::map<std::string, int>& stats,
                 int num_sequences) {
    std::cout << "\n" << std::string(50, '─') << std::endl;
    std::cout << "📊 RESUMEN DEL ALINEAMIENTO" << std::endl;
    std::cout << std::string(50, '─') << std::endl;
    std::cout << "⏱️  Tiempo total: " << std::fixed << std::setprecision(3) 
              << duration.count() << " segundos" << std::endl;
    std::cout << "🔢 Secuencias procesadas: " << num_sequences << std::endl;
    std::cout << "📏 Longitud final: " << stats.at("final_length") << " posiciones" << std::endl;
    std::cout << "🕳️  Gaps insertados: " << stats.at("total_gaps") << std::endl;
    
    if (stats.at("final_length") > 0) {
        double gap_percentage = (static_cast<double>(stats.at("total_gaps")) / 
                               (num_sequences * stats.at("final_length"))) * 100.0;
        std::cout << "📈 Porcentaje de gaps: " << std::fixed << std::setprecision(1) 
                  << gap_percentage << "%" << std::endl;
    }
    
    std::cout << std::string(50, '─') << std::endl;
    std::cout << "✅ Alineamiento completado exitosamente!" << std::endl;
}

bool validateInputFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "❌ Error: No se puede abrir el archivo de entrada: " << filename << std::endl;
        return false;
    }
    
    // Verificar que el archivo no esté vacío
    file.seekg(0, std::ios::end);
    if (file.tellg() == 0) {
        std::cerr << "❌ Error: El archivo de entrada está vacío: " << filename << std::endl;
        return false;
    }
    
    file.close();
    return true;
}

bool validateOutputPath(const std::string& filename) {
    std::ofstream test_file(filename);
    if (!test_file.is_open()) {
        std::cerr << "❌ Error: No se puede escribir en el archivo de salida: " << filename << std::endl;
        return false;
    }
    test_file.close();
    
    // Remover el archivo de prueba
    std::remove(filename.c_str());
    return true;
}

int main(int argc, char* argv[]) {
    printHeader();
    
    // Verificar argumentos de línea de comandos
    if (argc != 3) {
        printUsage(argv[0]);
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    
    // Validar archivos
    if (!validateInputFile(input_file)) {
        return 1;
    }
    
    if (!validateOutputPath(output_file)) {
        return 1;
    }
    
    try {
        // Medir tiempo de ejecución
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Paso 1: Leer secuencias del archivo FASTA
        std::cout << "\n📁 Leyendo archivo de entrada: " << input_file << std::endl;
        auto sequences = FastaIO::readFasta(input_file);
        
        if (sequences.empty()) {
            std::cerr << "❌ Error: No se pudieron leer secuencias del archivo." << std::endl;
            return 1;
        }
        
        if (sequences.size() < 2) {
            std::cerr << "❌ Error: Se necesitan al menos 2 secuencias para el alineamiento." << std::endl;
            return 1;
        }
        
        // Mostrar estadísticas de entrada
        FastaIO::printSequenceStats(sequences, "Secuencias de entrada");
        
        // Paso 2: Crear alineador y realizar MSA
        MSAAligner aligner;
        std::cout << "\n🚀 Iniciando proceso de alineamiento..." << std::endl;
        
        auto aligned_sequences = aligner.alignSequences(sequences);
        
        if (aligned_sequences.empty()) {
            std::cerr << "❌ Error: Fallo en el proceso de alineamiento." << std::endl;
            return 1;
        }
        
        // Paso 3: Mostrar árbol guía
        aligner.printGuideTree();
        
        // Paso 4: Guardar secuencias alineadas
        std::cout << "\n💾 Guardando secuencias alineadas en: " << output_file << std::endl;
        FastaIO::writeFasta(aligned_sequences, output_file, true);
        
        // Mostrar estadísticas de salida
        FastaIO::printSequenceStats(aligned_sequences, "Secuencias alineadas");
        
        // Calcular tiempo total
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
        
        // Mostrar resumen final
        auto stats = aligner.getAlignmentStats();
        printSummary(duration, stats, static_cast<int>(sequences.size()));
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ Error inesperado: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "\n❌ Error desconocido durante la ejecución." << std::endl;
        return 1;
    }
}