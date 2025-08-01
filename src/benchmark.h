#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "alignment.h"
#include "io.h"
#include <string>
#include <vector>
#include <chrono>
#include <map>

/**
 * Estructura para almacenar los resultados de un benchmark
 */
struct BenchmarkResult {
    // Métricas de rendimiento
    double execution_time_ms;      // Tiempo de ejecución en milisegundos
    size_t memory_usage_mb;        // Uso de memoria en MB
    
    // Métricas del alineamiento
    int num_sequences;             // Número de secuencias procesadas
    int original_avg_length;       // Longitud promedio original
    int final_length;              // Longitud final del alineamiento
    int total_gaps;                // Total de gaps insertados
    double gap_percentage;         // Porcentaje de gaps
    
    // Métricas de calidad (si se proporciona referencia)
    double accuracy_score;         // Puntuación de precisión vs referencia
    bool has_reference;            // Si se comparó con referencia
    
    // Metadatos
    std::string dataset_name;      // Nombre del dataset
    std::string timestamp;         // Momento de ejecución
    
    BenchmarkResult() : execution_time_ms(0.0), memory_usage_mb(0), 
                       num_sequences(0), original_avg_length(0), 
                       final_length(0), total_gaps(0), gap_percentage(0.0),
                       accuracy_score(0.0), has_reference(false) {}
};

/**
 * Clase para ejecutar y gestionar benchmarks del alineador MSA
 */
class Benchmark {
public:
    /**
     * Constructor
     */
    Benchmark();
    
    /**
     * Destructor
     */
    ~Benchmark() = default;
    
    /**
     * Ejecuta un benchmark sobre un dataset específico
     * @param dataset_path Ruta al archivo FASTA del dataset
     * @param output_path Ruta donde guardar el resultado (opcional)
     * @return Resultado del benchmark
     */
    BenchmarkResult runSingleBenchmark(const std::string& dataset_path,
                                      const std::string& output_path = "");
    
    /**
     * Ejecuta benchmarks sobre múltiples datasets
     * @param dataset_paths Vector de rutas a los datasets
     * @return Vector de resultados de benchmarks
     */
    std::vector<BenchmarkResult> runMultipleBenchmarks(const std::vector<std::string>& dataset_paths);
    
    /**
     * Compara un alineamiento con una referencia conocida
     * @param alignment_path Ruta al alineamiento generado
     * @param reference_path Ruta al alineamiento de referencia
     * @return Puntuación de precisión (0.0 = mala, 1.0 = perfecta)
     */
    double compareWithReference(const std::string& alignment_path,
                               const std::string& reference_path);
    
    /**
     * Genera un reporte detallado de los benchmarks
     * @param results Vector de resultados
     * @param output_file Archivo donde guardar el reporte
     */
    void generateReport(const std::vector<BenchmarkResult>& results,
                       const std::string& output_file = "");
    
    /**
     * Ejecuta un benchmark de escalabilidad
     * @param base_sequences Secuencias base
     * @param max_sequences Número máximo de secuencias a probar
     * @param step Incremento en cada prueba
     * @return Vector de resultados de escalabilidad
     */
    std::vector<BenchmarkResult> runScalabilityBenchmark(const std::vector<Sequence>& base_sequences,
                                                        int max_sequences = 100,
                                                        int step = 10);
    
    /**
     * Crea datasets sintéticos para pruebas
     * @param num_sequences Número de secuencias
     * @param base_length Longitud base de las secuencias
     * @param mutation_rate Tasa de mutación (0.0 - 1.0)
     * @param output_path Archivo donde guardar el dataset
     */
    void createSyntheticDataset(int num_sequences, int base_length,
                               double mutation_rate, const std::string& output_path);
    
    /**
     * Exporta resultados a formato CSV
     * @param results Vector de resultados
     * @param csv_file Archivo CSV de salida
     */
    void exportToCSV(const std::vector<BenchmarkResult>& results,
                     const std::string& csv_file);

private:
    MSAAligner aligner;
    
    /**
     * Obtiene el uso actual de memoria del proceso
     * @return Uso de memoria en MB
     */
    size_t getCurrentMemoryUsage();
    
    /**
     * Calcula estadísticas básicas de un conjunto de secuencias
     * @param sequences Vector de secuencias
     * @return Mapa con estadísticas (longitud promedio, etc.)
     */
    std::map<std::string, double> calculateSequenceStats(const std::vector<Sequence>& sequences);
    
    /**
     * Genera una secuencia mutada a partir de una base
     * @param base_sequence Secuencia base
     * @param mutation_rate Tasa de mutación
     * @return Secuencia mutada
     */
    std::string mutateSequence(const std::string& base_sequence, double mutation_rate);
    
    /**
     * Calcula la precisión comparando dos alineamientos
     * @param alignment1 Primer alineamiento
     * @param alignment2 Segundo alineamiento
     * @return Puntuación de precisión
     */
    double calculateAlignmentAccuracy(const std::vector<Sequence>& alignment1,
                                     const std::vector<Sequence>& alignment2);
    
    /**
     * Obtiene timestamp actual como string
     * @return Timestamp formateado
     */
    std::string getCurrentTimestamp();
};

#endif // BENCHMARK_H