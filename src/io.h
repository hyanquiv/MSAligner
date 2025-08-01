#ifndef IO_H
#define IO_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

/**
 * Estructura para representar una secuencia biológica
 */
struct Sequence {
    std::string header;     // Encabezado de la secuencia (sin '>')
    std::string sequence;   // Secuencia de nucleótidos/aminoácidos
    
    Sequence() = default;
    Sequence(const std::string& h, const std::string& s) : header(h), sequence(s) {}
};

/**
 * Clase para manejo de entrada/salida de archivos FASTA
 */
class FastaIO {
public:
    /**
     * Lee secuencias de un archivo FASTA
     * @param filename Nombre del archivo FASTA
     * @return Vector de secuencias leídas
     */
    static std::vector<Sequence> readFasta(const std::string& filename);
    
    /**
     * Escribe secuencias a un archivo FASTA
     * @param sequences Vector de secuencias a escribir
     * @param filename Nombre del archivo de salida
     * @param aligned Indica si las secuencias están alineadas (para formato)
     */
    static void writeFasta(const std::vector<Sequence>& sequences, 
                          const std::string& filename, 
                          bool aligned = true);
    
    /**
     * Valida el formato de una secuencia
     * @param sequence Secuencia a validar
     * @return true si es válida, false en caso contrario
     */
    static bool validateSequence(const std::string& sequence);
    
    /**
     * Imprime estadísticas de las secuencias
     * @param sequences Vector de secuencias
     * @param title Título para las estadísticas
     */
    static void printSequenceStats(const std::vector<Sequence>& sequences, 
                                  const std::string& title = "Secuencias");

private:
    /**
     * Limpia una línea removiendo espacios en blanco y caracteres de control
     * @param line Línea a limpiar
     * @return Línea limpia
     */
    static std::string cleanLine(const std::string& line);
};

#endif // IO_H