#include "io.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <iomanip> // Agregar este encabezado para std::setprecision
        
std::vector<Sequence> FastaIO::readFasta(const std::string& filename) {
    std::vector<Sequence> sequences;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo " << filename << std::endl;
        return sequences;
    }
    
    std::string line;
    std::string current_header;
    std::string current_sequence;
    bool in_sequence = false;
    
    while (std::getline(file, line)) {
        line = cleanLine(line);
        
        if (line.empty()) {
            continue;
        }
        
        if (line[0] == '>') {
            // Si ya teníamos una secuencia anterior, la guardamos
            if (in_sequence && !current_header.empty()) {
                if (validateSequence(current_sequence)) {
                    sequences.emplace_back(current_header, current_sequence);
                } else {
                    std::cerr << "Advertencia: Secuencia inválida ignorada: " 
                              << current_header << std::endl;
                }
            }
            
            // Nueva secuencia
            current_header = line.substr(1); // Quitar el '>'
            current_sequence.clear();
            in_sequence = true;
            
        } else if (in_sequence) {
            current_sequence += line;
        }
    }
    
    // Guardar la última secuencia
    if (in_sequence && !current_header.empty()) {
        if (validateSequence(current_sequence)) {
            sequences.emplace_back(current_header, current_sequence);
        } else {
            std::cerr << "Advertencia: Secuencia inválida ignorada: " 
                      << current_header << std::endl;
        }
    }
    
    file.close();
    
    if (sequences.empty()) {
        std::cerr << "Error: No se encontraron secuencias válidas en " << filename << std::endl;
    } else {
        std::cout << "✓ Leídas " << sequences.size() << " secuencias de " << filename << std::endl;
    }
    
    return sequences;
}

void FastaIO::writeFasta(const std::vector<Sequence>& sequences, 
                        const std::string& filename, 
                        bool aligned) {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: No se pudo crear el archivo " << filename << std::endl;
        return;
    }
    
    const int line_width = aligned ? 80 : 80; // Ancho de línea para el formato FASTA
    
    for (const auto& seq : sequences) {
        file << ">" << seq.header << std::endl;
        
        // Escribir la secuencia en líneas de ancho fijo
        for (size_t i = 0; i < seq.sequence.length(); i += line_width) {
            size_t end = std::min(i + line_width, seq.sequence.length());
            file << seq.sequence.substr(i, end - i) << std::endl;
        }
    }
    
    file.close();
    std::cout << "✓ Guardadas " << sequences.size() << " secuencias en " << filename << std::endl;
}

bool FastaIO::validateSequence(const std::string& sequence) {
    if (sequence.empty()) {
        return false;
    }
    
    // Contar caracteres válidos vs inválidos
    size_t valid_chars = 0;
    const std::string valid_nucleotides = "ATCGRYSWKMBDHVN-";
    const std::string valid_amino_acids = "ABCDEFGHIKLMNPQRSTVWXYZ*-";
    
    for (char c : sequence) {
        char upper_c = std::toupper(c);
        if (valid_nucleotides.find(upper_c) != std::string::npos ||
            valid_amino_acids.find(upper_c) != std::string::npos) {
            valid_chars++;
        }
    }
    
    // Al menos 80% de caracteres válidos
    return (static_cast<double>(valid_chars) / sequence.length()) >= 0.8;
}

void FastaIO::printSequenceStats(const std::vector<Sequence>& sequences, 
                                const std::string& title) {
    if (sequences.empty()) {
        std::cout << "No hay secuencias para mostrar estadísticas." << std::endl;
        return;
    }
    
    std::cout << "\n=== " << title << " ===" << std::endl;
    std::cout << "Número de secuencias: " << sequences.size() << std::endl;
    
    // Calcular estadísticas de longitud
    size_t min_length = sequences[0].sequence.length();
    size_t max_length = sequences[0].sequence.length();
    size_t total_length = 0;
    
    for (const auto& seq : sequences) {
        size_t len = seq.sequence.length();
        min_length = std::min(min_length, len);
        max_length = std::max(max_length, len);
        total_length += len;
    }
    
    double avg_length = static_cast<double>(total_length) / sequences.size();
    
    std::cout << "Longitud mínima: " << min_length << std::endl;
    std::cout << "Longitud máxima: " << max_length << std::endl;
    std::cout << "Longitud promedio: " << std::fixed << std::setprecision(1) << avg_length << std::endl;
    
    // Mostrar algunas secuencias de ejemplo
    std::cout << "\nEjemplos de secuencias:" << std::endl;
    for (size_t i = 0; i < std::min(size_t(3), sequences.size()); ++i) {
        std::string preview = sequences[i].sequence.substr(0, 50);
        if (sequences[i].sequence.length() > 50) {
            preview += "...";
        }
        std::cout << "  " << sequences[i].header << ": " << preview << std::endl;
    }
    
    std::cout << std::endl;
}

std::string FastaIO::cleanLine(const std::string& line) {
    std::string cleaned;
    cleaned.reserve(line.length());
    
    for (char c : line) {
        if (c != '\r' && c != '\n' && c != '\t') {
            cleaned += c;
        }
    }
    
    // Quitar espacios al inicio y final
    size_t start = cleaned.find_first_not_of(' ');
    if (start == std::string::npos) {
        return "";
    }
    
    size_t end = cleaned.find_last_not_of(' ');
    return cleaned.substr(start, end - start + 1);
}