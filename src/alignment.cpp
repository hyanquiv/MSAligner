﻿#include "alignment.h"
#include <algorithm>
#include <climits>
#include <iostream>
#include <iomanip>
#include <cmath>

// Definición de constantes estáticas
const std::string MSAAligner::DNA_ALPHABET = "ATCG";
const std::string MSAAligner::PROTEIN_ALPHABET = "ARNDCQEGHILKMFPSTWYV";
const int MSAAligner::ALPHABET_SIZE = 4; // Usaremos DNA por simplicidad

MSAAligner::MSAAligner() 
    : match_score(2), mismatch_score(-1), gap_penalty(-2), gap_extension_penalty(-1),
      total_gaps(0), final_length(0), guide_tree(nullptr) {
}

std::vector<Sequence> MSAAligner::alignSequences(const std::vector<Sequence>& sequences) {
    if (sequences.size() < 2) {
        std::cerr << "Error: Se necesitan al menos 2 secuencias para el alineamiento." << std::endl;
        return sequences;
    }

    std::cout << "\nIniciando alineamiento multiple de secuencias..." << std::endl;
    std::cout << "Numero de secuencias: " << sequences.size() << std::endl;

    // Reiniciar estadisticas
    total_gaps = 0;
    final_length = 0;

    // Paso 1: Calcular matriz de distancias
    std::cout << "Calculando matriz de distancias..." << std::endl;
    auto distance_matrix = calculateDistanceMatrix(sequences);

    // Paso 2: Construir arbol guia
    std::cout << "Construyendo arbol guia con UPGMA..." << std::endl;
    guide_tree = buildGuideTree(sequences, distance_matrix);

    // Paso 3: Alineamiento progresivo
    std::cout << "Realizando alineamiento progresivo..." << std::endl;
    Profile final_profile = progressiveAlignment(sequences, guide_tree);

    // Paso 4: Convertir perfil a secuencias
    std::cout << "Generando secuencias alineadas..." << std::endl;
    std::vector<int> sequence_order;
    for (int i = 0; i < static_cast<int>(sequences.size()); ++i) {
        sequence_order.push_back(i);
    }

    auto aligned_sequences = profileToSequences(final_profile, sequences, sequence_order);

    // Actualizar estadisticas
    if (!aligned_sequences.empty()) {
        final_length = aligned_sequences[0].sequence.length();

        // Contar gaps totales
        for (const auto& seq : aligned_sequences) {
            total_gaps += std::count(seq.sequence.begin(), seq.sequence.end(), '-');
        }
    }

    std::cout << "Alineamiento completado!" << std::endl;
    std::cout << "Longitud final: " << final_length << std::endl;
    std::cout << "Gaps totales insertados: " << total_gaps << std::endl;

    return aligned_sequences;
}

std::vector<std::vector<double>> MSAAligner::calculateDistanceMatrix(const std::vector<Sequence>& sequences) {
    size_t n = sequences.size();
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 0.0));
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double distance = calculateSequenceDistance(sequences[i].sequence, sequences[j].sequence);
            matrix[i][j] = matrix[j][i] = distance;
        }
    }
    
    return matrix;
}

double MSAAligner::calculateSequenceDistance(const std::string& seq1, const std::string& seq2) {
    if (seq1.empty() || seq2.empty()) {
        return 1.0; // Máxima distancia
    }
    
    size_t min_length = std::min(seq1.length(), seq2.length());
    size_t max_length = std::max(seq1.length(), seq2.length());
    
    if (min_length == 0) {
        return 1.0;
    }
    
    // Contar coincidencias en las posiciones superpuestas
    size_t matches = 0;
    for (size_t i = 0; i < min_length; ++i) {
        if (std::toupper(seq1[i]) == std::toupper(seq2[i])) {
            matches++;
        }
    }
    
    // Calcular identidad considerando la diferencia de longitud
    double identity = static_cast<double>(matches) / max_length;
    return 1.0 - identity; // Convertir identidad a distancia
}

std::shared_ptr<TreeNode> MSAAligner::buildGuideTree(const std::vector<Sequence>& sequences,
                                                     const std::vector<std::vector<double>>& distance_matrix) {
    size_t n = sequences.size();
    
    // Crear nodos hoja
    std::vector<std::shared_ptr<TreeNode>> nodes;
    for (size_t i = 0; i < n; ++i) {
        auto node = std::make_shared<TreeNode>(static_cast<int>(i));
        node->sequences.push_back(static_cast<int>(i));
        nodes.push_back(node);
    }

    // Algoritmo UPGMA
    while (nodes.size() > 1) {
        size_t min_i = 0, min_j = 1;
        double min_distance = std::numeric_limits<double>::max();

        // Buscar los nodos más cercanos
        for (size_t i = 0; i < nodes.size(); ++i) {
            for (size_t j = i + 1; j < nodes.size(); ++j) {
                double dist = 0.0;
                for (int si : nodes[i]->sequences) {
                    for (int sj : nodes[j]->sequences) {
                        dist += distance_matrix[si][sj];
                    }
                }
                dist /= (nodes[i]->sequences.size() * nodes[j]->sequences.size());

                if (dist < min_distance) {
                    min_distance = dist;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        // Crear nuevo nodo interno
        auto new_node = std::make_shared<TreeNode>();
        new_node->distance = min_distance / 2.0;
        new_node->left = nodes[min_i];
        new_node->right = nodes[min_j];

        // Combinar las secuencias
        new_node->sequences = nodes[min_i]->sequences;
        new_node->sequences.insert(new_node->sequences.end(),
                                   nodes[min_j]->sequences.begin(),
                                   nodes[min_j]->sequences.end());

        // Reemplazar el nodo min_i con el nuevo nodo y eliminar min_j
        if (min_i < min_j) {
            nodes.erase(nodes.begin() + min_j);
            nodes.erase(nodes.begin() + min_i);
        } else {
            nodes.erase(nodes.begin() + min_i);
            nodes.erase(nodes.begin() + min_j);
        }

        nodes.push_back(new_node);
    }

    return nodes.empty() ? nullptr : nodes.front();
}

Profile MSAAligner::progressiveAlignment(const std::vector<Sequence>& sequences,
                                       const std::shared_ptr<TreeNode>& node) {
    if (!node) {
        return Profile();
    }
    
    // Nodo hoja - crear perfil de una sola secuencia
    if (!node->sequences.empty() && !node->left && !node->right) {
        int seq_idx = node->sequences[0];
        return createProfile(sequences[seq_idx].sequence);
    }
    
    // Nodo interno - alinear subperfiles
    if (node->left && node->right) {
        Profile left_profile = progressiveAlignment(sequences, node->left);
        Profile right_profile = progressiveAlignment(sequences, node->right);
        return alignProfiles(left_profile, right_profile);
    }
    
    return Profile();
}

std::pair<std::string, std::string> MSAAligner::pairwiseAlignment(const std::string& seq1,
                                                                const std::string& seq2) {
    size_t m = seq1.length();
    size_t n = seq2.length();
    
    std::vector<std::vector<int>> dp = initializeDPMatrix(m, n);
    fillDPMatrix(dp, seq1, seq2, m, n);
    
    return reconstructAlignment(dp, seq1, seq2, m, n);
}

std::vector<std::vector<int>> MSAAligner::initializeDPMatrix(size_t m, size_t n) {
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    
    for (size_t i = 0; i <= m; ++i) {
        dp[i][0] = static_cast<int>(i) * gap_penalty;
    }
    for (size_t j = 0; j <= n; ++j) {
        dp[0][j] = static_cast<int>(j) * gap_penalty;
    }
    
    return dp;
}

void MSAAligner::fillDPMatrix(std::vector<std::vector<int>>& dp, 
                             const std::string& seq1, const std::string& seq2,
                             size_t m, size_t n) {
    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            int match_score_val = calculateMatchScore(seq1[i-1], seq2[j-1]);
            int match = dp[i-1][j-1] + match_score_val;
            int delete_op = dp[i-1][j] + gap_penalty;
            int insert_op = dp[i][j-1] + gap_penalty;
            
            dp[i][j] = std::max({match, delete_op, insert_op});
        }
    }
}

int MSAAligner::calculateMatchScore(char c1, char c2) {
    return (std::toupper(c1) == std::toupper(c2)) ? match_score : mismatch_score;
}

std::pair<std::string, std::string> MSAAligner::reconstructAlignment(
    const std::vector<std::vector<int>>& dp,
    const std::string& seq1, const std::string& seq2,
    size_t m, size_t n) {
    
    std::string aligned_seq1, aligned_seq2;
    size_t i = m, j = n;
    
    while (i > 0 || j > 0) {
        AlignmentStep step = determineAlignmentStep(dp, seq1, seq2, i, j);
        
        switch (step) {
            case AlignmentStep::MATCH:
                aligned_seq1 = seq1[i-1] + aligned_seq1;
                aligned_seq2 = seq2[j-1] + aligned_seq2;
                i--; j--;
                break;
            case AlignmentStep::DELETE:
                aligned_seq1 = seq1[i-1] + aligned_seq1;
                aligned_seq2 = '-' + aligned_seq2;
                i--;
                break;
            case AlignmentStep::INSERT:
                aligned_seq1 = '-' + aligned_seq1;
                aligned_seq2 = seq2[j-1] + aligned_seq2;
                j--;
                break;
        }
    }
    
    return {aligned_seq1, aligned_seq2};
}

AlignmentStep MSAAligner::determineAlignmentStep(
    const std::vector<std::vector<int>>& dp,
    const std::string& seq1, const std::string& seq2,
    size_t i, size_t j) {
    
    if (i > 0 && j > 0 && isMatchStep(dp, seq1, seq2, i, j)) {
        return AlignmentStep::MATCH;
    }
    if (i > 0 && isDeleteStep(dp, i, j)) {
        return AlignmentStep::DELETE;
    }
    return AlignmentStep::INSERT;
}

bool MSAAligner::isMatchStep(const std::vector<std::vector<int>>& dp,
                            const std::string& seq1, const std::string& seq2,
                            size_t i, size_t j) {
    int match_score_val = calculateMatchScore(seq1[i-1], seq2[j-1]);
    return dp[i][j] == dp[i-1][j-1] + match_score_val;
}

bool MSAAligner::isDeleteStep(const std::vector<std::vector<int>>& dp,
                             size_t i, size_t j) {
    return dp[i][j] == dp[i-1][j] + gap_penalty;
}

Profile MSAAligner::alignSequenceToProfile(const std::string& sequence, const Profile& profile) {
    std::string consensus = generateConsensusFromProfile(profile);
    auto aligned_pair = pairwiseAlignment(consensus, sequence);
    
    Profile new_profile = initializeCombinedProfile(aligned_pair, profile);
    fillCombinedProfile(new_profile, aligned_pair, profile, sequence);
    
    return new_profile;
}

std::string MSAAligner::generateConsensusFromProfile(const Profile& profile) {
    std::string consensus;
    for (int pos = 0; pos < profile.length; ++pos) {
        consensus += findBestCharacterAtPosition(profile, pos);
    }
    return consensus;
}

char MSAAligner::findBestCharacterAtPosition(const Profile& profile, int pos) {
    char best_char = 'A';
    double best_freq = 0.0;
    
    for (int base = 0; base < ALPHABET_SIZE; ++base) {
        if (profile.frequencies[pos][base] > best_freq) {
            best_freq = profile.frequencies[pos][base];
            best_char = getAlphabetChar(base);
        }
    }
    return best_char;
}

Profile MSAAligner::initializeCombinedProfile(const std::pair<std::string, std::string>& aligned_pair, 
                                            const Profile& profile) {
    Profile new_profile;
    new_profile.length = aligned_pair.first.length();
    new_profile.num_sequences = profile.num_sequences + 1;
    new_profile.frequencies.resize(new_profile.length, std::vector<double>(ALPHABET_SIZE, 0.0));
    new_profile.gap_frequencies.resize(new_profile.length, 0.0);
    return new_profile;
}

void MSAAligner::fillCombinedProfile(Profile& new_profile, 
                                   const std::pair<std::string, std::string>& aligned_pair,
                                   const Profile& original_profile,
                                   const std::string& sequence) {
    for (int pos = 0; pos < new_profile.length; ++pos) {
        char profile_char = aligned_pair.first[pos];
        char seq_char = aligned_pair.second[pos];
        
        addOriginalProfileFrequencies(new_profile, original_profile, profile_char, pos, aligned_pair.first);
        addNewSequenceFrequencies(new_profile, seq_char, pos);
        normalizeFrequenciesAtPosition(new_profile, pos);
    }
}

void MSAAligner::addOriginalProfileFrequencies(Profile& new_profile, const Profile& original_profile,
                                              char profile_char, int pos, const std::string& aligned_consensus) {
    if (profile_char == '-') return;
    
    int orig_pos = findOriginalPosition(aligned_consensus, pos);
    if (orig_pos >= 0 && orig_pos < original_profile.length) {
        copyFrequenciesFromOriginal(new_profile, original_profile, pos, orig_pos);
    }
}

int MSAAligner::findOriginalPosition(const std::string& aligned_consensus, int aligned_pos) {
    int orig_pos = 0;
    for (int i = 0; i <= aligned_pos; ++i) {
        if (aligned_consensus[i] != '-') {
            if (i == aligned_pos) return orig_pos;
            orig_pos++;
        }
    }
    return -1;
}

void MSAAligner::copyFrequenciesFromOriginal(Profile& new_profile, const Profile& original_profile,
                                           int new_pos, int orig_pos) {
    for (int base = 0; base < ALPHABET_SIZE; ++base) {
        new_profile.frequencies[new_pos][base] += 
            original_profile.frequencies[orig_pos][base] * original_profile.num_sequences;
    }
    new_profile.gap_frequencies[new_pos] += 
        original_profile.gap_frequencies[orig_pos] * original_profile.num_sequences;
}

void MSAAligner::addNewSequenceFrequencies(Profile& new_profile, char seq_char, int pos) {
    if (seq_char == '-') {
        new_profile.gap_frequencies[pos] += 1.0;
    } else {
        int base_idx = getAlphabetIndex(seq_char);
        if (base_idx >= 0) {
            new_profile.frequencies[pos][base_idx] += 1.0;
        }
    }
}

void MSAAligner::normalizeFrequenciesAtPosition(Profile& new_profile, int pos) {
    for (int base = 0; base < ALPHABET_SIZE; ++base) {
        new_profile.frequencies[pos][base] /= new_profile.num_sequences;
    }
    new_profile.gap_frequencies[pos] /= new_profile.num_sequences;
}

Profile MSAAligner::alignProfiles(const Profile& profile1, const Profile& profile2) {
    // Simplificación: convertir perfiles a secuencias consenso y alinear
    std::string consensus1, consensus2;
    
    for (int pos = 0; pos < profile1.length; ++pos) {
        char best_char = 'A';
        double best_freq = 0.0;
        
        for (int base = 0; base < ALPHABET_SIZE; ++base) {
            if (profile1.frequencies[pos][base] > best_freq) {
                best_freq = profile1.frequencies[pos][base];
                best_char = getAlphabetChar(base);
            }
        }
        consensus1 += best_char;
    }
    
    for (int pos = 0; pos < profile2.length; ++pos) {
        char best_char = 'A';
        double best_freq = 0.0;
        
        for (int base = 0; base < ALPHABET_SIZE; ++base) {
            if (profile2.frequencies[pos][base] > best_freq) {
                best_freq = profile2.frequencies[pos][base];
                best_char = getAlphabetChar(base);
            }
        }
        consensus2 += best_char;
    }
    
    auto aligned_pair = pairwiseAlignment(consensus1, consensus2);
    
    // Crear perfil combinado
    Profile combined_profile;
    combined_profile.length = aligned_pair.first.length();
    combined_profile.num_sequences = profile1.num_sequences + profile2.num_sequences;
    combined_profile.frequencies.resize(combined_profile.length, std::vector<double>(ALPHABET_SIZE, 0.0));
    combined_profile.gap_frequencies.resize(combined_profile.length, 0.0);
    
    // Combinar los perfiles basándose en el alineamiento
    int pos1 = 0, pos2 = 0;
    for (int pos = 0; pos < combined_profile.length; ++pos) {
        char char1 = aligned_pair.first[pos];
        char char2 = aligned_pair.second[pos];
        
        // Agregar contribución del primer perfil
        if (char1 != '-' && pos1 < profile1.length) {
            for (int base = 0; base < ALPHABET_SIZE; ++base) {
                combined_profile.frequencies[pos][base] += 
                    profile1.frequencies[pos1][base] * profile1.num_sequences;
            }
            combined_profile.gap_frequencies[pos] += 
                profile1.gap_frequencies[pos1] * profile1.num_sequences;
            pos1++;
        }
        
        // Agregar contribución del segundo perfil
        if (char2 != '-' && pos2 < profile2.length) {
            for (int base = 0; base < ALPHABET_SIZE; ++base) {
                combined_profile.frequencies[pos][base] += 
                    profile2.frequencies[pos2][base] * profile2.num_sequences;
            }
            combined_profile.gap_frequencies[pos] += 
                profile2.gap_frequencies[pos2] * profile2.num_sequences;
            pos2++;
        }
        
        // Normalizar frecuencias
        for (int base = 0; base < ALPHABET_SIZE; ++base) {
            combined_profile.frequencies[pos][base] /= combined_profile.num_sequences;
        }
        combined_profile.gap_frequencies[pos] /= combined_profile.num_sequences;
    }
    
    return combined_profile;
}

std::vector<Sequence> MSAAligner::profileToSequences(const Profile& profile,
                                                   const std::vector<Sequence>& sequences,
                                                   const std::vector<int>& sequence_order) {
    (void)sequence_order; // Suprimir advertencia de parámetro no utilizado
    std::vector<Sequence> aligned_sequences;
    
    // Simplificación: generar secuencias alineadas basadas en el perfil
    // En una implementación completa, se mantendría el rastro de cada secuencia individual
    
    for (size_t i = 0; i < sequences.size(); ++i) {
        Sequence aligned_seq;
        aligned_seq.header = sequences[i].header;
        
        // Para esta implementación simplificada, alineamos cada secuencia individual al perfil
        std::string consensus;
        for (int pos = 0; pos < profile.length; ++pos) {
            char best_char = 'A';
            double best_freq = 0.0;
            
            for (int base = 0; base < ALPHABET_SIZE; ++base) {
                if (profile.frequencies[pos][base] > best_freq) {
                    best_freq = profile.frequencies[pos][base];
                    best_char = getAlphabetChar(base);
                }
            }
            consensus += best_char;
        }
        
        auto aligned_pair = pairwiseAlignment(sequences[i].sequence, consensus);
        aligned_seq.sequence = aligned_pair.first;
        
        aligned_sequences.push_back(aligned_seq);
    }
    
    return aligned_sequences;
}

Profile MSAAligner::createProfile(const std::string& sequence) {
    Profile profile;
    profile.length = sequence.length();
    profile.num_sequences = 1;
    profile.frequencies.resize(profile.length, std::vector<double>(ALPHABET_SIZE, 0.0));
    profile.gap_frequencies.resize(profile.length, 0.0);
    
    for (int pos = 0; pos < profile.length; ++pos) {
        char c = sequence[pos];
        if (c == '-') {
            profile.gap_frequencies[pos] = 1.0;
        } else {
            int base_idx = getAlphabetIndex(c);
            if (base_idx >= 0) {
                profile.frequencies[pos][base_idx] = 1.0;
            }
        }
    }
    
    return profile;
}

std::map<std::string, int> MSAAligner::getAlignmentStats() const {
    std::map<std::string, int> stats;
    stats["total_gaps"] = total_gaps;
    stats["final_length"] = final_length;
    return stats;
}

void MSAAligner::printGuideTree() const {
    if (!guide_tree) {
        std::cout << "No hay arbol guia disponible." << std::endl;
        return;
    }

    std::cout << "\nArbol Guia (UPGMA):" << std::endl;
    printTreeNode(guide_tree, 0);
    std::cout << std::endl;
}

void MSAAligner::printTreeNode(const std::shared_ptr<TreeNode>& node, int depth) const {
    if (!node) return;
    
    std::string indent(depth * 2, ' ');
    
    if (node->id >= 0) {
        // Nodo hoja
        std::cout << indent << "├─ Secuencia " << node->id << " (dist: " 
                  << std::fixed << std::setprecision(3) << node->distance << ")" << std::endl;
    } else {
        // Nodo interno
        std::cout << indent << "├─ Nodo interno (dist: " 
                  << std::fixed << std::setprecision(3) << node->distance << ")" << std::endl;
        if (node->left) {
            printTreeNode(node->left, depth + 1);
        }
        if (node->right) {
            printTreeNode(node->right, depth + 1);
        }
    }
}

int MSAAligner::getAlphabetIndex(char c) const {
    char upper_c = std::toupper(c);
    size_t pos = DNA_ALPHABET.find(upper_c);
    return (pos != std::string::npos) ? static_cast<int>(pos) : -1;
}

char MSAAligner::getAlphabetChar(int index) const {
    if (index >= 0 && index < static_cast<int>(DNA_ALPHABET.length())) {
        return DNA_ALPHABET[index];
    }
    return 'N'; // Carácter desconocido
}