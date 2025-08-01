#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "io.h"
#include <vector>
#include <string>
#include <map>
#include <memory>

/**
 * Estructura para representar un nodo en el árbol guía
 */
struct TreeNode {
    int id;                                    // ID del nodo (-1 para nodos internos)
    double distance;                           // Distancia acumulada
    std::shared_ptr<TreeNode> left;           // Hijo izquierdo
    std::shared_ptr<TreeNode> right;          // Hijo derecho
    std::vector<int> sequences;               // IDs de secuencias en este nodo
    
    TreeNode(int id = -1) : id(id), distance(0.0), left(nullptr), right(nullptr) {}
};

/**
 * Estructura para representar un perfil de alineamiento
 */
struct Profile {
    std::vector<std::vector<double>> frequencies; // Frecuencias de cada base/aminoácido por posición
    std::vector<double> gap_frequencies;          // Frecuencias de gaps por posición
    int length;                                   // Longitud del perfil
    int num_sequences;                            // Número de secuencias en el perfil
    
    Profile() : length(0), num_sequences(0) {}
};

/**
 * Clase principal para el alineamiento múltiple de secuencias
 */
class MSAAligner {
public:
    /**
     * Constructor
     */
    MSAAligner();
    
    /**
     * Destructor
     */
    ~MSAAligner() = default;
    
    /**
     * Realiza el alineamiento múltiple de secuencias
     * @param sequences Vector de secuencias no alineadas
     * @return Vector de secuencias alineadas
     */
    std::vector<Sequence> alignSequences(const std::vector<Sequence>& sequences);
    
    /**
     * Obtiene estadísticas del último alineamiento
     * @return Mapa con estadísticas (gaps, longitud final, etc.)
     */
    std::map<std::string, int> getAlignmentStats() const;
    
    /**
     * Imprime el árbol guía en consola
     */
    void printGuideTree() const;

private:
    // Matrices de puntuación y parámetros
    int match_score;
    int mismatch_score;
    int gap_penalty;
    int gap_extension_penalty;
    
    // Estadísticas del alineamiento
    int total_gaps;
    int final_length;
    std::shared_ptr<TreeNode> guide_tree;
    
    /**
     * Calcula la matriz de distancias entre todas las secuencias
     * @param sequences Vector de secuencias
     * @return Matriz de distancias
     */
    std::vector<std::vector<double>> calculateDistanceMatrix(const std::vector<Sequence>& sequences);
    
    /**
     * Calcula la distancia entre dos secuencias usando identidad porcentual
     * @param seq1 Primera secuencia
     * @param seq2 Segunda secuencia
     * @return Distancia entre las secuencias (0.0 = idénticas, 1.0 = completamente diferentes)
     */
    double calculateSequenceDistance(const std::string& seq1, const std::string& seq2);
    
    /**
     * Construye el árbol guía usando UPGMA
     * @param sequences Vector de secuencias originales
     * @param distance_matrix Matriz de distancias
     * @return Nodo raíz del árbol guía
     */
    std::shared_ptr<TreeNode> buildGuideTree(const std::vector<Sequence>& sequences,
                                           const std::vector<std::vector<double>>& distance_matrix);
    
    /**
     * Realiza el alineamiento progresivo siguiendo el árbol guía
     * @param sequences Secuencias originales
     * @param node Nodo actual del árbol
     * @return Perfil del alineamiento en este nodo
     */
    Profile progressiveAlignment(const std::vector<Sequence>& sequences,
                               const std::shared_ptr<TreeNode>& node);
    
    /**
     * Alinea dos secuencias usando Needleman-Wunsch
     * @param seq1 Primera secuencia
     * @param seq2 Segunda secuencia
     * @return Par de secuencias alineadas
     */
    std::pair<std::string, std::string> pairwiseAlignment(const std::string& seq1,
                                                         const std::string& seq2);
    
    /**
     * Alinea una secuencia con un perfil
     * @param sequence Secuencia a alinear
     * @param profile Perfil existente
     * @return Nuevo perfil que incluye la secuencia alineada
     */
    Profile alignSequenceToProfile(const std::string& sequence, const Profile& profile);
    
    /**
     * Alinea dos perfiles
     * @param profile1 Primer perfil
     * @param profile2 Segundo perfil
     * @return Perfil combinado
     */
    Profile alignProfiles(const Profile& profile1, const Profile& profile2);
    
    /**
     * Convierte un perfil final a secuencias alineadas
     * @param profile Perfil final del alineamiento
     * @param sequences Secuencias originales
     * @param sequence_order Orden de las secuencias en el perfil
     * @return Vector de secuencias alineadas
     */
    std::vector<Sequence> profileToSequences(const Profile& profile,
                                           const std::vector<Sequence>& sequences,
                                           const std::vector<int>& sequence_order);
    
    /**
     * Crea un perfil a partir de una sola secuencia
     * @param sequence Secuencia base
     * @return Perfil creado
     */
    Profile createProfile(const std::string& sequence);
    
    /**
     * Imprime un nodo del árbol recursivamente
     * @param node Nodo a imprimir
     * @param depth Profundidad actual
     */
    void printTreeNode(const std::shared_ptr<TreeNode>& node, int depth = 0) const;
    
    /**
     * Encuentra el índice del carácter en el alfabeto
     * @param c Carácter a buscar
     * @return Índice en el alfabeto (0-3 para DNA, 0-19 para proteínas)
     */
    int getAlphabetIndex(char c) const;
    
    /**
     * Obtiene el carácter del alfabeto por índice
     * @param index Índice en el alfabeto
     * @return Carácter correspondiente
     */
    char getAlphabetChar(int index) const;
    
    // Constantes
    static const std::string DNA_ALPHABET;
    static const std::string PROTEIN_ALPHABET;
    static const int ALPHABET_SIZE;
};

#endif // ALIGNMENT_H