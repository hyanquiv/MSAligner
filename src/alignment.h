#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "io.h"
#include <vector>
#include <string>
#include <map>
#include <memory>

/**
 * Enumeración para los pasos del alineamiento
 */
enum class AlignmentStep {
    MATCH,
    DELETE,
    INSERT
};

/**
 * Estructura para representar un nodo en el �rbol gu�a
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
    std::vector<std::vector<double>> frequencies; // Frecuencias de cada base/amino�cido por posici�n
    std::vector<double> gap_frequencies;          // Frecuencias de gaps por posici�n
    int length;                                   // Longitud del perfil
    int num_sequences;                            // N�mero de secuencias en el perfil
    
    Profile() : length(0), num_sequences(0) {}
};

/**
 * Clase principal para el alineamiento m�ltiple de secuencias
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
     * Realiza el alineamiento m�ltiple de secuencias
     * @param sequences Vector de secuencias no alineadas
     * @return Vector de secuencias alineadas
     */
    std::vector<Sequence> alignSequences(const std::vector<Sequence>& sequences);
    
    /**
     * Obtiene estad�sticas del �ltimo alineamiento
     * @return Mapa con estad�sticas (gaps, longitud final, etc.)
     */
    std::map<std::string, int> getAlignmentStats() const;
    
    /**
     * Imprime el �rbol gu�a en consola
     */
    void printGuideTree() const;

private:
    // Matrices de puntuaci�n y par�metros
    int match_score;
    int mismatch_score;
    int gap_penalty;
    int gap_extension_penalty;
    
    // Estad�sticas del alineamiento
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
     * @return Distancia entre las secuencias (0.0 = id�nticas, 1.0 = completamente diferentes)
     */
    double calculateSequenceDistance(const std::string& seq1, const std::string& seq2);
    
    /**
     * Construye el �rbol gu�a usando UPGMA
     * @param sequences Vector de secuencias originales
     * @param distance_matrix Matriz de distancias
     * @return Nodo ra�z del �rbol gu�a
     */
    std::shared_ptr<TreeNode> buildGuideTree(const std::vector<Sequence>& sequences,
                                           const std::vector<std::vector<double>>& distance_matrix);
    
    /**
     * Realiza el alineamiento progresivo siguiendo el �rbol gu�a
     * @param sequences Secuencias originales
     * @param node Nodo actual del �rbol
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
     * Imprime un nodo del �rbol recursivamente
     * @param node Nodo a imprimir
     * @param depth Profundidad actual
     */
    void printTreeNode(const std::shared_ptr<TreeNode>& node, int depth = 0) const;
    
    /**
     * Encuentra el �ndice del car�cter en el alfabeto
     * @param c Car�cter a buscar
     * @return �ndice en el alfabeto (0-3 para DNA, 0-19 para prote�nas)
     */
    int getAlphabetIndex(char c) const;
    
    /**
     * Obtiene el car�cter del alfabeto por �ndice
     * @param index �ndice en el alfabeto
     * @return Car�cter correspondiente
     */
    char getAlphabetChar(int index) const;
    
    /**
     * Inicializa la matriz de programación dinámica
     */
    std::vector<std::vector<int>> initializeDPMatrix(size_t m, size_t n);
    
    /**
     * Llena la matriz de programación dinámica
     */
    void fillDPMatrix(std::vector<std::vector<int>>& dp, 
                     const std::string& seq1, const std::string& seq2,
                     size_t m, size_t n);
    
    /**
     * Calcula el puntaje de coincidencia entre dos caracteres
     */
    int calculateMatchScore(char c1, char c2);
    
    /**
     * Reconstruye el alineamiento a partir de la matriz DP
     */
    std::pair<std::string, std::string> reconstructAlignment(
        const std::vector<std::vector<int>>& dp,
        const std::string& seq1, const std::string& seq2,
        size_t m, size_t n);
    
    /**
     * Determina el próximo paso en la reconstrucción del alineamiento
     */
    AlignmentStep determineAlignmentStep(
        const std::vector<std::vector<int>>& dp,
        const std::string& seq1, const std::string& seq2,
        size_t i, size_t j);
    
    /**
     * Verifica si el paso actual es una coincidencia/desajuste
     */
    bool isMatchStep(const std::vector<std::vector<int>>& dp,
                    const std::string& seq1, const std::string& seq2,
                    size_t i, size_t j);
    
    /**
     * Verifica si el paso actual es una eliminación
     */
    bool isDeleteStep(const std::vector<std::vector<int>>& dp,
                     size_t i, size_t j);
    
    std::string generateConsensusFromProfile(const Profile& profile);
    char findBestCharacterAtPosition(const Profile& profile, int pos);
    Profile initializeCombinedProfile(const std::pair<std::string, std::string>& aligned_pair, 
                                    const Profile& profile);
    void fillCombinedProfile(Profile& new_profile, 
                           const std::pair<std::string, std::string>& aligned_pair,
                           const Profile& original_profile,
                           const std::string& sequence);
    void addOriginalProfileFrequencies(Profile& new_profile, const Profile& original_profile,
                                      char profile_char, int pos, const std::string& aligned_consensus);
    int findOriginalPosition(const std::string& aligned_consensus, int aligned_pos);
    void copyFrequenciesFromOriginal(Profile& new_profile, const Profile& original_profile,
                                   int new_pos, int orig_pos);
    void addNewSequenceFrequencies(Profile& new_profile, char seq_char, int pos);
    void normalizeFrequenciesAtPosition(Profile& new_profile, int pos);
    
    // Constantes
    static const std::string DNA_ALPHABET;
    static const std::string PROTEIN_ALPHABET;
    static const int ALPHABET_SIZE;
};

#endif // ALIGNMENT_H