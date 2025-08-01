# 🧬 Alineador Múltiple de Secuencias (MSA)

Implementación en C++ de un alineador múltiple de secuencias que utiliza algoritmo progresivo con árbol guía UPGMA.

## 🎯 Características

- ✅ Lectura de archivos FASTA estándar
- ✅ Algoritmo UPGMA para construcción del árbol guía
- ✅ Alineamiento progresivo con programación dinámica
- ✅ Soporte para secuencias de ADN y proteínas
- ✅ Estadísticas detalladas del alineamiento
- ✅ Interfaz de línea de comandos simple
- ✅ Código modular y bien documentado

## 🏗️ Estructura del Proyecto

```
MSAAligner/
├── CMakeLists.txt          # Configuración de compilación
├── README.md               # Este archivo
├── src/
│   ├── main.cpp           # Programa principal
│   ├── alignment.h        # Definiciones del alineador
│   ├── alignment.cpp      # Implementación del alineador
│   ├── io.h              # Definiciones de E/S FASTA
│   └── io.cpp            # Implementación de E/S FASTA
└── build/                 # Directorio de compilación (generado)
```

## 🛠️ Compilación

### Requisitos

- CMake 3.10 o superior
- Compilador C++17 compatible (GCC, Clang, MSVC)
- Sistema operativo: Linux, macOS, Windows

### Pasos de compilación

```bash
# Clonar o descargar el proyecto
git clone <repository_url>
cd MSAAligner

# Crear directorio de compilación
mkdir build
cd build

# Configurar con CMake
cmake ..

# Compilar
make

# En Windows con Visual Studio:
# cmake --build . --config Release
```

### Compilación alternativa (sin CMake)

```bash
# Compilación directa con g++
g++ -std=c++17 -O3 -Wall -Wextra \
    src/main.cpp src/alignment.cpp src/io.cpp \
    -o alineador
```

## 🚀 Uso

### Sintaxis básica

```bash
./alineador <archivo_entrada.fasta> <archivo_salida.fasta>
```

### Ejemplo

```bash
# Alinear secuencias del archivo input.fasta
./alineador sequences.fasta aligned_sequences.fasta
```

### Formato de entrada

El archivo de entrada debe ser un FASTA estándar con múltiples secuencias:

```fasta
>Secuencia_1
ATCGATCGATCG
>Secuencia_2
ATCGAACGATCG
>Secuencia_3
ATCGATCGAACG
```

### Formato de salida

El archivo de salida contendrá las secuencias alineadas con gaps (`-`):

```fasta
>Secuencia_1
ATCGATCGATCG
>Secuencia_2
ATCGAACGATCG
>Secuencia_3
ATCGATCGAACG
```

## 🔬 Algoritmo

El alineador implementa un algoritmo MSA progresivo en tres pasos:

### 1. Matriz de Distancias
- Calcula distancias par a par entre todas las secuencias
- Utiliza identidad porcentual como métrica
- Complejidad: O(n² × L²) donde n = número de secuencias, L = longitud promedio

### 2. Árbol Guía (UPGMA)
- Construye un árbol filogenético usando UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
- Determina el orden óptimo de alineamiento
- Complejidad: O(n³)

### 3. Alineamiento Progresivo
- Alinea secuencias siguiendo el árbol guía
- Utiliza programación dinámica (Needleman-Wunsch)
- Combina perfiles de secuencias gradualmente
- Complejidad: O(n × L²)

## 📊 Parámetros de Alineamiento

| Parámetro | Valor | Descripción |
|-----------|-------|-------------|
| Match Score | +2 | Puntuación por coincidencia |
| Mismatch Score | -1 | Penalización por no coincidencia |
| Gap Penalty | -2 | Penalización por abrir gap |
| Gap Extension | -1 | Penalización por extender gap |

## 📈 Estadísticas de Salida

El programa proporciona información detallada:

- ✅ Tiempo de ejecución total
- ✅ Número de secuencias procesadas
- ✅ Longitud final del alineamiento
- ✅ Número total de gaps insertados
- ✅ Porcentaje de gaps en el alineamiento
- ✅ Visualización del árbol guía

## 🔧 Personalización

### Modificar parámetros de puntuación

Edita los valores en el constructor de `MSAAligner` (archivo `alignment.cpp`):

```cpp
MSAAligner::MSAAligner() 
    : match_score(2),           // Cambiar aquí
      mismatch_score(-1),       // Cambiar aquí
      gap_penalty(-2),          // Cambiar aquí
      gap_extension_penalty(-1) // Cambiar aquí
```

### Cambiar alfabeto de secuencias

Modifica las constantes en `alignment.h`:

```cpp
const std::string DNA_ALPHABET = "ATCG";           // Para ADN
const std::string PROTEIN_ALPHABET = "ARNDCQEGH..."; // Para proteínas
```

## 🧪 Casos de Prueba

### Secuencias de ejemplo

Crea un archivo `test.fasta`:

```fasta
>seq1
ATCGATCGATCG
>seq2
ATCGAACGATCG
>seq3
ATCGATCGAACG
>seq4
ATCGAACGAACG
```

Ejecuta:

```bash
./alineador test.fasta result.fasta
```

## ⚠️ Limitaciones

- **Complejidad computacional**: O(n³ + n×L²) - no óptimo para datasets muy grandes
- **Memoria**: Almacena matrices completas en memoria
- **Algoritmo simple**: Implementación básica, no incluye optimizaciones avanzadas
- **Gap penalties**: Modelo lineal simple, no considera estructura secundaria

## 🔮 Mejoras Futuras

- [ ] Algoritmo de alineamiento más sofisticado (T-Coffee, MUSCLE)
- [ ] Soporte para perfiles HMM
- [ ] Paralelización con OpenMP
- [ ] Interfaz gráfica
- [ ] Formato de salida en otros formatos (Clustal, Phylip)
- [ ] Análisis de calidad del alineamiento
- [ ] Soporte para secuencias muy largas (streaming)

## 📝 Licencia

Este proyecto está bajo licencia MIT. Ver archivo LICENSE para más detalles.

## 🤝 Contribuciones

Las contribuciones son bienvenidas. Por favor:

1. Fork el proyecto
2. Crea una rama para tu feature (`git checkout -b feature/AmazingFeature`)
3. Commit tus cambios (`git commit -m 'Add some AmazingFeature'`)
4. Push a la rama (`git push origin feature/AmazingFeature`)
5. Abre un Pull Request

## 📞 Contacto

Para preguntas, sugerencias o reportar bugs, abre un issue en el repositorio.

---

**¡Gracias por usar el Alineador Múltiple de Secuencias!** 🧬✨