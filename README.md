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

- **Visual Studio Community 2022** (Windows)  
- Compilador C++17 compatible  
- Sistema operativo compatible: **Windows (recomendado con MSVC)**, Linux, macOS  
- CMake (recomendado para portabilidad entre plataformas)

### Pasos de compilación

#### 🪟 En Windows con Visual Studio Community 2022 (MSVC)

1. Abre la terminal "x64 Native Tools Command Prompt for VS 2022"
2. Clona el repositorio y crea un directorio de compilación:

```bash
git clone <repository_url>
cd MSAAligner
mkdir build
cd build
```

3. Configura el proyecto con CMake:

```bash
cmake .. -G "Visual Studio 17 2022"
```

4. Compila el proyecto en modo Release:

```bash
cmake --build . --config Release
```

Esto generará el ejecutable `alineador.exe` en el subdirectorio `Release`.

#### 🐧 En Linux/macOS con g++

```bash
# Compilación directa sin CMake (requiere g++)
g++ -std=c++17 -O3 -Wall -Wextra     src/main.cpp src/alignment.cpp src/io.cpp     -o alineador
```

O bien con CMake:

```bash
mkdir build
cd build
cmake ..
make
```

## 🚀 Uso

### Sintaxis básica

```bash
./alineador <archivo_entrada.fasta> <archivo_salida.fasta>
```

### Ejemplo

```bash
./alineador sequences.fasta aligned_sequences.fasta
```

### Formato de entrada

```fasta
>Secuencia_1
ATCGATCGATCG
>Secuencia_2
ATCGAACGATCG
>Secuencia_3
ATCGATCGAACG
```

### Formato de salida

```fasta
>Secuencia_1
ATCGATCGATCG
>Secuencia_2
ATCGAACGATCG
>Secuencia_3
ATCGATCGAACG
```

## 🔬 Algoritmo

El alineador implementa un enfoque progresivo en tres pasos:

1. **Matriz de Distancias**  
2. **Árbol Guía (UPGMA)**  
3. **Alineamiento Progresivo (Needleman-Wunsch)**  

## 📊 Parámetros de Alineamiento

| Parámetro         | Valor | Descripción                         |
|------------------|-------|-------------------------------------|
| Match Score      | +2    | Puntuación por coincidencia         |
| Mismatch Score   | -1    | Penalización por no coincidencia    |
| Gap Penalty      | -2    | Penalización por abrir un gap       |
| Gap Extension    | -1    | Penalización por extender un gap    |

## 📈 Estadísticas de Salida

- Tiempo de ejecución total  
- Número de secuencias procesadas  
- Longitud del alineamiento  
- Total de gaps insertados  
- Porcentaje de gaps  
- Representación del árbol guía  

## 🔧 Personalización

Puedes modificar parámetros de puntuación en `alignment.cpp`, y cambiar el alfabeto (ADN o proteínas) en `alignment.h`.

## 🧪 Casos de Prueba

Crea un archivo `test.fasta` con varias secuencias, luego ejecuta:

```bash
./alineador test.fasta result.fasta
```

## ⚠️ Limitaciones

- Alto costo computacional para muchos datos  
- Uso de memoria elevado (sin optimizaciones avanzadas)  
- Modelo de penalización lineal  
- No considera estructura secundaria  

## 🔮 Mejoras Futuras

- [ ] Integración de algoritmos más avanzados (T-Coffee, MUSCLE)  
- [ ] Soporte para modelos HMM  
- [ ] Paralelización con OpenMP  
- [ ] Interfaz gráfica  
- [ ] Exportación en múltiples formatos (Clustal, Phylip)  
- [ ] Alineamiento en streaming para secuencias largas  

## 📝 Licencia

Este proyecto está bajo licencia MIT. Ver archivo LICENSE.

## 🤝 Contribuciones

Las contribuciones son bienvenidas. Por favor:

1. Fork del repositorio  
2. Crea una rama (`feature/NuevaFuncion`)  
3. Haz commit (`git commit -m 'Agrega nueva función'`)  
4. Push a tu rama  
5. Abre un Pull Request  

## 📊 Evaluación y Benchmarking

### Sistema de Benchmarking Integrado

El proyecto incluye un framework completo de evaluación:

```bash
# Compilar sistema de benchmarks
g++ -std=c++17 -O3 -Wall -Wextra src/benchmark_main.cpp src/benchmark.cpp src/alignment.cpp src/io.cpp -o benchmark

# Ejecutar benchmarks individuales
./benchmark single dataset.fasta
./benchmark scalability base.fasta 100 10
./benchmark synthetic 25 200 0.05 output.fasta

# Script automatizado Python
python3 scripts/run_benchmarks.py --all
python3 scripts/run_benchmarks.py --category small
```

### Casos de Uso Evaluados

| Tipo | Descripción | Rendimiento |
|------|-------------|-------------|
| **ADN Conservado** | Genes ortólogos entre especies | ✅ Excelente |
| **ADN Divergente** | Familias génicas con variación | ✅ Muy bueno |
| **Proteínas** | Enzimas y proteínas funcionales | ✅ Óptimo |
| **Escalabilidad** | 10-10,000 secuencias | ⚠️ Límite ~1000 |

### Métricas de Rendimiento

- **Datasets pequeños** (<100 seqs): <0.1s
- **Datasets medianos** (100-1000 seqs): 0.5-5s  
- **Uso de memoria**: 1-20 MB (escala con dataset)
- **Calidad**: 1-30% gaps (dependiente de similitud)

📋 **Ver análisis completo**: [EVALUACION_Y_RESULTADOS.md](EVALUACION_Y_RESULTADOS.md)

## 📞 Contacto

Para dudas, sugerencias o errores, abre un *issue* en el repositorio.

---

**¡Gracias por usar el Alineador Múltiple de Secuencias!** 🧬✨
