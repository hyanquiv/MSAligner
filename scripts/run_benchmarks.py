#!/usr/bin/env python3
"""
Script para ejecutar benchmarks automatizados del MSA Aligner
Ejecuta múltiples benchmarks y genera reportes detallados
"""

import os
import sys
import subprocess
import json
import csv
import argparse
from datetime import datetime
from pathlib import Path

class MSABenchmarkRunner:
    def __init__(self, executable_path="./alineador"):
        self.executable_path = executable_path
        self.benchmark_dir = Path("benchmarks")
        self.results_dir = self.benchmark_dir / "results"
        self.datasets_dir = self.benchmark_dir / "datasets"
        
        # Crear directorios si no existen
        self.results_dir.mkdir(parents=True, exist_ok=True)
    
    def find_datasets(self, category=None):
        """Encuentra todos los datasets disponibles"""
        datasets = []
        
        if category:
            search_dir = self.datasets_dir / category
            if search_dir.exists():
                datasets.extend(search_dir.glob("*.fasta"))
        else:
            # Buscar en todas las categorías
            for subdir in ["small", "medium", "large"]:
                cat_dir = self.datasets_dir / subdir
                if cat_dir.exists():
                    datasets.extend(cat_dir.glob("*.fasta"))
        
        return [str(d) for d in datasets]
    
    def run_single_benchmark(self, dataset_path, output_dir=None):
        """Ejecuta un benchmark individual"""
        if output_dir is None:
            output_dir = self.results_dir
        
        dataset_name = Path(dataset_path).stem
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = output_dir / f"{dataset_name}_aligned_{timestamp}.fasta"
        
        print(f"Ejecutando benchmark: {dataset_path}")
        print(f"Salida: {output_file}")
        
        try:
            # Ejecutar el alineador y capturar salida
            cmd = [self.executable_path, dataset_path, str(output_file)]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0:
                print("✅ Benchmark completado exitosamente")
                return {
                    "dataset": dataset_path,
                    "output": str(output_file),
                    "stdout": result.stdout,
                    "stderr": result.stderr,
                    "success": True,
                    "timestamp": timestamp
                }
            else:
                print("❌ Error en benchmark")
                print("STDERR:", result.stderr)
                return {
                    "dataset": dataset_path,
                    "success": False,
                    "error": result.stderr,
                    "timestamp": timestamp
                }
                
        except subprocess.TimeoutExpired:
            print("⏰ Timeout - Benchmark cancelado")
            return {
                "dataset": dataset_path,
                "success": False,
                "error": "Timeout",
                "timestamp": timestamp
            }
        except Exception as e:
            print(f"💥 Excepción: {e}")
            return {
                "dataset": dataset_path,
                "success": False,
                "error": str(e),
                "timestamp": timestamp
            }
    
    def run_multiple_benchmarks(self, datasets, output_dir=None):
        """Ejecuta múltiples benchmarks"""
        if output_dir is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_dir = self.results_dir / f"batch_{timestamp}"
            output_dir.mkdir(exist_ok=True)
        
        results = []
        
        print("\\n🚀 Iniciando {} benchmarks...".format(len(datasets)))
        print("=" * 60)
        
        for i, dataset in enumerate(datasets, 1):
            print(f"\\n[{i}/{len(datasets)}] {Path(dataset).name}")
            print("-" * 40)
            
            result = self.run_single_benchmark(dataset, output_dir)
            results.append(result)
        
        return results, output_dir
    
    def generate_summary_report(self, results, output_dir):
        """Genera reporte resumen de benchmarks"""
        report_file = output_dir / "benchmark_summary.txt"
        csv_file = output_dir / "benchmark_results.csv"
        json_file = output_dir / "benchmark_results.json"
        
        self._generate_text_report(results, report_file)
        self._generate_csv_report(results, csv_file)
        self._generate_json_report(results, json_file)
        
        print("\\n📊 Reportes generados:")
        print(f"  📄 Resumen: {report_file}")
        print(f"  📈 CSV: {csv_file}")
        print(f"  🔧 JSON: {json_file}")
    
    def _generate_text_report(self, results, report_file):
        """Genera reporte de texto"""
        successful = [r for r in results if r['success']]
        failed = [r for r in results if not r['success']]
        
        with open(report_file, 'w') as f:
            self._write_report_header(f, results)
            self._write_summary_stats(f, successful, failed)
            self._write_failed_benchmarks(f, failed)
            self._write_successful_benchmarks(f, successful)
    
    def _write_report_header(self, f, results):
        """Escribe el encabezado del reporte"""
        f.write("REPORTE DE BENCHMARKS - MSA Aligner\\n")
        f.write("=" * 50 + "\\n\\n")
        f.write(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
        f.write(f"Total de benchmarks: {len(results)}\\n")
    
    def _write_summary_stats(self, f, successful, failed):
        """Escribe estadísticas de resumen"""
        f.write(f"Exitosos: {len(successful)}\\n")
        f.write(f"Fallidos: {len(failed)}\\n\\n")
    
    def _write_failed_benchmarks(self, f, failed):
        """Escribe sección de benchmarks fallidos"""
        if not failed:
            return
        
        f.write("BENCHMARKS FALLIDOS:\\n")
        f.write("-" * 20 + "\\n")
        for fail in failed:
            f.write(f"  {Path(fail['dataset']).name}: {fail.get('error', 'Error desconocido')}\\n")
        f.write("\\n")
    
    def _write_successful_benchmarks(self, f, successful):
        """Escribe sección de benchmarks exitosos"""
        f.write("BENCHMARKS EXITOSOS:\\n")
        f.write("-" * 20 + "\\n")
        for success in successful:
            f.write(f"  {Path(success['dataset']).name}\\n")
            self._write_metrics_if_available(f, success)
            f.write("\\n")
    
    def _write_metrics_if_available(self, f, success):
        """Escribe métricas si están disponibles en stdout"""
        if 'stdout' not in success:
            return
        
        stdout = success['stdout']
        if "Tiempo total:" not in stdout:
            return
        
        lines = stdout.split('\\n')
        keywords = ["Tiempo total:", "Secuencias procesadas:", "Gaps insertados:"]
        for line in lines:
            if any(keyword in line for keyword in keywords):
                f.write(f"    {line.strip()}\\n")
    
    def _generate_csv_report(self, results, csv_file):
        """Genera reporte CSV"""
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Dataset', 'Success', 'Timestamp', 'Output_File', 'Error'])
            
            for result in results:
                writer.writerow([
                    Path(result['dataset']).name,
                    result['success'],
                    result['timestamp'],
                    result.get('output', ''),
                    result.get('error', '')
                ])
    
    def _generate_json_report(self, results, json_file):
        """Genera reporte JSON"""
        with open(json_file, 'w') as f:
            json.dump(results, f, indent=2)
    
    def run_scalability_test(self, base_dataset, max_sequences=50, step=10):
        """Ejecuta test de escalabilidad"""
        print("\\n🔬 Test de escalabilidad con {}".format(base_dataset))
        print("Desde {} hasta {} secuencias".format(step, max_sequences))
        
        # Este método requeriría integración con el código C++ de benchmark
        # Por ahora, solo muestra el concepto
        print("⚠️  Test de escalabilidad requiere integración con benchmark.cpp")
    
    def create_synthetic_datasets(self, sizes, output_dir=None):
        """Crea datasets sintéticos para pruebas"""
        if output_dir is None:
            output_dir = self.datasets_dir / "synthetic"
            output_dir.mkdir(exist_ok=True)
        
        print("\\n🧪 Creando datasets sintéticos...")
        
        for size in sizes:
            print("  Creando dataset con {} secuencias...".format(size))
            # Este método también requeriría integración con benchmark.cpp
            print("    ⚠️  Requiere implementación en C++")


def main():
    parser = argparse.ArgumentParser(
        description="Ejecutor de benchmarks para MSA Aligner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  python run_benchmarks.py --all                    # Todos los datasets
  python run_benchmarks.py --category small         # Solo datasets pequeños
  python run_benchmarks.py --dataset entrada.fasta  # Dataset específico
  python run_benchmarks.py --scalability            # Test de escalabilidad
        """
    )
    
    parser.add_argument('--executable', '-e', default='./alineador',
                       help='Ruta al ejecutable del alineador')
    parser.add_argument('--all', action='store_true',
                       help='Ejecutar benchmarks en todos los datasets')
    parser.add_argument('--category', '-c', choices=['small', 'medium', 'large'],
                       help='Ejecutar benchmarks en una categoría específica')
    parser.add_argument('--dataset', '-d',
                       help='Ejecutar benchmark en un dataset específico')
    parser.add_argument('--scalability', action='store_true',
                       help='Ejecutar test de escalabilidad')
    parser.add_argument('--create-synthetic', nargs='+', type=int, metavar='SIZE',
                       help='Crear datasets sintéticos con los tamaños especificados')
    parser.add_argument('--output', '-o',
                       help='Directorio de salida personalizado')
    
    args = parser.parse_args()
    
    # Verificar que el ejecutable existe
    if not os.path.exists(args.executable):
        print("❌ Error: No se encuentra el ejecutable {}".format(args.executable))
        print("💡 Asegúrate de compilar el proyecto primero:")
        print("   g++ -std=c++17 -O3 -Wall -Wextra src/MSAligner.cpp src/alignment.cpp src/io.cpp -o alineador")
        sys.exit(1)
    
    runner = MSABenchmarkRunner(args.executable)
    
    if args.create_synthetic:
        runner.create_synthetic_datasets(args.create_synthetic, args.output)
        return
    
    if args.scalability:
        # Buscar un dataset base para escalabilidad
        datasets = runner.find_datasets("small")
        if datasets:
            runner.run_scalability_test(datasets[0])
        else:
            print("❌ No se encontraron datasets para test de escalabilidad")
        return
    
    # Determinar qué datasets ejecutar
    datasets = []
    
    if args.dataset:
        if os.path.exists(args.dataset):
            datasets = [args.dataset]
        else:
            print("❌ Error: Dataset {} no encontrado".format(args.dataset))
            sys.exit(1)
    elif args.category:
        datasets = runner.find_datasets(args.category)
    elif args.all:
        datasets = runner.find_datasets()
    else:
        print("❌ Error: Especifica --all, --category, --dataset o --scalability")
        parser.print_help()
        sys.exit(1)
    
    if not datasets:
        print("❌ No se encontraron datasets para ejecutar")
        sys.exit(1)
    
    print("\\n📋 Datasets encontrados ({}):".format(len(datasets)))
    for dataset in datasets:
        print("  📁 {}".format(dataset))
    
    # Ejecutar benchmarks
    results, output_dir = runner.run_multiple_benchmarks(datasets, args.output)
    
    # Generar reportes
    runner.generate_summary_report(results, output_dir)
    
    print("\\n✅ Benchmarks completados!")
    print("📂 Resultados en: {}".format(output_dir))


if __name__ == "__main__":
    main()