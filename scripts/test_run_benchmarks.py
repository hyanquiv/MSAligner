#!/usr/bin/env python3
"""
Test suite for run_benchmarks.py
Provides comprehensive coverage for the MSABenchmarkRunner class
"""

import unittest
import tempfile
import shutil
import json
import csv
from pathlib import Path
from unittest.mock import Mock, patch, mock_open
import subprocess

# Import the module to test
import run_benchmarks


class TestMSABenchmarkRunner(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.runner = run_benchmarks.MSABenchmarkRunner("./test_executable")
        self.runner.benchmark_dir = self.temp_dir
        self.runner.results_dir = self.temp_dir / "results"
        self.runner.datasets_dir = self.temp_dir / "datasets"
        
        # Create test directory structure
        self.runner.results_dir.mkdir(parents=True, exist_ok=True)
        self.runner.datasets_dir.mkdir(parents=True, exist_ok=True)
        
        # Create test datasets
        for category in ["small", "medium", "large"]:
            cat_dir = self.runner.datasets_dir / category
            cat_dir.mkdir(exist_ok=True)
            for i in range(2):
                test_file = cat_dir / f"test{i}.fasta"
                test_file.write_text(">seq1\nATCG\n>seq2\nGCTA\n")
    
    def tearDown(self):
        """Clean up test fixtures"""
        shutil.rmtree(self.temp_dir)
    
    def test_find_datasets_with_category(self):
        """Test finding datasets for a specific category"""
        datasets = self.runner.find_datasets("small")
        self.assertEqual(len(datasets), 2)
        self.assertTrue(all("small" in str(d) for d in datasets))
    
    def test_find_datasets_all_categories(self):
        """Test finding datasets for all categories"""
        datasets = self.runner.find_datasets()
        self.assertEqual(len(datasets), 6)  # 2 files per 3 categories
    
    def test_find_datasets_nonexistent_category(self):
        """Test finding datasets for nonexistent category"""
        datasets = self.runner.find_datasets("nonexistent")
        self.assertEqual(len(datasets), 0)
    
    @patch('subprocess.run')
    def test_run_single_benchmark_success(self, mock_run):
        """Test successful single benchmark run"""
        mock_result = Mock()
        mock_result.returncode = 0
        mock_result.stdout = "Tiempo total: 1.5s\nSecuencias procesadas: 10"
        mock_result.stderr = ""
        mock_run.return_value = mock_result
        
        test_dataset = str(self.runner.datasets_dir / "small" / "test0.fasta")
        result = self.runner.run_single_benchmark(test_dataset)
        
        self.assertTrue(result["success"])
        self.assertEqual(result["dataset"], test_dataset)
        self.assertIn("output", result)
        self.assertIn("timestamp", result)
    
    @patch('subprocess.run')
    def test_run_single_benchmark_failure(self, mock_run):
        """Test failed single benchmark run"""
        mock_result = Mock()
        mock_result.returncode = 1
        mock_result.stderr = "Error: Invalid input"
        mock_run.return_value = mock_result
        
        test_dataset = str(self.runner.datasets_dir / "small" / "test0.fasta")
        result = self.runner.run_single_benchmark(test_dataset)
        
        self.assertFalse(result["success"])
        self.assertEqual(result["error"], "Error: Invalid input")
    
    @patch('subprocess.run')
    def test_run_single_benchmark_timeout(self, mock_run):
        """Test benchmark timeout"""
        mock_run.side_effect = subprocess.TimeoutExpired("cmd", 300)
        
        test_dataset = str(self.runner.datasets_dir / "small" / "test0.fasta")
        result = self.runner.run_single_benchmark(test_dataset)
        
        self.assertFalse(result["success"])
        self.assertEqual(result["error"], "Timeout")
    
    @patch('subprocess.run')
    def test_run_single_benchmark_exception(self, mock_run):
        """Test benchmark with general exception"""
        mock_run.side_effect = Exception("General error")
        
        test_dataset = str(self.runner.datasets_dir / "small" / "test0.fasta")
        result = self.runner.run_single_benchmark(test_dataset)
        
        self.assertFalse(result["success"])
        self.assertEqual(result["error"], "General error")
    
    @patch.object(run_benchmarks.MSABenchmarkRunner, 'run_single_benchmark')
    def test_run_multiple_benchmarks(self, mock_single):
        """Test running multiple benchmarks"""
        mock_single.return_value = {
            "dataset": "test.fasta",
            "success": True,
            "timestamp": "20231201_120000"
        }
        
        datasets = ["test1.fasta", "test2.fasta"]
        results, output_dir = self.runner.run_multiple_benchmarks(datasets)
        
        self.assertEqual(len(results), 2)
        self.assertEqual(mock_single.call_count, 2)
        self.assertTrue(output_dir.exists())
    
    def test_generate_text_report_with_mixed_results(self):
        """Test text report generation with mixed success/failure results"""
        results = [
            {
                "dataset": "test1.fasta",
                "success": True,
                "stdout": "Tiempo total: 1.5s\nSecuencias procesadas: 10"
            },
            {
                "dataset": "test2.fasta", 
                "success": False,
                "error": "Memory error"
            }
        ]
        
        report_file = self.temp_dir / "test_report.txt"
        self.runner._generate_text_report(results, report_file)
        
        self.assertTrue(report_file.exists())
        content = report_file.read_text()
        self.assertIn("REPORTE DE BENCHMARKS", content)
        self.assertIn("Exitosos: 1", content)
        self.assertIn("Fallidos: 1", content)
        self.assertIn("test2.fasta: Memory error", content)
    
    def test_write_report_header(self):
        """Test report header writing"""
        results = [{"test": "data"}, {"test": "data2"}]
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            self.runner._write_report_header(f, results)
            temp_file = f.name
        
        with open(temp_file, 'r') as f:
            content = f.read()
        
        self.assertIn("REPORTE DE BENCHMARKS - MSA Aligner", content)
        self.assertIn("Total de benchmarks: 2", content)
        Path(temp_file).unlink()
    
    def test_write_summary_stats(self):
        """Test summary statistics writing"""
        successful = [{"success": True}] * 3
        failed = [{"success": False}] * 2
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            self.runner._write_summary_stats(f, successful, failed)
            temp_file = f.name
        
        with open(temp_file, 'r') as f:
            content = f.read()
        
        self.assertIn("Exitosos: 3", content)
        self.assertIn("Fallidos: 2", content)
        Path(temp_file).unlink()
    
    def test_write_failed_benchmarks_empty(self):
        """Test writing failed benchmarks section when no failures"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            self.runner._write_failed_benchmarks(f, [])
            temp_file = f.name
        
        with open(temp_file, 'r') as f:
            content = f.read()
        
        self.assertEqual(content, "")  # Should write nothing
        Path(temp_file).unlink()
    
    def test_write_failed_benchmarks_with_failures(self):
        """Test writing failed benchmarks section with failures"""
        failed = [
            {"dataset": "test1.fasta", "error": "Error 1"},
            {"dataset": "test2.fasta", "error": "Error 2"}
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            self.runner._write_failed_benchmarks(f, failed)
            temp_file = f.name
        
        with open(temp_file, 'r') as f:
            content = f.read()
        
        self.assertIn("BENCHMARKS FALLIDOS", content)
        self.assertIn("test1.fasta: Error 1", content)
        self.assertIn("test2.fasta: Error 2", content)
        Path(temp_file).unlink()
    
    def test_write_successful_benchmarks(self):
        """Test writing successful benchmarks section"""
        successful = [
            {"dataset": "test1.fasta"},
            {"dataset": "test2.fasta", "stdout": "Tiempo total: 1.5s"}
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            self.runner._write_successful_benchmarks(f, successful)
            temp_file = f.name
        
        with open(temp_file, 'r') as f:
            content = f.read()
        
        self.assertIn("BENCHMARKS EXITOSOS", content)
        self.assertIn("test1.fasta", content)
        self.assertIn("test2.fasta", content)
        Path(temp_file).unlink()
    
    def test_write_metrics_if_available_no_stdout(self):
        """Test metrics writing when no stdout available"""
        success = {"dataset": "test.fasta"}
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            self.runner._write_metrics_if_available(f, success)
            temp_file = f.name
        
        with open(temp_file, 'r') as f:
            content = f.read()
        
        self.assertEqual(content, "")  # Should write nothing
        Path(temp_file).unlink()
    
    def test_write_metrics_if_available_with_metrics(self):
        """Test metrics writing when metrics are available"""
        success = {
            "dataset": "test.fasta",
            "stdout": "Tiempo total: 1.5s\nSecuencias procesadas: 10\nGaps insertados: 5"
        }
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            self.runner._write_metrics_if_available(f, success)
            temp_file = f.name
        
        with open(temp_file, 'r') as f:
            content = f.read()
        
        self.assertIn("Tiempo total: 1.5s", content)
        self.assertIn("Secuencias procesadas: 10", content)
        self.assertIn("Gaps insertados: 5", content)
        Path(temp_file).unlink()
    
    def test_generate_csv_report(self):
        """Test CSV report generation"""
        results = [
            {
                "dataset": "test1.fasta",
                "success": True,
                "timestamp": "20231201_120000",
                "output": "output1.fasta"
            },
            {
                "dataset": "test2.fasta",
                "success": False,
                "timestamp": "20231201_120001",
                "error": "Test error"
            }
        ]
        
        csv_file = self.temp_dir / "test_report.csv"
        self.runner._generate_csv_report(results, csv_file)
        
        self.assertTrue(csv_file.exists())
        
        with open(csv_file, 'r', newline='') as f:
            reader = csv.reader(f)
            rows = list(reader)
        
        self.assertEqual(len(rows), 3)  # Header + 2 data rows
        self.assertEqual(rows[0], ['Dataset', 'Success', 'Timestamp', 'Output_File', 'Error'])
        self.assertEqual(rows[1][0], 'test1.fasta')
        self.assertEqual(rows[2][0], 'test2.fasta')
    
    def test_generate_json_report(self):
        """Test JSON report generation"""
        results = [
            {"dataset": "test1.fasta", "success": True},
            {"dataset": "test2.fasta", "success": False}
        ]
        
        json_file = self.temp_dir / "test_report.json"
        self.runner._generate_json_report(results, json_file)
        
        self.assertTrue(json_file.exists())
        
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        self.assertEqual(len(data), 2)
        self.assertEqual(data[0]["dataset"], "test1.fasta")
        self.assertEqual(data[1]["dataset"], "test2.fasta")
    
    @patch('builtins.print')
    def test_generate_summary_report_integration(self, mock_print):
        """Test complete summary report generation"""
        results = [
            {"dataset": "test1.fasta", "success": True, "timestamp": "20231201_120000"},
            {"dataset": "test2.fasta", "success": False, "error": "Test error", "timestamp": "20231201_120001"}
        ]
        
        output_dir = self.temp_dir / "test_output"
        output_dir.mkdir()
        
        self.runner.generate_summary_report(results, output_dir)
        
        # Check all report files were created
        self.assertTrue((output_dir / "benchmark_summary.txt").exists())
        self.assertTrue((output_dir / "benchmark_results.csv").exists())
        self.assertTrue((output_dir / "benchmark_results.json").exists())
        
        # Check print was called for file paths
        mock_print.assert_called()


class TestMainFunction(unittest.TestCase):
    """Test the main function and argument parsing"""
    
    @patch('os.path.exists')
    @patch('sys.exit')
    @patch('builtins.print')
    def test_main_executable_not_found(self, mock_print, mock_exit, mock_exists):
        """Test main function when executable doesn't exist"""
        mock_exists.return_value = False
        
        with patch('sys.argv', ['run_benchmarks.py', '--all']):
            run_benchmarks.main()
        
        mock_exit.assert_called_with(1)
        mock_print.assert_called()
    
    @patch('os.path.exists')
    @patch('sys.exit')
    @patch('builtins.print')
    def test_main_no_arguments(self, mock_print, mock_exit, mock_exists):
        """Test main function with no arguments"""
        mock_exists.return_value = True
        
        with patch('sys.argv', ['run_benchmarks.py']):
            run_benchmarks.main()
        
        mock_exit.assert_called_with(1)
    
    @patch('os.path.exists')
    @patch.object(run_benchmarks.MSABenchmarkRunner, 'create_synthetic_datasets')
    def test_main_create_synthetic(self, mock_create, mock_exists):
        """Test main function with create synthetic option"""
        mock_exists.return_value = True
        
        with patch('sys.argv', ['run_benchmarks.py', '--create-synthetic', '10', '20']):
            run_benchmarks.main()
        
        mock_create.assert_called_once()


if __name__ == '__main__':
    unittest.main()