import pathlib
import traceback
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QGroupBox, QDoubleSpinBox, QComboBox, QProgressBar, QFileDialog
)
from PySide6.QtCore import QThread, Signal

from .. import docking

class DockingWorker(QThread):
    """Worker thread for running docking in background with progress updates."""
    progress_updated = Signal(int, int)  # (current, total)
    docking_finished = Signal(list)  # docking results
    docking_error = Signal(str)  # error message
    
    def __init__(self, docking_instance):
        super().__init__()
        self.docking_instance = docking_instance
    
    def run(self):
        """Run docking in background thread using Docking.run_docking."""
        try:
            # Relay progress updates from Docking class to GUI signals
            def progress_callback(current, total):
                self.progress_updated.emit(current, total)

            # Execution
            # Get results
            # The docking instance stores results internally, but run_docking also returns them.
            results = self.docking_instance.run_docking(progress_callback=progress_callback)
            
            # Emit completion
            self.docking_finished.emit(results)
            
        except Exception as e:
            traceback.print_exc()
            self.docking_error.emit(str(e))

class DockingTab(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self._setup_ui()
        
    def _setup_ui(self):
        main_layout = QVBoxLayout(self)
        
        # Create GroupBox for Docking
        docking_group = QGroupBox("Docking")
        layout = QVBoxLayout(docking_group)
        
        # Single row: All docking parameters centered
        # Row 1: Basic numeric parameters
        row1_layout = QHBoxLayout()
        row1_layout.addStretch()
        
        row1_layout.addWidget(QLabel("Exhaustiveness:"))
        self._exhaustiveness_input = QLineEdit()
        self._exhaustiveness_input.setText("32")
        self._exhaustiveness_input.setMaximumWidth(80)
        row1_layout.addWidget(self._exhaustiveness_input)
        
        row1_layout.addWidget(QLabel("Number of Modes:"))
        self._num_modes_input = QDoubleSpinBox()
        self._num_modes_input.setRange(1, 100)
        self._num_modes_input.setValue(10)
        self._num_modes_input.setDecimals(0)
        self._num_modes_input.setMaximumWidth(80)
        row1_layout.addWidget(self._num_modes_input)
        
        row1_layout.addWidget(QLabel("Energy Range:"))
        self._energy_range_input = QDoubleSpinBox()
        self._energy_range_input.setRange(1.0, 10.0)
        self._energy_range_input.setValue(3.0)
        self._energy_range_input.setSingleStep(0.1)
        self._energy_range_input.setDecimals(1)
        self._energy_range_input.setMaximumWidth(80)
        row1_layout.addWidget(self._energy_range_input)
        
        row1_layout.addStretch()
        layout.addLayout(row1_layout)
        
        # Row 2: Engine and Scoring
        row2_layout = QHBoxLayout()
        row2_layout.addStretch()
        
        row2_layout.addWidget(QLabel("Engine:"))
        self._engine_combo = QComboBox()
        self._engine_combo.addItems(["smina", "vina-new", "qvina2", "qvinaw", "vina-classic", "DSDP", "gnina"])
        self._engine_combo.setMaximumWidth(120)
        self._engine_combo.currentTextChanged.connect(self._on_engine_changed)
        row2_layout.addWidget(self._engine_combo)
        
        row2_layout.addWidget(QLabel("Scoring Function:"))
        self._scoring_combo = QComboBox()
        self._scoring_combo.addItems(["vina", "vinardo"])
        self._scoring_combo.setMaximumWidth(120)
        row2_layout.addWidget(self._scoring_combo)
        
        row2_layout.addWidget(QLabel("Flexible Residues:"))
        self._flexible_residues_input = QLineEdit()
        self._flexible_residues_input.setPlaceholderText("e.g. A:45,A:82,B:120")
        self._flexible_residues_input.setMaximumWidth(180)
        row2_layout.addWidget(self._flexible_residues_input)
        
        row2_layout.addStretch()
        layout.addLayout(row2_layout)
        
        # Row 3: Gnina-specific options (CNN mode)
        row3_layout = QHBoxLayout()
        row3_layout.addStretch()
        
        row3_layout.addWidget(QLabel("CNN Mode:"))
        self._cnn_mode_combo = QComboBox()
        self._cnn_mode_combo.addItems(["none", "rescore", "refinement", "all"])
        self._cnn_mode_combo.setCurrentText("rescore")  # Default to rescore
        self._cnn_mode_combo.setMaximumWidth(120)
        self._cnn_mode_combo.setEnabled(False)  # Disabled by default, enabled only for gnina
        row3_layout.addWidget(self._cnn_mode_combo)
        
        row3_layout.addStretch()
        layout.addLayout(row3_layout)
        
        # Button row: Start Docking button
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        
        self._start_docking_btn = QPushButton("Start Docking")
        self._start_docking_btn.clicked.connect(self.start_docking)
        button_layout.addWidget(self._start_docking_btn)
        
        self._load_results_btn = QPushButton("Load Results (from server)")
        self._load_results_btn.clicked.connect(self.load_results)
        button_layout.addWidget(self._load_results_btn)
        
        button_layout.addStretch()
        layout.addLayout(button_layout)
        
        # Progress bar row
        self._docking_progress = QProgressBar()
        self._docking_progress.setMinimum(0)
        self._docking_progress.setMaximum(100)
        self._docking_progress.setValue(0)
        self._docking_progress.setFormat("0% (0/0 ligands)")
        layout.addWidget(self._docking_progress)
        
        # Add GroupBox to main layout
        main_layout.addWidget(docking_group)
        main_layout.addStretch()
        
        # Initialize the enabled/disabled state based on default engine
        self._on_engine_changed(self._engine_combo.currentText())

    def _on_engine_changed(self, engine_name):
        """Handle engine selection changes."""
        if engine_name in ["vina-new", "smina", "gnina"]:
            # vina-new, smina, and gnina support custom scoring functions
            self._scoring_combo.setEnabled(True)
        else:
            # For other engines (qvina, vina-classic, DSDP), force 'vina' and disable selection
            index = self._scoring_combo.findText("vina")
            if index >= 0:
                self._scoring_combo.setCurrentIndex(index)
            self._scoring_combo.setEnabled(False)
        
        # Disable Energy Range for DSDP and gnina
        if engine_name in ["DSDP", "gnina"]:
            self._energy_range_input.setEnabled(False)
        else:
            self._energy_range_input.setEnabled(True)
        
        # Only smina and gnina support flexible residues
        if engine_name in ["smina", "gnina"]:
            self._flexible_residues_input.setEnabled(True)
        else:
            self._flexible_residues_input.setEnabled(False)
            
        # Only enable "Load Results (No Run)" button when engine is DSDP or gnina
        if engine_name in ["DSDP", "gnina"]:
            self._load_results_btn.setEnabled(True)
            if engine_name == "DSDP":
                self._exhaustiveness_input.setText("384")
            else:
                self._exhaustiveness_input.setText("32")
        else:
            self._load_results_btn.setEnabled(False)
        
        # CNN mode is only available for gnina
        if engine_name == "gnina":
            self._cnn_mode_combo.setEnabled(True)
        else:
            self._cnn_mode_combo.setEnabled(False)
            
    def start_docking(self):
        """Start the molecular docking process."""
        try:
            # Validate output directory
            if not self.main_window._output_directory:
                self.main_window._show_warning("Error", "Please select an output directory first.")
                return
            
            # Validate protein parser and PDBQT
            if not self.main_window._protein_parser:
                self.main_window._show_warning("Error", "Please process the protein first (Protein tab).")
                return
            
            protein_pdbqt = self.main_window._protein_parser.get_generated_pdbqt_path()
            if not protein_pdbqt:
                self.main_window._show_warning("Error", "Protein PDBQT file not found. Please process the protein first.")
                return
            
            # Validate ligand processing results
            if not self.main_window._ligand_pdbqt_files:
                self.main_window._show_warning("Error", "No ligand PDBQT files found. Please process ligands first (Ligand tab).")
                return
            
            if not self.main_window._ligand_smiles_list:
                self.main_window._show_warning("Error", "No ligand SMILES found. Please process ligands first.")
                return
            
            # Validate PDBQT and SMILES counts match
            if len(self.main_window._ligand_pdbqt_files) != len(self.main_window._ligand_smiles_list):
                self.main_window._show_warning("Error", 
                    f"Ligand data mismatch: {len(self.main_window._ligand_pdbqt_files)} PDBQT files but "
                    f"{len(self.main_window._ligand_smiles_list)} SMILES. Please re-process ligands.")
                return
            
            # Validate docking range
            docking_range = self.main_window._protein_parser.get_generated_docking_range()
            if not docking_range:
                self.main_window._show_warning("Error", 
                    "Docking range not set. Please select a docking region mode in Protein tab.")
                return
            
            # Get docking parameters
            try:
                exhaustiveness = int(self._exhaustiveness_input.text())
                num_modes = int(self._num_modes_input.value())
                energy_range = self._energy_range_input.value()
            except ValueError:
                self.main_window._show_warning("Error", 
                    "Invalid docking parameters. Please enter valid numbers.")
                return
            
            scoring = self._scoring_combo.currentText()
            engine = self._engine_combo.currentText()
            
            # Get flexible residues (used for smina and gnina)
            flexible_residues = ""
            if engine in ["smina", "gnina"]:
                flexible_residues = self._flexible_residues_input.text().strip()
            
            # Get CNN scoring mode (only for gnina)
            cnn_scoring = "rescore"  # Default
            if engine == "gnina":
                cnn_scoring = self._cnn_mode_combo.currentText()
            
            # Prepare output directory
            output_dir = pathlib.Path(self.main_window._output_directory)
            result_dir = output_dir / "result"
            result_dir.mkdir(parents=True, exist_ok=True)
            
            # Create Docking object
            print("Creating Docking object...")
            self.main_window._docking_instance = docking.Docking(
                protein_pdbqt=str(protein_pdbqt),
                ligand_pdbqts=[str(p) for p in self.main_window._ligand_pdbqt_files],
                ligand_smiles=self.main_window._ligand_smiles_list,
                docking_range=docking_range,
                exhaustiveness=exhaustiveness,
                num_modes=num_modes,
                scoring=scoring,
                output_dir=str(result_dir),
                energy_range=energy_range,
                engine=engine,
                flexible_residues=flexible_residues,
                cnn_scoring=cnn_scoring
            )
            
            print(f"✓ Docking object created")
            print(f"  Protein: {protein_pdbqt}")
            print(f"  Ligands: {len(self.main_window._ligand_pdbqt_files)}")
            print(f"  Exhaustiveness: {exhaustiveness}")
            print(f"  Num modes: {num_modes}")
            print(f"  Scoring: {scoring}")
            print(f"  Energy range: {energy_range}")
            print(f"  Engine: {engine}")
            if engine in ["smina", "gnina"] and flexible_residues:
                print(f"  Flexible residues: {flexible_residues}")
            if engine == "gnina":
                print(f"  CNN mode: {cnn_scoring}")
            print(f"  Output: {result_dir}")
            
            # Special handling for DSDP engine
            if engine == "DSDP":
                job_path = self.main_window._docking_instance.generate_dsdp_job_package()
                self.main_window._show_info(
                    "DSDP Job Generated",
                    f"Job package generated at:\n{job_path}\n\n"
                    "1. Copy this folder to your Linux machine.\n"
                    "2. Run 'sh run_dsdp.sh' inside the folder.\n"
                    "3. When finished, copy the folder back to Windows.\n"
                    "4. Use the 'Load Results' button in this tab to load the results."
                )
                return

            # Special handling for gnina engine
            if engine == "gnina":
                job_path = self.main_window._docking_instance.generate_gnina_job_package()
                self.main_window._show_info(
                    "Gnina Job Generated",
                    f"Job package generated at:\n{job_path}\n\n"
                    "1. Copy this folder to your Linux machine.\n"
                    "2. Run 'sh run_gnina.sh' inside the folder.\n"
                    "3. When finished, copy the folder back to Windows.\n"
                    "4. Use the 'Load Results' button in this tab to load the results."
                )
                return

            # Initialize progress bar
            total_ligands = len(self.main_window._ligand_pdbqt_files)
            self._docking_progress.setMaximum(total_ligands)
            self._docking_progress.setValue(0)
            self._docking_progress.setFormat(f"0% (0/{total_ligands} ligands)")
            
            # Disable Start Docking button during docking
            self._start_docking_btn.setEnabled(False)
            
            # Create and configure worker thread
            print("\nStarting docking process...")
            self.main_window._docking_worker = DockingWorker(self.main_window._docking_instance)
            self.main_window._docking_worker.progress_updated.connect(self._on_docking_progress)
            self.main_window._docking_worker.docking_finished.connect(self._on_docking_finished)
            self.main_window._docking_worker.docking_error.connect(self._on_docking_error)
            self.main_window._docking_worker.start()
        
        except Exception as e:
            self.main_window._show_error("Error", f"Error setting up docking:\n{str(e)}")
            traceback.print_exc()

    def _on_docking_progress(self, current, total):
        """Handle docking progress updates."""
        self._docking_progress.setValue(current)
        percentage = int((current / total) * 100) if total > 0 else 0
        self._docking_progress.setFormat(f"{percentage}% ({current}/{total} ligands)")
        print(f"Progress: {current}/{total} ligands ({percentage}%)")
    
    def _on_docking_finished(self, results):
        """Handle docking completion."""
        print(f"✓ Docking completed!")
        
        # Re-enable Start Docking button
        self._start_docking_btn.setEnabled(True)
        
        # Update results table - delegating to ResultsTab
        if hasattr(self.main_window, '_results_tab'):
             self.main_window._results_tab.update_results_table()
        
        # Count successful docking results and report failures explicitly.
        successful_docks = len([r for r in results if r.get('success', False)])
        failed_docks = [r for r in results if not r.get('success', False)]
        result_dir = pathlib.Path(self.main_window._output_directory) / 'result'

        if failed_docks:
            failure_lines = []
            for result in failed_docks[:3]:
                ligand_name = pathlib.Path(result.get('ligand', 'ligand')).name
                diagnostic = result.get('stderr') or result.get('error') or ''
                diagnostic_lines = [
                    line.strip() for line in str(diagnostic).splitlines()
                    if line.strip()
                ]
                reason = diagnostic_lines[-1] if diagnostic_lines else (
                    f"return code {result.get('return_code', 'unknown')}"
                )
                failure_lines.append(f"{ligand_name}: {reason}")

            more = "" if len(failed_docks) <= 3 else "\n..."
            self.main_window._show_warning(
                "Docking completed with errors",
                f"Docking finished with failures.\n\n"
                f"Total ligands: {len(results)}\n"
                f"Successful dockings: {successful_docks}\n"
                f"Failed dockings: {len(failed_docks)}\n\n"
                f"Diagnostics:\n{chr(10).join(failure_lines)}{more}\n\n"
                f"Full logs: {result_dir}"
            )
        else:
            self.main_window._show_info(
                "Success",
                f"Docking completed successfully!\n\n"
                f"Total ligands: {len(results)}\n"
                f"Successful dockings: {successful_docks}\n"
                f"Output directory: {result_dir}"
            )
    
    def _on_docking_error(self, error_msg):
        """Handle docking errors."""
        print(f"Error during docking: {error_msg}")
        
        # Re-enable Start Docking button
        self._start_docking_btn.setEnabled(True)
        
        self.main_window._show_error("Error", f"Error during docking:\n{error_msg}")

    def load_results(self):
        """Load docking results from an existing run (e.g. DSDP from Linux)."""
        try:
             # Validate prerequisites matching start_docking
            if not self.main_window._protein_parser or not self.main_window._ligand_pdbqt_files:
                 self.main_window._show_warning("Error", "Please make sure Protein and Ligands are processed and loaded first.")
                 return

            # Select results directory
            result_dir = QFileDialog.getExistingDirectory(
                self, 
                "Select Docking Results Directory (containing index.csv)",
                self.main_window._last_dir
            )
            
            if not result_dir:
                return
                
            self.main_window._last_dir = result_dir
            
            if not (pathlib.Path(result_dir) / "index.csv").exists():
                self.main_window._show_warning("Error", "Please select the folder containing index.csv, usually the dsdp_job/gnina_job folder.")
                return

            # Re-initialize docking instance if needed, or update if exists.
            # We strictly need to match the current GUI state (protein/ligands) to the loaded results.
            # So we create a new docking instance with current GUI parameters but dummy engine settings
            # just to use the loading capability.
            
            # We need minimal valid parameters to init Docking
            try:
                # Try to read docking config first to get engine and flexible residues
                import json
                config_path = pathlib.Path(result_dir) / "docking_config.json"
                engine = "vina-classic"  # Safe default
                flexible_residues = ""
                scoring = "vina"
                cnn_scoring = "rescore"  # Default CNN scoring mode
                
                if config_path.exists():
                    try:
                        with open(config_path, "r", encoding="utf-8") as f:
                            config_data = json.load(f)
                        engine = config_data.get("engine", "vina-classic")
                        flexible_residues = config_data.get("flexible_residues", "")
                        scoring = config_data.get("scoring", "vina")
                        cnn_scoring = config_data.get("cnn_scoring", "rescore")
                        print(f"Loaded config: engine={engine}, flexible_residues={flexible_residues}, scoring={scoring}, cnn_scoring={cnn_scoring}")
                    except Exception as e:
                        print(f"Warning: Failed to read config: {e}")
                
                # Use current settings from GUI
                protein_pdbqt = self.main_window._protein_parser.get_generated_pdbqt_path()
                docking_range = self.main_window._protein_parser.get_generated_docking_range()
                exhaustiveness = int(self._exhaustiveness_input.text())
                num_modes = int(self._num_modes_input.value())
                energy_range = self._energy_range_input.value()
                
                # Use smina as engine if flexible residues are present (supports flex loading)
                # Note: We need an engine that supports flexible_residues parameter
                if flexible_residues and engine not in ["smina", "gnina"]:
                    engine = "smina"  # Use smina for flexible residue support
                
                self.main_window._docking_instance = docking.Docking(
                    protein_pdbqt=str(protein_pdbqt),
                    ligand_pdbqts=[str(p) for p in self.main_window._ligand_pdbqt_files],
                    ligand_smiles=self.main_window._ligand_smiles_list,
                    docking_range=docking_range,
                    exhaustiveness=exhaustiveness,
                    num_modes=num_modes,
                    scoring=scoring,
                    output_dir=result_dir, # Temporary, won't be used for writing
                    energy_range=energy_range,
                    engine=engine,
                    flexible_residues=flexible_residues,
                    cnn_scoring=cnn_scoring
                )
                
                # Load results
                results = self.main_window._docking_instance.load_dsdp_results(result_dir)
                
                # Update UI
                self._on_docking_finished(results)
                
            except Exception as e:
                self.main_window._show_error("Error", f"Failed to load results:\n{str(e)}")
                traceback.print_exc()
                
        except Exception as e:
            traceback.print_exc()
