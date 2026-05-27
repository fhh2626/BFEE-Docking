import pathlib
import traceback
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QGroupBox, QDoubleSpinBox, QSpinBox, QProgressBar, QTableWidget,
    QTableWidgetItem, QFileDialog, QCheckBox
)
from PySide6.QtCore import Qt, QThread, Signal

from .. import drug_file_parser

class LigandProcessingWorker(QThread):
    """Worker thread for processing ligand files in background with progress updates."""
    progress_updated = Signal(int, int, str)  # (current, total, step_description)
    processing_finished = Signal(list, list)  # (smiles_list, pdbqt_files)
    processing_error = Signal(str)  # error message
    
    def __init__(self, ligand_files, output_dir, ph, reconsider_protonation=True, max_workers=None):
        super().__init__()
        self.ligand_files = ligand_files
        self.output_dir = output_dir
        self.ph = ph
        self.reconsider_protonation = reconsider_protonation
        self.max_workers = max_workers
    
    def run(self):
        """Run ligand processing using DrugFileParser's public methods."""
        try:
            # Create DrugFileParser
            ligand_parser = drug_file_parser.DrugFileParser(
                input_data=self.ligand_files,
                output_dir=str(self.output_dir)
            )
            
            # Run the processing pipeline
            # The pipeline handles progress updates internally and ensures correct sequence
            if self.reconsider_protonation:
                ligand_parser.process_ligand_pipeline(
                    ph=self.ph,
                    progress_callback=self.progress_updated.emit,
                    max_workers=self.max_workers
                )
            else:
                ligand_parser.process_ligand_naive_pipeline(
                    progress_callback=self.progress_updated.emit,
                    max_workers=self.max_workers
                )
            
            # Get results
            smiles_list = ligand_parser.get_generated_smiles()
            pdbqt_files = [str(p) for p in ligand_parser.get_generated_pdbqt_paths()]
            
            # Emit completion
            self.processing_finished.emit(smiles_list, pdbqt_files)
            
        except Exception as e:
            traceback.print_exc()
            self.processing_error.emit(str(e))

class LigandTab(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self._setup_ui()
        
    def _setup_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(15)
        
        # Create GroupBox for Process Ligand
        process_group = QGroupBox("Process Ligand")
        layout = QVBoxLayout(process_group)
        layout.setSpacing(15)
        
        # First row: Ligand file selection
        ligand_file_layout = QHBoxLayout()
        ligand_file_layout.addWidget(QLabel("Ligand:"))
        
        self._ligand_file_input = QLineEdit()
        self._ligand_file_input.setPlaceholderText("Select ligand file(s)...")
        self._ligand_file_input.setReadOnly(True)
        ligand_file_layout.addWidget(self._ligand_file_input)
        
        browse_ligand_btn = QPushButton("Browse")
        browse_ligand_btn.clicked.connect(self._browse_ligand_files)
        ligand_file_layout.addWidget(browse_ligand_btn)
        
        layout.addLayout(ligand_file_layout)
        
        # Second row: pH and Process button
        ph_process_layout = QHBoxLayout()
        
        self._reconsider_protonation_chk = QCheckBox("Reconsider Protonation States")
        self._reconsider_protonation_chk.setChecked(True)
        ph_process_layout.addWidget(self._reconsider_protonation_chk)
        
        ph_process_layout.addWidget(QLabel("pH:"))
        
        self._ligand_ph_spinbox = QDoubleSpinBox()
        self._ligand_ph_spinbox.setRange(0.0, 14.0)
        self._ligand_ph_spinbox.setValue(7.0)
        self._ligand_ph_spinbox.setSingleStep(0.1)
        self._ligand_ph_spinbox.setDecimals(1)
        ph_process_layout.addWidget(self._ligand_ph_spinbox)
        
        # Connect checkbox to spinbox enabled state
        self._reconsider_protonation_chk.toggled.connect(self._ligand_ph_spinbox.setEnabled)
        
        ph_process_layout.addWidget(QLabel("Max Processes:"))
        
        self._max_processes_spinbox = QSpinBox()
        self._max_processes_spinbox.setRange(1, 64)
        self._max_processes_spinbox.setValue(4)
        self._max_processes_spinbox.setToolTip("Maximum number of parallel processes for PDBQT conversion")
        ph_process_layout.addWidget(self._max_processes_spinbox)
        
        ph_process_layout.addStretch()
        
        self._process_ligand_btn = QPushButton("Process")
        self._process_ligand_btn.clicked.connect(self.process_ligand)
        ph_process_layout.addWidget(self._process_ligand_btn)
        
        layout.addLayout(ph_process_layout)
        
        # Third row: Progress bar
        self._ligand_progress = QProgressBar()
        self._ligand_progress.setMinimum(0)
        self._ligand_progress.setMaximum(100)  # Will be updated dynamically based on number of molecules
        self._ligand_progress.setValue(0)
        self._ligand_progress.setFormat("%p%")
        layout.addWidget(self._ligand_progress)
        
        # Add GroupBox to main layout
        main_layout.addWidget(process_group)
        
        # Create GroupBox for Ligands table
        ligands_group = QGroupBox("Ligands")
        ligands_layout = QVBoxLayout(ligands_group)
        
        # Create table widget
        self._ligands_table = QTableWidget()
        self._ligands_table.setColumnCount(3)
        self._ligands_table.setHorizontalHeaderLabels(["Number", "SMILES", "PDBQT"])
        
        # Set table properties
        self._ligands_table.horizontalHeader().setStretchLastSection(True)
        self._ligands_table.setSelectionBehavior(QTableWidget.SelectRows)
        self._ligands_table.setEditTriggers(QTableWidget.NoEditTriggers)
        
        ligands_layout.addWidget(self._ligands_table)
        main_layout.addWidget(ligands_group)
        
        main_layout.addStretch()

    def _browse_ligand_files(self):
        """Browse for ligand files (multiple selection allowed)."""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Select Ligand File(s)",
            self.main_window._last_dir,
            "Ligand Files (*.pdb *.mol2 *.sdf *.smi);;All Files (*.*)"
        )
        
        if file_paths:
            self.main_window._ligand_files = file_paths
            # Update last used directory
            self.main_window._last_dir = str(pathlib.Path(file_paths[0]).parent)
            # Display file names in the input field
            file_names = [pathlib.Path(f).name for f in file_paths]
            if len(file_names) <= 3:
                display_text = ", ".join(file_names)
            else:
                display_text = f"{', '.join(file_names[:3])}, ... ({len(file_names)} files)"
            self._ligand_file_input.setText(display_text)

    def process_ligand(self):
        """Process ligand files with DrugFileParser in a background thread."""
        try:
            # Validate output directory
            if not self.main_window._output_directory:
                self.main_window._show_warning("Error", "Please select an output directory first.")
                return
            
            # Validate ligand files
            if not self.main_window._ligand_files:
                self.main_window._show_warning("Error", "Please select ligand file(s).")
                return
            
            # Get pH
            ph = self._ligand_ph_spinbox.value()
            
            # Prepare output directory
            output_dir = pathlib.Path(self.main_window._output_directory)
            ligand_dir = output_dir / "ligand"
            ligand_dir.mkdir(parents=True, exist_ok=True)
            
            # Reset progress bar
            self._ligand_progress.setValue(0)
            self._ligand_progress.setFormat("0% - Initializing...")
            
            # Disable the process button during processing
            self._process_ligand_btn.setEnabled(False)
            
            # Get max workers
            max_workers = self._max_processes_spinbox.value()
            
            # Get reconsideration status
            reconsider_protonation = self._reconsider_protonation_chk.isChecked()
            
            # Create and start worker thread
            print("Starting ligand processing...")
            self.main_window._ligand_processing_worker = LigandProcessingWorker(
                ligand_files=self.main_window._ligand_files,
                output_dir=ligand_dir,
                ph=ph,
                reconsider_protonation=reconsider_protonation,
                max_workers=max_workers
            )
            
            # Connect signals
            self.main_window._ligand_processing_worker.progress_updated.connect(self._on_ligand_progress_updated)
            self.main_window._ligand_processing_worker.processing_finished.connect(self._on_ligand_processing_finished)
            self.main_window._ligand_processing_worker.processing_error.connect(self._on_ligand_processing_error)
            
            # Start the worker
            self.main_window._ligand_processing_worker.start()
        
        except Exception as e:
            self.main_window._show_error("Error", f"Error starting ligand processing:\n{str(e)}")
            print(f"Error: {e}")
            traceback.print_exc()
            self._process_ligand_btn.setEnabled(True)

    def _on_ligand_progress_updated(self, current, total, description):
        """Handle ligand processing progress updates."""
        self._ligand_progress.setMaximum(total)
        self._ligand_progress.setValue(current)
        # Calculate percentage
        percentage = int((current / total * 100)) if total > 0 else 0
        self._ligand_progress.setFormat(f"{percentage}% - {description}")
        print(f"Progress: {description} ({current}/{total}, {percentage}%)")
    
    def _on_ligand_processing_finished(self, smiles_list, pdbqt_files):
        """Handle ligand processing completion."""
        # Store the results in member variables for later use in docking
        self.main_window._ligand_smiles_list = smiles_list
        self.main_window._ligand_pdbqt_files = pdbqt_files
        
        # Update ligands table
        self._update_ligands_table(smiles_list, pdbqt_files)
        
        # Re-enable the process button
        self._process_ligand_btn.setEnabled(True)
        
        # Get output directory
        output_dir = pathlib.Path(self.main_window._output_directory)
        ligand_dir = output_dir / "ligand"
        
        print(f"✓ Generated {len(smiles_list)} protonated SMILES")
        print(f"✓ Generated {len(pdbqt_files)} PDBQT files")
        
        self.main_window._show_info("Success",
            f"Ligand processing completed successfully!\n\n"
            f"Input files: {len(self.main_window._ligand_files)}\n"
            f"Generated SMILES: {len(smiles_list)}\n"
            f"Generated PDBQT files: {len(pdbqt_files)}\n\n"
            f"Output directory: {ligand_dir}")
    
    def _on_ligand_processing_error(self, error_msg):
        """Handle ligand processing errors."""
        self.main_window._show_error("Error", f"Error processing ligand:\n{error_msg}")
        print(f"Error: {error_msg}")
        # Re-enable the process button
        self._process_ligand_btn.setEnabled(True)
        # Reset progress bar
        self._ligand_progress.setValue(0)

    def _update_ligands_table(self, smiles_list, pdbqt_files):
        """Update the ligands table with SMILES and PDBQT file information."""
        # Clear existing rows
        self._ligands_table.setRowCount(0)
        
        # Ensure we have data to display
        if not pdbqt_files:
            return
        
        # Set row count to match number of PDBQT files
        self._ligands_table.setRowCount(len(pdbqt_files))
        
        # Populate table
        for idx, pdbqt_path in enumerate(pdbqt_files):
            # Number column
            number_item = QTableWidgetItem(str(idx + 1))
            number_item.setTextAlignment(Qt.AlignCenter)
            self._ligands_table.setItem(idx, 0, number_item)
            
            # SMILES column (if available)
            if idx < len(smiles_list):
                smiles_item = QTableWidgetItem(smiles_list[idx])
                self._ligands_table.setItem(idx, 1, smiles_item)
            else:
                self._ligands_table.setItem(idx, 1, QTableWidgetItem("N/A"))
            
            # PDBQT column (filename only)
            pdbqt_filename = pathlib.Path(pdbqt_path).name
            pdbqt_item = QTableWidgetItem(pdbqt_filename)
            self._ligands_table.setItem(idx, 2, pdbqt_item)
        
        # Resize columns to contents
        self._ligands_table.resizeColumnsToContents()
