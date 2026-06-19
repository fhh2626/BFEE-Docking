import pathlib
import traceback
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QGroupBox, QDoubleSpinBox, QSpinBox, QTableWidget, QTableWidgetItem,
    QFileDialog
)
from PySide6.QtCore import Qt

class ResultsTab(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self._setup_ui()
        
    def _setup_ui(self):
        main_layout = QVBoxLayout(self)
        
        # Results GroupBox
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout(results_group)
        
        # Create results table
        self._results_table = QTableWidget()
        self._results_table.setColumnCount(6)
        self._results_table.setHorizontalHeaderLabels([
            "Number", "SMILES", "Best Affinity", "Affinity", 
            "Molecule Number", "Filename"
        ])
        
        # Set table properties
        self._results_table.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self._results_table.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self._results_table.setSelectionBehavior(QTableWidget.SelectRows)
        self._results_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self._results_table.cellClicked.connect(self._on_table_row_clicked)
        
        results_layout.addWidget(self._results_table)
        main_layout.addWidget(results_group)
        
        # Visualization GroupBox
        viz_group = QGroupBox("Visualization")
        viz_layout = QVBoxLayout(viz_group)
        
        # First row: Sort, Cluster, Similarity Threshold
        first_row_layout = QHBoxLayout()
        
        # Sort button
        self._sort_results_btn = QPushButton("Sort")
        self._sort_results_btn.clicked.connect(self._sort_results_table)
        first_row_layout.addWidget(self._sort_results_btn)

        # Cluster button
        self._cluster_btn = QPushButton("Cluster")
        self._cluster_btn.clicked.connect(self._cluster_results)
        first_row_layout.addWidget(self._cluster_btn)
        
        # Similarity Threshold
        first_row_layout.addWidget(QLabel("Similarity Threshold:"))
        self._similarity_threshold_input = QDoubleSpinBox()
        self._similarity_threshold_input.setRange(0.0, 1.0)
        self._similarity_threshold_input.setValue(0.7)
        self._similarity_threshold_input.setSingleStep(0.1)
        self._similarity_threshold_input.setDecimals(2)
        #self._similarity_threshold_input.setMaximumWidth(60)
        first_row_layout.addWidget(self._similarity_threshold_input)
        
        viz_layout.addLayout(first_row_layout)
        
        # Second row: Molecule, Structure, Show in Pymol
        second_row_layout = QHBoxLayout()
        
        second_row_layout.addWidget(QLabel("Molecule:"))
        self._viz_molecule_input = QLineEdit()
        self._viz_molecule_input.setText("0")
        #self._viz_molecule_input.setMaximumWidth(80)
        self._viz_molecule_input.textChanged.connect(self._update_structure_spinbox_range)
        second_row_layout.addWidget(self._viz_molecule_input)
        
        second_row_layout.addWidget(QLabel("Structure:"))
        self._viz_structure_spinbox = QSpinBox()
        self._viz_structure_spinbox.setMinimum(0)
        self._viz_structure_spinbox.setValue(0)
        #self._viz_structure_spinbox.setMaximumWidth(80)
        second_row_layout.addWidget(self._viz_structure_spinbox)
        
        self._show_in_pymol_btn = QPushButton("Show in Pymol")
        self._show_in_pymol_btn.clicked.connect(self._show_in_pymol)
        second_row_layout.addWidget(self._show_in_pymol_btn)
        
        self._output_md_prepper_btn = QPushButton("Output MD Prepper Files")
        self._output_md_prepper_btn.clicked.connect(self._output_md_prepper_files)
        second_row_layout.addWidget(self._output_md_prepper_btn)
        
        viz_layout.addLayout(second_row_layout)
    
        # Third row: Get Redocking Range button and Centers/Lengths display
        third_row_layout = QHBoxLayout()
        
        self._get_redocking_range_btn = QPushButton("Get Redocking Range")
        self._get_redocking_range_btn.clicked.connect(self._get_redocking_range)
        third_row_layout.addWidget(self._get_redocking_range_btn)

        self._apply_redocking_btn = QPushButton("Apply to Docking")
        self._apply_redocking_btn.clicked.connect(self._apply_redocking_range)
        third_row_layout.addWidget(self._apply_redocking_btn)
        
        third_row_layout.addWidget(QLabel("Centers:"))
        self._redocking_centers_input = QLineEdit()
        self._redocking_centers_input.setReadOnly(True)
        self._redocking_centers_input.setPlaceholderText("x, y, z")
        third_row_layout.addWidget(self._redocking_centers_input)
        
        third_row_layout.addWidget(QLabel("Lengths:"))
        self._redocking_lengths_input = QLineEdit()
        self._redocking_lengths_input.setReadOnly(True)
        self._redocking_lengths_input.setPlaceholderText("x, y, z")
        third_row_layout.addWidget(self._redocking_lengths_input)
        
        viz_layout.addLayout(third_row_layout)
        
        main_layout.addWidget(viz_group)
        main_layout.addStretch()

    def _on_table_row_clicked(self, row, column):
        """Update molecule input when a table row is clicked."""
        # Get the item in the "Molecule Number" column (index 4)
        item = self._results_table.item(row, 4)
        if item:
            molecule_number = item.text()
            self._viz_molecule_input.setText(molecule_number)

    def update_results_table(self):
        """Update the results table with docking results."""
        if not self.main_window._docking_instance:
            self._results_table.setRowCount(0)
            return
        
        # Get combined results from docking instance
        combined_results = self.main_window._docking_instance.get_combined_results()
        
        # Populate table using shared helper
        self._populate_results_table_rows(combined_results)
        
        # Update structure spinbox maximum value for molecule 0
        if combined_results:
             num_structures = len(combined_results[0]['affinity'])
             self._viz_structure_spinbox.setMaximum(max(0, num_structures - 1))
    
    def _populate_results_table_rows(self, results: list[dict]) -> None:
        """Populate the results table with the given result data."""
        # Clear existing rows and set new row count
        self._results_table.setRowCount(0)
        self._results_table.setRowCount(len(results))
        
        # Populate table
        for idx, result in enumerate(results):
            # Number column
            number_item = QTableWidgetItem(str(result['number']))
            number_item.setTextAlignment(Qt.AlignCenter)
            self._results_table.setItem(idx, 0, number_item)
            
            # SMILES column
            smiles_item = QTableWidgetItem(result['smiles'])
            self._results_table.setItem(idx, 1, smiles_item)
            
            # Best Affinity column
            best_val = result['best_affinity']
            if isinstance(best_val, (int, float)):
                 best_item = QTableWidgetItem(f"{best_val:.2f}")
            else:
                 best_item = QTableWidgetItem(str(best_val))
            best_item.setTextAlignment(Qt.AlignCenter)
            self._results_table.setItem(idx, 2, best_item)
            
            # Affinity (list) column
            affinities = result['affinity']
            if affinities:
                affinity_str = ", ".join([f"{a:.2f}" for a in affinities])
                affinity_item = QTableWidgetItem(f"[{affinity_str}]")
            else:
                affinity_item = QTableWidgetItem("[]")
            self._results_table.setItem(idx, 3, affinity_item)
            
            # Molecule Number column (0-based)
            mol_num_item = QTableWidgetItem(str(result['molecule_number']))
            mol_num_item.setTextAlignment(Qt.AlignCenter)
            self._results_table.setItem(idx, 4, mol_num_item)

            # Filename column
            filename_item = QTableWidgetItem(result.get('filename', ''))
            self._results_table.setItem(idx, 5, filename_item)
        
        # Resize columns to contents
        self._results_table.resizeColumnsToContents()

    def _update_structure_spinbox_range(self):
        """Update the structure spinbox maximum value based on selected molecule."""
        if not self.main_window._docking_instance:
            return
        
        try:
            molecule_idx = int(self._viz_molecule_input.text())
        except ValueError:
            # Invalid input, reset structure to 0
            self._viz_structure_spinbox.setValue(0)
            return
        
        # Get affinity values
        affinity_values = self.main_window._docking_instance.get_affinity_values()
        
        if not affinity_values or molecule_idx < 0 or molecule_idx >= len(affinity_values):
            # Invalid molecule index, set structure max to 0
            self._viz_structure_spinbox.setMaximum(0)
            self._viz_structure_spinbox.setValue(0)
            return
        
        # Set structure max based on the number of poses for this molecule
        num_structures = len(affinity_values[molecule_idx])
        self._viz_structure_spinbox.setMaximum(max(0, num_structures - 1))
        
        # Reset structure to 0 if current value exceeds new maximum
        if self._viz_structure_spinbox.value() > num_structures - 1:
            self._viz_structure_spinbox.setValue(0)

    def _show_in_pymol(self):
        """Show selected docking result in PyMOL."""
        if not self.main_window._docking_instance:
            self.main_window._show_warning("Error", "No docking results available. Please run docking first.")
            return
        
        try:
            molecule_idx = int(self._viz_molecule_input.text())
            structure_idx = self._viz_structure_spinbox.value()
        except ValueError:
            self.main_window._show_warning("Error", "Please enter a valid molecule number.")
            return
        
        try:
            # Call open_results with the user-specified indices
            self.main_window._docking_instance.open_results(molecule_idx, structure_idx)
        except Exception as e:
            self.main_window._show_error("Error", f"Error opening PyMOL:\n{str(e)}")

    def _output_md_prepper_files(self):
        """Output MD prepper files for the selected molecule and structure."""
        if not self.main_window._docking_instance:
            self.main_window._show_warning("Error", "No docking results available. Please run docking first.")
            return
        
        try:
            molecule_idx = int(self._viz_molecule_input.text())
            structure_idx = self._viz_structure_spinbox.value()
        except ValueError:
            self.main_window._show_warning("Error", "Please enter a valid molecule number.")
            return
        
        try:
            # Get the output directory from the common settings
            output_dir = pathlib.Path(self.main_window._output_dir_input.text())
            if not output_dir or not output_dir.exists():
                self.main_window._show_warning("Error", "Please specify a valid output directory.")
                return
            
            # Create a subdirectory for MD prepper files
            md_prepper_dir = output_dir / "md_prepper"
            md_prepper_dir.mkdir(parents=True, exist_ok=True)
            
            # Call output_MD_prepper_files
            result_files = self.main_window._docking_instance.output_MD_prepper_files(
                ligand_index=molecule_idx,
                pose_index=structure_idx,
                output_dir=md_prepper_dir
            )
            
            # Build success message
            message_parts = [f"MD prepper files saved to:\n{md_prepper_dir}\n"]
            message_parts.append(f"\nGenerated files:")
            for key, path in result_files.items():
                message_parts.append(f"  - {path.name}")
            
            self.main_window._show_info("Success", "\n".join(message_parts))
            print(f"✓ Output MD prepper files to: {md_prepper_dir}")
            for key, path in result_files.items():
                print(f"  - {key}: {path}")
                
        except Exception as e:
            self.main_window._show_error("Error", f"Error outputting MD prepper files:\n{str(e)}")
            traceback.print_exc()

    def _get_redocking_range(self):
        """Get redocking range for selected molecule and structure."""
        if not self.main_window._docking_instance:
            self.main_window._show_warning("Error", "No docking results available. Please run docking first.")
            return
        
        try:
            molecule_idx = int(self._viz_molecule_input.text())
            structure_idx = self._viz_structure_spinbox.value()
        except ValueError:
            self.main_window._show_warning("Error", "Please enter a valid molecule number.")
            return
        
        try:
            # Call get_redocking_range from docking instance
            center, length = self.main_window._docking_instance.get_redocking_range(molecule_idx, structure_idx)
            
            # Format the results as comma-separated strings with 2 decimal places
            center_str = f"{center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}"
            length_str = f"{length[0]:.2f}, {length[1]:.2f}, {length[2]:.2f}"
            
            # Display in LineEdits
            self._redocking_centers_input.setText(center_str)
            self._redocking_lengths_input.setText(length_str)
            
            print(f"✓ Redocking range for molecule {molecule_idx}, structure {structure_idx}:")
            print(f"  Centers: [{center_str}]")
            print(f"  Lengths: [{length_str}]")
            
        except Exception as e:
            self.main_window._show_error("Error", f"Error getting redocking range:\n{str(e)}")

    def _apply_redocking_range(self):
        """Apply the calculated redocking range to the docking configuration."""
        center_text = self._redocking_centers_input.text()
        length_text = self._redocking_lengths_input.text()
        
        if not center_text or not length_text:
            self.main_window._show_warning("Error", "No redocking range calculated. Please click 'Get Redocking Range' first.")
            return

        try:
            # Parse values (simple string split)
            center_parts = [x.strip() for x in center_text.split(',')]
            length_parts = [x.strip() for x in length_text.split(',')]
            
            if len(center_parts) != 3 or len(length_parts) != 3:
                 raise ValueError("Invalid format in redocking fields")

            # Update Manual Input Fields in Protein Tab
            # Accessing protein_tab via main_window
            self.main_window._protein_tab._center_x_input.setText(center_parts[0])
            self.main_window._protein_tab._center_y_input.setText(center_parts[1])
            self.main_window._protein_tab._center_z_input.setText(center_parts[2])
            
            self.main_window._protein_tab._length_x_input.setText(length_parts[0])
            self.main_window._protein_tab._length_y_input.setText(length_parts[1])
            self.main_window._protein_tab._length_z_input.setText(length_parts[2])
            
            # Switch to "Specify Manually" mode
            self.main_window._protein_tab._radio_specify_manually.setChecked(True)
            
            self.main_window._show_info("Success", 
                "Redocking range applied to Docking configuration.\n"
                "Mode switched to 'Specify Manually' with updated coordinates.")

        except Exception as e:
            self.main_window._show_error("Error", f"Failed to apply redocking range: {e}")

    def _sort_results_table(self):
        """Sort results table by Best Affinity using backend logic."""
        if not self.main_window._docking_instance:
            return
        
        # Get sorted results directly from backend and populate table
        sorted_results = self.main_window._docking_instance.get_sorted_combined_results()
        self._populate_results_table_rows(sorted_results)

    def _cluster_results(self):
        """Cluster results and save to CSV."""
        if not self.main_window._docking_instance:
            self.main_window._show_warning("Error", "No docking results available. Please run docking first.")
            return

        # Get combined results from docking instance
        combined_results = self.main_window._docking_instance.get_combined_results()
        if not combined_results:
             self.main_window._show_warning("Error", "No results to cluster.")
             return

        # Ask for save path
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Clustered Results",
            self.main_window._last_dir,
            "CSV Files (*.csv)"
        )

        if not file_path:
            return

        try:
            # Get similarity threshold from input
            similarity_threshold = self._similarity_threshold_input.value()
            
            # Perform clustering
            self.main_window._docking_instance.cluster_results(similarity_threshold=similarity_threshold)
            
            # Save to CSV
            self.main_window._docking_instance.save_results_to_csv(file_path)
            
            self.main_window._show_info("Success", f"Results clustered and saved to:\n{file_path}")
            
        except Exception as e:
            self.main_window._show_error("Error", f"Error during clustering/saving:\n{str(e)}")
            print(f"Error: {e}")
