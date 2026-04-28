import pathlib
import traceback
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QRadioButton, QGroupBox, QDoubleSpinBox, QComboBox, QCheckBox
)
from PySide6.QtCore import Qt

import pdb_parser
from gui_utils import get_app_icon

class ProteinTab(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self._setup_ui()
        
    def _setup_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(15)
        
        # Create GroupBox for Process Protein
        process_group = QGroupBox("Process Protein")
        layout = QVBoxLayout(process_group)
        layout.setSpacing(15)
        
        # Radio buttons for input selection
        radio_layout = QHBoxLayout()
        radio_layout.addStretch()
        
        self._radio_pdb_id = QRadioButton("From PDB ID")
        self._radio_from_file = QRadioButton("From File")
        self._radio_pdb_id.setChecked(True)  # Default selection
        
        # Connect radio buttons to toggle function
        self._radio_pdb_id.toggled.connect(self._toggle_protein_input)
        self._radio_from_file.toggled.connect(self._toggle_protein_input)
        
        radio_layout.addWidget(self._radio_pdb_id)
        radio_layout.addWidget(self._radio_from_file)
        radio_layout.addStretch()
        
        layout.addLayout(radio_layout)
        
        # PDB ID input
        pdb_id_layout = QHBoxLayout()
        self._pdb_id_label = QLabel("PDB ID:")
        pdb_id_layout.addWidget(self._pdb_id_label)
        
        self._pdb_id_input = QLineEdit()
        self._pdb_id_input.setPlaceholderText("e.g., 1BBZ")
        pdb_id_layout.addWidget(self._pdb_id_input)
        
        layout.addLayout(pdb_id_layout)
        
        # Protein file input
        protein_file_layout = QHBoxLayout()
        self._protein_label = QLabel("Protein:")
        protein_file_layout.addWidget(self._protein_label)
        
        self._protein_file_input = QLineEdit()
        self._protein_file_input.setPlaceholderText("Select protein file...")
        protein_file_layout.addWidget(self._protein_file_input)
        
        self._browse_protein_btn = QPushButton("Browse")
        self._browse_protein_btn.clicked.connect(self._browse_protein_file)
        protein_file_layout.addWidget(self._browse_protein_btn)
        
        layout.addLayout(protein_file_layout)
        
        # pH and Chains input
        ph_chains_layout = QHBoxLayout()
        
        # pH
        ph_label = QLabel("pH:")
        ph_chains_layout.addWidget(ph_label)
        
        self._ph_spinbox = QDoubleSpinBox()
        self._ph_spinbox.setRange(0.0, 14.0)
        self._ph_spinbox.setValue(7.0)
        self._ph_spinbox.setSingleStep(0.1)
        self._ph_spinbox.setDecimals(1)
        ph_chains_layout.addWidget(self._ph_spinbox)
        
        # Chains
        chains_label = QLabel("Chains (Empty for All):")
        ph_chains_layout.addWidget(chains_label)
        
        self._chains_input = QLineEdit()
        self._chains_input.setText("A")
        ph_chains_layout.addWidget(self._chains_input)
        
        # Consider Ligand in Protonation checkbox
        self._consider_ligand_cb = QCheckBox("Consider Ligand in Protonation")
        self._consider_ligand_cb.setChecked(True)  # Default checked
        ph_chains_layout.addWidget(self._consider_ligand_cb)
        
        layout.addLayout(ph_chains_layout)

        
        # Preserve options (checkboxes for metal ions and waters) - Row 1
        preserve_options_row1 = QHBoxLayout()
        
        self._preserve_metal_cb = QCheckBox("Preserve Metal Ions")
        self._preserve_metal_cb.setChecked(True)  # Default checked
        preserve_options_row1.addWidget(self._preserve_metal_cb)
        
        self._preserve_metal_coord_waters_cb = QCheckBox("Preserve Metal-Coordinated Waters")
        self._preserve_metal_coord_waters_cb.toggled.connect(self._on_preserve_metal_coord_waters_toggled)
        preserve_options_row1.addWidget(self._preserve_metal_coord_waters_cb)
        
        preserve_options_row1.addStretch()
        layout.addLayout(preserve_options_row1)
        
        # Preserve options - Row 2
        preserve_options_row2 = QHBoxLayout()
        
        self._preserve_waters_2hbonds_cb = QCheckBox("Preserve Waters Forming 2 H-bonds with Proteins")
        self._preserve_waters_2hbonds_cb.toggled.connect(self._on_preserve_2hbonds_toggled)
        preserve_options_row2.addWidget(self._preserve_waters_2hbonds_cb)
        
        self._preserve_waters_3hbonds_cb = QCheckBox("Preserve Waters Forming 3 H-bonds with Proteins")
        preserve_options_row2.addWidget(self._preserve_waters_3hbonds_cb)
        
        preserve_options_row2.addStretch()
        layout.addLayout(preserve_options_row2)
        
        # Action buttons
        buttons_layout = QHBoxLayout()
        buttons_layout.addStretch()
        
        self._process_protein_btn = QPushButton("Process")
        self._process_protein_btn.clicked.connect(self.process_protein)
        buttons_layout.addWidget(self._process_protein_btn)
        
        self._view_pymol_btn = QPushButton("View in Pymol")
        self._view_pymol_btn.clicked.connect(self._view_in_pymol)
        buttons_layout.addWidget(self._view_pymol_btn)
        
        self._renew_btn = QPushButton("Renew")
        self._renew_btn.clicked.connect(self._renew_protein)
        buttons_layout.addWidget(self._renew_btn)
        
        buttons_layout.addStretch()
        
        layout.addLayout(buttons_layout)
        
        # Add GroupBox to main layout
        main_layout.addWidget(process_group)
        
        # Create GroupBox for Specify Docking Region
        docking_region_group = QGroupBox("Specify Docking Region")
        docking_layout = QVBoxLayout(docking_region_group)
        docking_layout.setSpacing(15)
        
        # First row: Radio buttons and Margin
        first_row_layout = QHBoxLayout()
        first_row_layout.addStretch()
        
        self._radio_blind_docking = QRadioButton("Blind Docking")
        self._radio_from_pocketeer = QRadioButton("From Pocketeer")
        self._radio_from_ligand = QRadioButton("From Ligand in PDB")
        self._radio_specify_manually = QRadioButton("Specify Manually")
        # Note: All radio buttons start unchecked
        
        # Connect radio buttons to toggle function
        self._radio_blind_docking.toggled.connect(self._toggle_docking_region)
        self._radio_from_pocketeer.toggled.connect(self._toggle_docking_region)
        self._radio_from_ligand.toggled.connect(self._toggle_docking_region)
        self._radio_specify_manually.toggled.connect(self._toggle_docking_region)
        
        # Connect radio buttons to apply docking range when selected
        self._radio_blind_docking.toggled.connect(self._apply_docking_range)
        self._radio_from_pocketeer.toggled.connect(self._on_pocketeer_selected)
        self._radio_from_ligand.toggled.connect(self._apply_docking_range)
        self._radio_specify_manually.toggled.connect(self._apply_docking_range)
        
        first_row_layout.addWidget(self._radio_blind_docking)
        first_row_layout.addWidget(self._radio_from_pocketeer)
        first_row_layout.addWidget(self._radio_from_ligand)
        first_row_layout.addWidget(self._radio_specify_manually)
        
        # Margin
        first_row_layout.addWidget(QLabel("Margin:"))
        self._margin_spinbox = QDoubleSpinBox()
        self._margin_spinbox.setRange(0.0, 100.0)
        self._margin_spinbox.setValue(6.0)
        self._margin_spinbox.setSingleStep(0.5)
        self._margin_spinbox.setDecimals(1)
        self._margin_spinbox.setMaximumWidth(80)
        # Connect to update docking region when margin changes
        self._margin_spinbox.valueChanged.connect(self._on_margin_changed)
        first_row_layout.addWidget(self._margin_spinbox)
        
        #first_row_layout.addStretch()
        docking_layout.addLayout(first_row_layout)
        
        # Second row: Select Pocket (for Pocketeer)
        pocket_row_layout = QHBoxLayout()
        pocket_row_layout.addWidget(QLabel("Select Pocket:"))
        self._pocket_combo = QComboBox()
        self._pocket_combo.addItem("No pockets detected")
        # Connect to update docking region when pocket selection changes
        self._pocket_combo.currentIndexChanged.connect(self._on_pocket_selection_changed)
        pocket_row_layout.addWidget(self._pocket_combo)
        docking_layout.addLayout(pocket_row_layout)
        
        # Third row: Select Ligand
        ligand_row_layout = QHBoxLayout()
        ligand_row_layout.addWidget(QLabel("Select Ligand:"))
        self._ligand_combo = QComboBox()
        self._ligand_combo.addItem("No ligands available")
        # Connect to update docking region when ligand selection changes
        self._ligand_combo.currentTextChanged.connect(self._on_ligand_selection_changed)
        ligand_row_layout.addWidget(self._ligand_combo)
        docking_layout.addLayout(ligand_row_layout)
        
        # Centers row
        centers_row_layout = QHBoxLayout()
        centers_row_layout.addWidget(QLabel("Centers:"))
        
        centers_row_layout.addWidget(QLabel("X:"))
        self._center_x_input = QLineEdit()
        self._center_x_input.setText("0")
        # Connect to update docking region when center changes
        self._center_x_input.editingFinished.connect(self._on_manual_range_changed)
        centers_row_layout.addWidget(self._center_x_input)
        
        centers_row_layout.addWidget(QLabel("Y:"))
        self._center_y_input = QLineEdit()
        self._center_y_input.setText("0")
        self._center_y_input.editingFinished.connect(self._on_manual_range_changed)
        centers_row_layout.addWidget(self._center_y_input)
        
        centers_row_layout.addWidget(QLabel("Z:"))
        self._center_z_input = QLineEdit()
        self._center_z_input.setText("0")
        self._center_z_input.editingFinished.connect(self._on_manual_range_changed)
        centers_row_layout.addWidget(self._center_z_input)
        docking_layout.addLayout(centers_row_layout)
        
        # Lengths row
        lengths_row_layout = QHBoxLayout()
        lengths_row_layout.addWidget(QLabel("Lengths:"))
        
        lengths_row_layout.addWidget(QLabel("X:"))
        self._length_x_input = QLineEdit()
        self._length_x_input.setText("20")
        # Connect to update docking region when length changes
        self._length_x_input.editingFinished.connect(self._on_manual_range_changed)
        lengths_row_layout.addWidget(self._length_x_input)
        
        lengths_row_layout.addWidget(QLabel("Y:"))
        self._length_y_input = QLineEdit()
        self._length_y_input.setText("20")
        self._length_y_input.editingFinished.connect(self._on_manual_range_changed)
        lengths_row_layout.addWidget(self._length_y_input)
        
        lengths_row_layout.addWidget(QLabel("Z:"))
        self._length_z_input = QLineEdit()
        self._length_z_input.setText("20")
        self._length_z_input.editingFinished.connect(self._on_manual_range_changed)
        lengths_row_layout.addWidget(self._length_z_input)
        docking_layout.addLayout(lengths_row_layout)
        
        # Detect Flexible Residues row
        detect_flex_layout = QHBoxLayout()
        detect_flex_layout.addWidget(QLabel("Detect Flexible Residues"))
        
        self._detect_flex_btn = QPushButton("Detect")
        self._detect_flex_btn.clicked.connect(self._detect_flexible_residues)
        detect_flex_layout.addWidget(self._detect_flex_btn)
        
        detect_flex_layout.addWidget(QLabel("Result:"))
        
        self._flex_result_input = QLineEdit()
        self._flex_result_input.setReadOnly(True)
        self._flex_result_input.setPlaceholderText("Detected residues will appear here...")
        detect_flex_layout.addWidget(self._flex_result_input)
        
        docking_layout.addLayout(detect_flex_layout)
        
        # View in Pymol and Add Chain Identifiers buttons for docking region
        docking_btn_layout = QHBoxLayout()
        docking_btn_layout.addStretch()
        
        self._add_chain_ids_btn = QPushButton("Re-add Chain Identifiers")
        self._add_chain_ids_btn.clicked.connect(self._add_chain_identifiers)
        docking_btn_layout.addWidget(self._add_chain_ids_btn)
        
        self._view_docking_region_btn = QPushButton("View in Pymol")
        self._view_docking_region_btn.clicked.connect(self._view_docking_region_in_pymol)
        docking_btn_layout.addWidget(self._view_docking_region_btn)
        
        docking_btn_layout.addStretch()
        docking_layout.addLayout(docking_btn_layout)
        
        # Add docking region GroupBox to main layout
        main_layout.addWidget(docking_region_group)
        main_layout.addStretch()
        
        # Initialize the enabled/disabled state
        self._toggle_protein_input()
        self._toggle_docking_region()


    def _toggle_protein_input(self):
        """Toggle enabled/disabled state based on radio button selection."""
        if self._radio_pdb_id.isChecked():
            # Enable PDB ID input
            self._pdb_id_label.setEnabled(True)
            self._pdb_id_input.setEnabled(True)
            
            # Disable file input
            self._protein_label.setEnabled(False)
            self._protein_file_input.setEnabled(False)
            self._browse_protein_btn.setEnabled(False)
        else:
            # Disable PDB ID input
            self._pdb_id_label.setEnabled(False)
            self._pdb_id_input.setEnabled(False)
            
            # Enable file input
            self._protein_label.setEnabled(True)
            self._protein_file_input.setEnabled(True)
            self._browse_protein_btn.setEnabled(True)

    def _browse_protein_file(self):
        """Browse for protein file."""
        from PySide6.QtWidgets import QFileDialog
        
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Protein File",
            self.main_window._last_dir,
            "PDB Files (*.pdb);;All Files (*.*)"
        )
        
        if file_path:
            self.main_window._protein_file = file_path
            self._protein_file_input.setText(file_path)
            self.main_window._last_dir = str(pathlib.Path(file_path).parent)

    def _on_preserve_2hbonds_toggled(self, checked):
        """Handle the 2 H-bonds checkbox toggle.
        
        When 2 H-bonds is checked, the 3 H-bonds option should be 
        automatically checked and disabled (grayed out), because 
        waters forming 2 H-bonds automatically includes waters forming 3+ H-bonds.
        """
        if checked:
            # When 2 H-bonds is checked, auto-check and disable 3 H-bonds
            self._preserve_waters_3hbonds_cb.setChecked(True)
            self._preserve_waters_3hbonds_cb.setEnabled(False)
        else:
            # When 2 H-bonds is unchecked, re-enable 3 H-bonds
            self._preserve_waters_3hbonds_cb.setEnabled(True)

    def _on_preserve_metal_coord_waters_toggled(self, checked):
        """Handle the Metal-Coordinated Waters checkbox toggle.
        
        When Metal-Coordinated Waters is checked, the Metal Ions option should be 
        automatically checked and disabled (grayed out), because 
        metal-coordinated waters require metal ions to be preserved first.
        """
        if checked:
            # When Metal-Coordinated Waters is checked, auto-check and disable Metal Ions
            self._preserve_metal_cb.setChecked(True)
            self._preserve_metal_cb.setEnabled(False)
        else:
            # When Metal-Coordinated Waters is unchecked, re-enable Metal Ions
            self._preserve_metal_cb.setEnabled(True)

    def process_protein(self):
        """Process the protein with selected parameters."""
        try:
            # Validate output directory
            if not self.main_window._output_directory:
                self.main_window._show_warning("Error", "Please select an output directory first.")
                return
            
            # Get parameters
            ph = self._ph_spinbox.value()
            chains = self._chains_input.text().strip()
            
            # Prepare paths
            output_dir = pathlib.Path(self.main_window._output_directory)
            receptor_dir = output_dir / "receptor"
            receptor_dir.mkdir(parents=True, exist_ok=True)
            output_pdb = receptor_dir / "receptor.pdb"
            
            # Create PDBParser based on radiobutton selection
            if self._radio_pdb_id.isChecked():
                pdb_id = self._pdb_id_input.text().strip()
                if not pdb_id:
                    self.main_window._show_warning("Error", "Please enter a PDB ID.")
                    return
                
                # For PDB ID, we need a temporary input file path
                input_pdb = receptor_dir / f"{pdb_id}.pdb"
                self.main_window._protein_parser = pdb_parser.PDBParser(
                    input_pdb=str(input_pdb),
                    output_pdb=str(output_pdb),
                    pdb_id=pdb_id,
                    unit_structure=False
                )
            else:
                if not self.main_window._protein_file:
                    self.main_window._show_warning("Error", "Please select a protein file.")
                    return
                
                self.main_window._protein_parser = pdb_parser.PDBParser(
                    input_pdb=self.main_window._protein_file,
                    output_pdb=str(output_pdb),
                    pdb_id=None,
                    unit_structure=False
                )
            
            # Get preserve options from checkboxes
            preserve_metal = self._preserve_metal_cb.isChecked()
            preserve_metal_coord_water = self._preserve_metal_coord_waters_cb.isChecked()
            
            # Determine preserve_coord_water parameter:
            # - If 2 H-bonds is checked, use 2 (which also includes 3 H-bonds waters)
            # - If only 3 H-bonds is checked, use 3
            # - Otherwise, use False
            if self._preserve_waters_2hbonds_cb.isChecked():
                preserve_coord_water = 2
            elif self._preserve_waters_3hbonds_cb.isChecked():
                preserve_coord_water = 3
            else:
                preserve_coord_water = False
            
            # Determine ligand_mol2 path based on checkbox
            # If checked, save original ligands to receptor/original_ligand.mol2
            if self._consider_ligand_cb.isChecked():
                ligand_mol2 = receptor_dir / "original_ligand.mol2"
            else:
                ligand_mol2 = None
            
            # Execute processing steps in order. If ligand-aware protonation
            # breaks pdb2pqr, warn once and retry without the ligand file.
            try:
                self.main_window._protein_parser.process_protein_pipeline(
                    ph=ph,
                    chains=chains,
                    preserve_metal=preserve_metal,
                    preserve_coord_water=preserve_coord_water,
                    preserve_metal_coord_water=preserve_metal_coord_water,
                    ligand_mol2=ligand_mol2
                )
            except pdb_parser.third_party_tools.PDB2PQRError as e:
                if not ligand_mol2:
                    raise

                self.main_window._show_warning(
                    "Warning",
                    "PDB2PQR failed while considering the ligand in protonation.\n\n"
                    f"Error details:\n{e}\n\n"
                    "After you click OK, processing will continue without "
                    "considering the ligand in protonation."
                )
                print(
                    "Warning: PDB2PQR failed with ligand-aware protonation; "
                    "retrying without ligand."
                )
                self._consider_ligand_cb.setChecked(False)
                ligand_mol2 = None
                self.main_window._protein_parser.process_protein_pipeline(
                    ph=ph,
                    chains=chains,
                    preserve_metal=preserve_metal,
                    preserve_coord_water=preserve_coord_water,
                    preserve_metal_coord_water=preserve_metal_coord_water,
                    ligand_mol2=None
                )

            
            # Check if protein has chain identifiers, if not, add them automatically
            if not self.main_window._protein_parser.has_chain_identifiers():
                self.main_window._show_warning("Warning", 
                    "The PDB file does not have Chain IDs.\n\n"
                    "Chain identifiers will be added automatically to ensure "
                    "compatibility with flexible residue detection and other features.")
                self.main_window._protein_parser.add_chain_identifiers()
                print("✓ Chain identifiers added automatically (PDB was missing chain IDs)")
            
            # Update ligand combo with available HETATM labels
            # Block signals to prevent triggering _on_ligand_selection_changed during update
            self._ligand_combo.blockSignals(True)
            self._ligand_combo.clear()
            hetatm_labels = self.main_window._protein_parser.get_hetatm_labels()
            if hetatm_labels:
                self._ligand_combo.addItems(hetatm_labels)
            else:
                self._ligand_combo.addItem("No ligands available")
            self._ligand_combo.blockSignals(False)
            
            self.main_window._show_info("Success", 
                f"Protein processing completed successfully!\n\n"
                f"Output PDB: {output_pdb}\n"
                f"Output PDBQT: {self.main_window._protein_parser.get_generated_pdbqt_path()}")
            
        except Exception as e:
            self.main_window._show_error("Error", f"Error processing protein:\n{str(e)}")
            print(f"Error: {e}")
            traceback.print_exc()

    def _view_in_pymol(self):
        """Open the processed protein in PyMOL."""
        if not self.main_window._protein_parser:
            self.main_window._show_warning("Error", "Please process the protein first.")
            return
        
        try:
            self.main_window._protein_parser.open_output_pdb_only()
        except Exception as e:
            self.main_window._show_error("Error", f"Error opening PyMOL:\n{str(e)}")

    def _renew_protein(self):
        """Renew/reset the protein processing."""
        if not self.main_window._protein_parser:
            self.main_window._show_warning("Error", "Please process the protein first.")
            return
        
        try:
            self.main_window._protein_parser.renew_output()
            print("✓ Output file renewed - reset to original input PDB")
            self.main_window._show_info("Success", "Output file has been renewed from the original input PDB.")
        except Exception as e:
            self.main_window._show_error("Error", f"Error renewing output:\n{str(e)}")

    def _toggle_docking_region(self):
        """Toggle enabled/disabled state for docking region based on radio button selection."""
        if self._radio_blind_docking.isChecked():
            # Blind Docking: only margin is enabled
            self._pocket_combo.setEnabled(False)
            self._ligand_combo.setEnabled(False)
            self._center_x_input.setEnabled(False)
            self._center_y_input.setEnabled(False)
            self._center_z_input.setEnabled(False)
            self._length_x_input.setEnabled(False)
            self._length_y_input.setEnabled(False)
            self._length_z_input.setEnabled(False)
            self._margin_spinbox.setEnabled(True)
        elif self._radio_from_pocketeer.isChecked():
            # From Pocketeer: pocket combo and margin are enabled
            self._pocket_combo.setEnabled(True)
            self._ligand_combo.setEnabled(False)
            self._center_x_input.setEnabled(False)
            self._center_y_input.setEnabled(False)
            self._center_z_input.setEnabled(False)
            self._length_x_input.setEnabled(False)
            self._length_y_input.setEnabled(False)
            self._length_z_input.setEnabled(False)
            self._margin_spinbox.setEnabled(True)
        elif self._radio_from_ligand.isChecked():
            # From Ligand: ligand combo and margin are enabled
            self._pocket_combo.setEnabled(False)
            self._ligand_combo.setEnabled(True)
            self._center_x_input.setEnabled(False)
            self._center_y_input.setEnabled(False)
            self._center_z_input.setEnabled(False)
            self._length_x_input.setEnabled(False)
            self._length_y_input.setEnabled(False)
            self._length_z_input.setEnabled(False)
            self._margin_spinbox.setEnabled(True)
        elif self._radio_specify_manually.isChecked():
            # Specify Manually: centers and lengths are enabled, margin disabled
            self._pocket_combo.setEnabled(False)
            self._ligand_combo.setEnabled(False)
            self._center_x_input.setEnabled(True)
            self._center_y_input.setEnabled(True)
            self._center_z_input.setEnabled(True)
            self._length_x_input.setEnabled(True)
            self._length_y_input.setEnabled(True)
            self._length_z_input.setEnabled(True)
            self._margin_spinbox.setEnabled(False)
        else:
            # No radio button selected: disable all
            self._pocket_combo.setEnabled(False)
            self._ligand_combo.setEnabled(False)
            self._center_x_input.setEnabled(False)
            self._center_y_input.setEnabled(False)
            self._center_z_input.setEnabled(False)
            self._length_x_input.setEnabled(False)
            self._length_y_input.setEnabled(False)
            self._length_z_input.setEnabled(False)
            self._margin_spinbox.setEnabled(False)

    def _apply_docking_range(self):
        """Apply the docking range based on the selected radio button."""
        # Only apply if a radio button is toggled to True (not False)
        sender = self.sender()
        if not sender.isChecked():
            return
        
        if not self.main_window._protein_parser:
            self.main_window._show_warning("Error", "Please process the protein first.")
            # Temporarily disable auto-exclusive to allow unchecking all radio buttons
            sender.setAutoExclusive(False)
            sender.blockSignals(True)
            sender.setChecked(False)
            sender.blockSignals(False)
            sender.setAutoExclusive(True)
            # Reset center and length fields to default values
            self._center_x_input.setText("0")
            self._center_y_input.setText("0")
            self._center_z_input.setText("0")
            self._length_x_input.setText("20")
            self._length_y_input.setText("20")
            self._length_z_input.setText("20")
            # Manually update UI state since we blocked signals
            self._toggle_docking_region()
            return
        
        try:
            if self._radio_blind_docking.isChecked():
                # Blind Docking
                margin = self._margin_spinbox.value()
                self.main_window._protein_parser.set_blind_docking_range(margin=margin)
                print(f"✓ Set blind docking range with margin: {margin}")
                self.main_window._show_info("Success", 
                    f"Blind docking range set with margin: {margin} Å")
            

            elif self._radio_from_ligand.isChecked():
                # From Ligand in PDB
                ligand_label = self._ligand_combo.currentText()
                if ligand_label == "No ligands available":
                    self.main_window._show_warning("Error", "No ligands available in the PDB file.")
                    sender.setChecked(False)
                    return
                
                margin = self._margin_spinbox.value()
                self.main_window._protein_parser.set_known_ligand_range(ligand_label, margin=margin)
                print(f"✓ Set docking range from ligand: {ligand_label}, margin: {margin}")
                self.main_window._show_info("Success", 
                    f"Docking range set from ligand: {ligand_label}\nMargin: {margin} Å")
            
            elif self._radio_specify_manually.isChecked():
                # Specify Manually
                try:
                    center_x = float(self._center_x_input.text())
                    center_y = float(self._center_y_input.text())
                    center_z = float(self._center_z_input.text())
                    length_x = float(self._length_x_input.text())
                    length_y = float(self._length_y_input.text())
                    length_z = float(self._length_z_input.text())
                except ValueError:
                    self.main_window._show_warning("Error", "Please enter valid numeric values for centers and lengths.")
                    sender.setChecked(False)
                    return
                
                # Validate that length values are positive
                if length_x <= 0 or length_y <= 0 or length_z <= 0:
                    self.main_window._show_warning("Error", "Length values must be positive (greater than 0).")
                    sender.setChecked(False)
                    return
                
                center = [center_x, center_y, center_z]
                side_length = [length_x, length_y, length_z]
                self.main_window._protein_parser.set_custom_docking_range(center=center, side_length=side_length)
                print(f"✓ Set custom docking range - Center: {center}, Size: {side_length}")
                self.main_window._show_info("Success", 
                    f"Custom docking range set\n"
                    f"Center: ({center_x}, {center_y}, {center_z})\n"
                    f"Size: ({length_x}, {length_y}, {length_z})")
        
        except Exception as e:
            self.main_window._show_error("Error", f"Error setting docking range:\n{str(e)}")
            print(f"Error: {e}")
            traceback.print_exc()
            sender.setChecked(False)

    def _on_pocketeer_selected(self, checked):
        """Handle Pocketeer radio button selection - detect pockets."""
        if not checked:
            return
        
        if not self.main_window._protein_parser:
            self.main_window._show_warning("Error", "Please process the protein first.")
            # Temporarily disable auto-exclusive to allow unchecking all radio buttons
            self._radio_from_pocketeer.setAutoExclusive(False)
            self._radio_from_pocketeer.blockSignals(True)
            self._radio_from_pocketeer.setChecked(False)
            self._radio_from_pocketeer.blockSignals(False)
            self._radio_from_pocketeer.setAutoExclusive(True)
            # Manually update UI state since we blocked signals
            self._toggle_docking_region()
            return
        
        try:
            # Detect pockets using pocketeer
            print("Detecting pockets with pocketeer...")
            pockets = self.main_window._protein_parser.detect_pockets()
            self.main_window._detected_pockets = self.main_window._protein_parser.get_detected_pockets() or []
            
            # Update pocket combo box
            # Block signals to prevent triggering _on_pocket_selection_changed during update
            self._pocket_combo.blockSignals(True)
            self._pocket_combo.clear()
            if self.main_window._detected_pockets:
                for i, pocket in enumerate(self.main_window._detected_pockets):
                    center, side_length = pocket
                    # Format: "Pocket N: Center (x, y, z), Size (lx, ly, lz)"
                    label = (f"Pocket {i+1}: Center ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}), "
                            f"Size ({side_length[0]:.1f}, {side_length[1]:.1f}, {side_length[2]:.1f})")
                    self._pocket_combo.addItem(label)
                self._pocket_combo.blockSignals(False)
                
                print(f"✓ Detected {len(self.main_window._detected_pockets)} pocket(s)")
                self.main_window._show_info("Success", f"Detected {len(self.main_window._detected_pockets)} pocket(s).")
                
                # Auto-select first pocket and update docking region
                if len(self.main_window._detected_pockets) > 0:
                    self._apply_pocket_as_docking_region(0)
            else:
                self._pocket_combo.addItem("No pockets detected")
                self._pocket_combo.blockSignals(False)
                self.main_window._show_warning("Warning", "No pockets detected in the protein structure.")
                
        except Exception as e:
            self.main_window._show_error("Error", f"Error detecting pockets:\n{str(e)}")
            print(f"Error: {e}")
            traceback.print_exc()
            self._radio_from_pocketeer.setChecked(False)

    def _on_pocket_selection_changed(self, index):
        """Handle pocket selection changes in the combo box."""
        # Only update if "From Pocketeer" mode is active
        if not self._radio_from_pocketeer.isChecked():
            return
        
        # Only update if protein parser exists and we have detected pockets
        if not self.main_window._protein_parser:
            return
        
        if index < 0 or not self.main_window._detected_pockets or index >= len(self.main_window._detected_pockets):
            return
        
        self._apply_pocket_as_docking_region(index)

    def _apply_pocket_as_docking_region(self, pocket_index):
        """Apply the selected pocket as the docking region."""
        if pocket_index < 0 or pocket_index >= len(self.main_window._detected_pockets):
            return
        
        try:
            import numpy as np
            
            pocket = self.main_window._detected_pockets[pocket_index]
            center, side_length = pocket
            
            # Ensure they are numpy arrays
            center = np.asarray(center, dtype=float)
            side_length = np.asarray(side_length, dtype=float)
            
            # Get current margin
            margin = self._margin_spinbox.value()
            
            # Apply margin to side length
            margin_arr = np.full(3, margin, dtype=float)
            expanded_side_length = side_length + margin_arr
            
            # Set the docking region with margin
            self.main_window._protein_parser.set_custom_docking_range(center=center, side_length=expanded_side_length)
            
            # Verify it was set
            docking_region = self.main_window._protein_parser.get_generated_docking_range()
            print(f"✓ Set docking region from Pocket {pocket_index + 1}")
            print(f"  Center: {center}")
            print(f"  Side length (with margin {margin}): {expanded_side_length}")
            if docking_region:
                print(f"  Verified: docking region is set")
            else:
                print(f"  WARNING: docking region is None after setting!")
            
        except Exception as e:
            print(f"Warning: Could not apply pocket as docking region: {e}")
            traceback.print_exc()

    def _on_ligand_selection_changed(self, ligand_label):
        """Handle ligand selection changes in the combo box."""
        # Only update if "From Ligand in PDB" mode is active
        if not self._radio_from_ligand.isChecked():
            return
        
        # Only update if protein parser exists and ligand is valid
        if not self.main_window._protein_parser:
            return
        
        if ligand_label == "No ligands available" or not ligand_label:
            return
        
        try:
            # Get current margin value
            margin = self._margin_spinbox.value()
            
            # Update docking range with the newly selected ligand
            self.main_window._protein_parser.set_known_ligand_range(ligand_label, margin=margin)
            print(f"✓ Updated docking range for ligand: {ligand_label}, margin: {margin}")
            
        except Exception as e:
            print(f"Warning: Could not update docking range for ligand {ligand_label}: {e}")

    def _on_margin_changed(self, margin):
        """Handle margin value changes."""
        # Only update if protein parser exists
        if not self.main_window._protein_parser:
            return
        
        try:
            if self._radio_blind_docking.isChecked():
                # Update blind docking range with new margin
                self.main_window._protein_parser.set_blind_docking_range(margin=margin)
                print(f"✓ Updated blind docking range with margin: {margin}")
            
            elif self._radio_from_pocketeer.isChecked():
                # Update pocket-based docking range with new margin
                pocket_index = self._pocket_combo.currentIndex()
                if pocket_index >= 0 and self.main_window._detected_pockets and pocket_index < len(self.main_window._detected_pockets):
                    self._apply_pocket_as_docking_region(pocket_index)
                    print(f"✓ Updated pocket docking range with margin: {margin}")
                
            elif self._radio_from_ligand.isChecked():
                # Update ligand-based docking range with new margin
                ligand_label = self._ligand_combo.currentText()
                if ligand_label != "No ligands available" and ligand_label:
                    self.main_window._protein_parser.set_known_ligand_range(ligand_label, margin=margin)
                    print(f"✓ Updated docking range for ligand: {ligand_label}, margin: {margin}")
                    
        except Exception as e:
            print(f"Warning: Could not update docking range with margin {margin}: {e}")

    def _on_manual_range_changed(self):
        """Handle manual center/length value changes."""
        # Only update if "Specify Manually" mode is active
        if not self._radio_specify_manually.isChecked():
            return
        
        # Only update if protein parser exists
        if not self.main_window._protein_parser:
            return
        
        try:
            # Parse center and length values
            center_x = float(self._center_x_input.text())
            center_y = float(self._center_y_input.text())
            center_z = float(self._center_z_input.text())
            length_x = float(self._length_x_input.text())
            length_y = float(self._length_y_input.text())
            length_z = float(self._length_z_input.text())
            
            # Validate that length values are positive
            if length_x <= 0 or length_y <= 0 or length_z <= 0:
                print("Warning: Length values must be positive (greater than 0)")
                return
            
            center = [center_x, center_y, center_z]
            side_length = [length_x, length_y, length_z]
            
            # Update docking range
            self.main_window._protein_parser.set_custom_docking_range(center=center, side_length=side_length)
            print(f"✓ Updated custom docking range - Center: {center}, Size: {side_length}")
            
        except ValueError:
            # Invalid input (non-numeric), silently ignore (user is still typing)
            pass
        except Exception as e:
            print(f"Warning: Could not update manual docking range: {e}")

    def _detect_flexible_residues(self):
        """Detect flexible residues within docking region + margin."""
        if not self.main_window._protein_parser:
            self.main_window._show_warning("Error", "Please process the protein first.")
            return
        
        # Check if a docking range has been set
        if not any([self._radio_blind_docking.isChecked(), 
                    self._radio_from_pocketeer.isChecked(),
                    self._radio_from_ligand.isChecked(), 
                    self._radio_specify_manually.isChecked()]):
            self.main_window._show_warning("Error", "Please select and apply a docking range mode first.")
            return
        
        try:
            # Call detect_flexible_residues from pdb_parser (search within selected range, margin=0.0)
            residues = self.main_window._protein_parser.detect_flexible_residues(margin=0.0)
            
            # Convert list to comma-separated string
            result_str = ",".join(residues)
            
            # Display in result input field
            self._flex_result_input.setText(result_str)
            
            print(f"✓ Detected {len(residues)} flexible residues within selected docking region")
            self.main_window._show_info("Success", 
                f"Detected {len(residues)} flexible residues:\n{result_str}")
            
        except ValueError as e:
            # Handle the case when protein doesn't have chain IDs
            self.main_window._show_error("Error", 
                f"Cannot detect flexible residues:\n{str(e)}\n\n"
                "Try clicking 'Add Chain Identifiers' first.")
            print(f"Error: {e}")
        except Exception as e:
            self.main_window._show_error("Error", f"Error detecting flexible residues:\n{str(e)}")
            print(f"Error: {e}")
            traceback.print_exc()

    def _add_chain_identifiers(self):
        """Add chain identifiers to the protein PDB file."""
        if not self.main_window._protein_parser:
            self.main_window._show_warning("Error", "Please process the protein first.")
            return
        
        try:
            # Call add_chain_identifiers from pdb_parser
            self.main_window._protein_parser.add_chain_identifiers()
            
            print("✓ Chain identifiers added to protein")
            self.main_window._show_info("Success", 
                "Chain identifiers have been added to the protein structure.")
            
        except Exception as e:
            self.main_window._show_error("Error", f"Error adding chain identifiers:\n{str(e)}")
            print(f"Error: {e}")
            traceback.print_exc()

    def _view_docking_region_in_pymol(self):
        """View the docking region in PyMOL."""
        if not self.main_window._protein_parser:
            self.main_window._show_warning("Error", "Please process the protein first.")
            return
        
        # Check if a docking range has been set
        if not any([self._radio_blind_docking.isChecked(), 
                    self._radio_from_pocketeer.isChecked(),
                    self._radio_from_ligand.isChecked(), 
                    self._radio_specify_manually.isChecked()]):
            self.main_window._show_warning("Error", "Please select and apply a docking range mode first.")
            return
        
        try:
            self.main_window._protein_parser.open_output_pdb_with_box(line_width=5.0)
        except Exception as e:
            self.main_window._show_error("Error", f"Error opening PyMOL with docking box:\n{str(e)}")
