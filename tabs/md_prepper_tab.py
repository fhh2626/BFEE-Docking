import pathlib
import traceback
import shutil
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QGroupBox, QDoubleSpinBox, QComboBox, QApplication, QFileDialog, QInputDialog
)

import md_prepper

class MDPrepperTab(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self._setup_ui()
        
    def _setup_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(15)
        
        # PDB Preparation GroupBox
        pdb_prep_group = QGroupBox("PDB Preparation")
        pdb_prep_layout = QVBoxLayout(pdb_prep_group)
        
        # VMD Path Row
        vmd_path_layout = QHBoxLayout()
        vmd_path_layout.addWidget(QLabel("VMD Path:"))
        
        self._vmd_path_input = QLineEdit()
        self._vmd_path_input.setPlaceholderText("Select VMD executable...")
        if self.main_window._vmd_path:
            self._vmd_path_input.setText(self.main_window._vmd_path)
            
        # Connect text change to save settings (optional, but good for manual entry)
        self._vmd_path_input.editingFinished.connect(self._save_vmd_path_from_input)
            
        vmd_path_layout.addWidget(self._vmd_path_input)
        
        browse_vmd_btn = QPushButton("Browse")
        browse_vmd_btn.clicked.connect(self._browse_vmd_path)
        vmd_path_layout.addWidget(browse_vmd_btn)
        
        pdb_prep_layout.addLayout(vmd_path_layout)
        
        # Second row: pH and Process buttons
        process_layout = QHBoxLayout()
        
        process_layout.addWidget(QLabel("pH:"))
        
        self._md_ph_spinbox = QDoubleSpinBox()
        self._md_ph_spinbox.setRange(0.0, 14.0)
        self._md_ph_spinbox.setValue(7.0)
        self._md_ph_spinbox.setSingleStep(0.1)
        self._md_ph_spinbox.setDecimals(1)
        process_layout.addWidget(self._md_ph_spinbox)
        
        process_layout.addStretch()
        
        self._md_process_protein_btn = QPushButton("Process Protein")
        self._md_process_protein_btn.clicked.connect(self._md_process_protein)
        process_layout.addWidget(self._md_process_protein_btn)
        
        self._md_process_ligand_btn = QPushButton("Process Ligand")
        self._md_process_ligand_btn.clicked.connect(self._md_process_ligand)
        process_layout.addWidget(self._md_process_ligand_btn)
        
        pdb_prep_layout.addLayout(process_layout)
        main_layout.addWidget(pdb_prep_group)

        # CHARMM PSF Generation GroupBox
        charmm_group = QGroupBox("CHARMM PSF Generation")
        charmm_layout = QVBoxLayout(charmm_group)
        charmm_layout.setSpacing(10)

        # Row 1: Ligand str File
        str_file_layout = QHBoxLayout()
        str_file_layout.addWidget(QLabel("Ligand str File:"))
        
        self._ligand_str_edit = QLineEdit()
        str_file_layout.addWidget(self._ligand_str_edit)
        
        browse_str_btn = QPushButton("Browse")
        browse_str_btn.clicked.connect(self._browse_ligand_str_file)
        str_file_layout.addWidget(browse_str_btn)
        
        charmm_layout.addLayout(str_file_layout)

        # Row 2: Water Edge and Ions
        row2_layout = QHBoxLayout()
        
        row2_layout.addWidget(QLabel("Water Padding:"))
        self._water_edge_spin = QDoubleSpinBox()
        self._water_edge_spin.setValue(15.0)
        self._water_edge_spin.setSingleStep(1.0)
        row2_layout.addWidget(self._water_edge_spin)
        
        row2_layout.addWidget(QLabel("Ion Type:"))
        self._ion_type_combo = QComboBox()
        self._ion_type_combo.addItems(["NaCl", "KCl", "CaCl2", "MgCl2"])
        row2_layout.addWidget(self._ion_type_combo)

        row2_layout.addWidget(QLabel("Ion Concentration:"))
        self._ion_concentration_spin = QDoubleSpinBox()
        self._ion_concentration_spin.setValue(0.1)
        self._ion_concentration_spin.setSingleStep(0.01)
        row2_layout.addWidget(self._ion_concentration_spin)

        row2_layout.addStretch()
        charmm_layout.addLayout(row2_layout)

        # Row 3: Generate PSF Button
        psf_btn_layout = QHBoxLayout()
        psf_btn_layout.addStretch()
        
        self._generate_psf_btn = QPushButton("Generate PSF")
        self._generate_psf_btn.clicked.connect(self._generate_charmm_psf)
        psf_btn_layout.addWidget(self._generate_psf_btn)
        
        psf_btn_layout.addStretch()
        charmm_layout.addLayout(psf_btn_layout)

        main_layout.addWidget(charmm_group)
        
        # MD File Generation GroupBox
        md_file_group = QGroupBox("MD File Generation")
        md_file_layout = QVBoxLayout(md_file_group)
        
        self._generate_md_config_btn = QPushButton("Generate MD and MM-GBSA Configuration Files")
        self._generate_md_config_btn.clicked.connect(self._generate_md_config_files)
        md_file_layout.addWidget(self._generate_md_config_btn)
        
        self._generate_enhanced_config_btn = QPushButton("Generate Enhanced Sampling Configuration Files for Binding Mode Optimization")
        self._generate_enhanced_config_btn.clicked.connect(self._generate_enhanced_sampling_config_files)
        md_file_layout.addWidget(self._generate_enhanced_config_btn)
        
        self._copy_rename_md_folder_btn = QPushButton("Copy&Paste and Rename MD File Folder")
        self._copy_rename_md_folder_btn.clicked.connect(self._on_copy_and_rename_md_folder_clicked)
        md_file_layout.addWidget(self._copy_rename_md_folder_btn)
        
        main_layout.addWidget(md_file_group)
        
        main_layout.addStretch()

    def copy_and_rename_md_prepper_folder(self, new_name: str) -> pathlib.Path:
        """
        Copy the 'md_prepper' folder and rename the copy to a user-specified name.
        
        Args:
            new_name: The new name for the copied folder (should be a valid folder name).
            
        Returns:
            The path to the newly created folder.
            
        Raises:
            ValueError: If output directory is not set or new_name is invalid.
            FileNotFoundError: If the md_prepper folder does not exist.
            FileExistsError: If a folder with the new name already exists.
        """
        # Validate output directory
        if not self.main_window._output_directory and not self.main_window._output_dir_input.text():
            raise ValueError("Output directory is not set. Please select an output directory first.")
        
        output_dir = pathlib.Path(
            self.main_window._output_directory 
            if self.main_window._output_directory 
            else self.main_window._output_dir_input.text()
        ).resolve()
        
        # Validate new_name
        if not new_name or not new_name.strip():
            raise ValueError("New folder name cannot be empty.")
        
        new_name = new_name.strip()
        
        # Check for invalid characters in folder name
        invalid_chars = ['<', '>', ':', '"', '/', '\\', '|', '?', '*']
        if any(char in new_name for char in invalid_chars):
            raise ValueError(f"Folder name contains invalid characters: {invalid_chars}")
        
        # Define source and target paths
        source_path = output_dir / "md_prepper"
        target_path = output_dir / new_name
        
        # Check if source exists
        if not source_path.exists():
            raise FileNotFoundError(f"The 'md_prepper' folder does not exist at: {source_path}")
        
        # Check if target already exists
        if target_path.exists():
            raise FileExistsError(f"A folder with the name '{new_name}' already exists at: {target_path}")
        
        # Perform the copy operation
        shutil.copytree(source_path, target_path)
        
        return target_path

    def _browse_vmd_path(self):
        """Browse for VMD executable."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select VMD Executable",
            self.main_window._last_dir,
            "Executables (*.exe);;All Files (*.*)"
        )
        
        if file_path:
            self.main_window._vmd_path = file_path
            self._vmd_path_input.setText(file_path)
            self.main_window._save_settings()

    def _save_vmd_path_from_input(self):
        """Save VMD path when manually edited."""
        path = self._vmd_path_input.text().strip()
        if path:
            self.main_window._vmd_path = path
            self.main_window._save_settings()

    def _md_process_protein(self):
        """Process protein using MDPrepper."""
        # Validate output directory
        if not self.main_window._output_directory and not self.main_window._output_dir_input.text():
            self.main_window._show_warning("Error", "Please select an output directory first.")
            return
        
        output_dir = pathlib.Path(self.main_window._output_directory if self.main_window._output_directory else self.main_window._output_dir_input.text()).resolve()
        
        # Validate VMD path
        if not self.main_window._vmd_path:
            self.main_window._show_warning("Error", "Please select VMD executable path.")
            return
            
        # Check source files existence in md_prepper subdirectory
        md_prepper_dir = output_dir / "md_prepper"
        rigid_protein_path = md_prepper_dir / "rigid_protein.pdb"
        ligand_pdb_path = md_prepper_dir / "ligand.pdb"
        
        if not rigid_protein_path.exists():
            self.main_window._show_warning("Error", f"Rigid protein file not found:\n{rigid_protein_path}\n\nPlease run 'Output MD Prepper Files' in Results tab first.")
            return
        
        # Determine flexible residues path if applicable
        flex_res_pdb_path = ""
        docking_type = "rigid"
        if self.main_window._docking_instance:
            last_type = self.main_window._docking_instance.get_last_docking_type()
            if last_type:
                docking_type = last_type
            
        if docking_type == "flexible":
            flex_res_path = md_prepper_dir / "flexible_residues.pdb"
            if flex_res_path.exists():
                flex_res_pdb_path = str(flex_res_path)
                print(f"Using flexible residues from: {flex_res_path}")
            else:
                print(f"Warning: Flexible docking indicated but file not found: {flex_res_path}")
        
        try:
            # Create MDPrepper instance
            prepper = md_prepper.MDPrepper(
                receptor_pdb_path=str(rigid_protein_path),
                ligand_pdb_path=str(ligand_pdb_path) if ligand_pdb_path.exists() else "",
                flex_res_pdb_path=flex_res_pdb_path,
                ligand_smiles="",
                ph=self._md_ph_spinbox.value(),
                output_dir=str(output_dir),
                vmd_path=self.main_window._vmd_path
            )
            
            # Run processing
            final_pdb = prepper.process_protein()
            
            self.main_window._show_info("Success", f"Protein processed successfully!\nOutput: {final_pdb}")
            
        except Exception as e:
            self.main_window._show_error("Error", f"Error processing protein:\n{str(e)}")
            traceback.print_exc()

    def _md_process_ligand(self):
        """Process ligand using MDPrepper (generate MOL2)."""
        # Validate output directory
        if not self.main_window._output_directory and not self.main_window._output_dir_input.text():
            self.main_window._show_warning("Error", "Please select an output directory first.")
            return
            
        output_dir = pathlib.Path(self.main_window._output_directory if self.main_window._output_directory else self.main_window._output_dir_input.text()).resolve()
        
        # Validate VMD path
        if not self.main_window._vmd_path:
            self.main_window._show_warning("Error", "Please select VMD executable path.")
            return
            
        # Check source files
        md_prepper_dir = output_dir / "md_prepper"
        ligand_pdb_path = md_prepper_dir / "ligand.pdb"
        ligand_smi_path = md_prepper_dir / "ligand.smi"
        
        if not ligand_pdb_path.exists():
            self.main_window._show_warning("Error", f"Ligand PDB file not found:\n{ligand_pdb_path}")
            return
            
        if not ligand_smi_path.exists():
            self.main_window._show_warning("Error", f"Ligand SMILES file not found:\n{ligand_smi_path}")
            return
            
        try:
            # Read SMILES
            with open(ligand_smi_path, 'r', encoding='utf-8') as f:
                smiles = f.read().strip()
                
            if not smiles:
                raise ValueError("SMILES file is empty")
            
            # Create MDPrepper instance
            rigid_protein_path = md_prepper_dir / "rigid_protein.pdb"
            
            prepper = md_prepper.MDPrepper(
                receptor_pdb_path=str(rigid_protein_path) if rigid_protein_path.exists() else "",
                ligand_pdb_path=str(ligand_pdb_path),
                flex_res_pdb_path="",
                ligand_smiles=smiles,
                ph=self._md_ph_spinbox.value(),
                output_dir=str(output_dir),
                vmd_path=self.main_window._vmd_path
            )
            
            # Run processing
            mol2_path = prepper.process_ligand()
            
            self.main_window._show_info("Success", f"Ligand processed successfully!\nOutput: {mol2_path}")
            
        except Exception as e:
            self.main_window._show_error("Error", f"Error processing ligand:\n{str(e)}")
            traceback.print_exc()

    def _browse_ligand_str_file(self):
        """Browse for Ligand STR file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Ligand STR File",
            self.main_window._last_dir,
            "Stream Files (*.str);;All Files (*.*)"
        )
        
        if file_path:
            self._ligand_str_edit.setText(file_path)
            self.main_window._last_dir = str(pathlib.Path(file_path).parent.resolve())

    def _generate_charmm_psf(self):
        """Generate CHARMM PSF file."""

        # Validate output directory
        if not self.main_window._output_directory and not self.main_window._output_dir_input.text():
            self.main_window._show_warning("Error", "Please select an output directory first.")
            return

        # Use the directory from input if self._output_directory is not set (e.g. manual entry)
        output_dir_str = self.main_window._output_directory if self.main_window._output_directory else self.main_window._output_dir_input.text()
        output_dir = pathlib.Path(output_dir_str).resolve()
        
        # 1. Define paths
        md_prepper_dir = output_dir / "md_prepper"
        fixed_files_dir = md_prepper_dir / "fixed_files"
        temp_files_dir = md_prepper_dir / "temp_files"
        md_files_dir = md_prepper_dir / "md_files"
        
        # Create output directories
        temp_files_dir.mkdir(parents=True, exist_ok=True)
        md_files_dir.mkdir(parents=True, exist_ok=True)
        
        # Input PDB paths
        protein_pdb_path = fixed_files_dir / "protein_final.pdb"
        ligand_pdb_path = fixed_files_dir / "ligand.pdb"
        
        # Check if files exist
        if not protein_pdb_path.exists():
             self.main_window._show_warning("Error", f"Protein PDB not found at:\n{protein_pdb_path}\nPlease run 'Process Protein' first.")
             return
        
        if not ligand_pdb_path.exists():
             self.main_window._show_warning("Error", f"Ligand PDB not found at:\n{ligand_pdb_path}\nPlease run 'Process Ligand' first.")
             return

        # 2. Get parameters from GUI
        ligand_str_path = self._ligand_str_edit.text().strip()
        padding = self._water_edge_spin.value()
        salt_type = self._ion_type_combo.currentText()
        salt_concentration = self._ion_concentration_spin.value()
        
        # Validate VMD path
        if not self.main_window._vmd_path:
            self.main_window._show_warning("Error", "Please select VMD executable path.")
            return

        try:
            # 3. Initialize MDPrepper
            prepper = md_prepper.MDPrepper(
                receptor_pdb_path=str(protein_pdb_path), # Placeholder
                ligand_pdb_path=str(ligand_pdb_path),
                flex_res_pdb_path="", 
                ligand_smiles="", # Not needed for PSF generation
                ph=7.0, # Not needed here
                output_dir=str(output_dir), # General output dir
                vmd_path=self.main_window._vmd_path
            )

            # 4. Call generate_system_psf
            self.main_window._show_info("Processing", "Generating System PSF... This may take a moment.")
            # Force UI update?
            QApplication.processEvents()

            ionized_psf, ionized_pdb = prepper.generate_system_psf(
                protein_pdb_path=str(protein_pdb_path),
                ligand_pdb_path=str(ligand_pdb_path),
                ligand_str_path=ligand_str_path if ligand_str_path else None,
                padding=padding,
                salt_type=salt_type,
                salt_concentration=salt_concentration,
                temp_output_path=str(temp_files_dir),
                output_path=str(md_files_dir)
            )

            # 5. Copy extra files
            third_party_dir = md_prepper.third_party_tools.get_third_party_dir()
            force_field_dir = third_party_dir / "force_field"
            
            files_to_copy = [
                "par_all36m_prot.prm",
                "par_all36_cgenff.prm",
                "toppar_water_ions.str"
            ]
            
            for file_name in files_to_copy:
                src = force_field_dir / file_name
                dst = md_files_dir / file_name
                if src.exists():
                    shutil.copy(src, dst)
                    print(f"Copied {file_name} to {md_files_dir}")
                else:
                    print(f"Warning: Force field file not found: {src}")

            # Copy ligand.str if it exists in temp_files
            ligand_str_src = (temp_files_dir / "ligand.str").resolve()
            if ligand_str_src.exists():
                target_str = (md_files_dir / "ligand.str").resolve()
                shutil.copy(ligand_str_src, target_str)
                print(f"Copied ligand.str to {md_files_dir}")

            # 6. Success message
            self.main_window._show_info("Success", 
                f"System PSF generated successfully!\n\n"
                f"Output Files located in:\n{md_files_dir}\n"
                f"  - {pathlib.Path(ionized_psf).name}\n"
                f"  - {pathlib.Path(ionized_pdb).name}")

        except Exception as e:
            self.main_window._show_error("Error", f"Error generating PSF:\n{str(e)}")
            traceback.print_exc()

    def _generate_md_config_files(self):
        """Generate MD configuration files (copy OpenMM scripts)."""
        # Validate output directory
        if not self.main_window._output_directory and not self.main_window._output_dir_input.text():
            self.main_window._show_warning("Error", "Please select an output directory first.")
            return

        # Use the directory from input if self._output_directory is not set
        output_dir_str = self.main_window._output_directory if self.main_window._output_directory else self.main_window._output_dir_input.text()
        output_dir = pathlib.Path(output_dir_str).resolve()
        
        # Define md_files directory inside md_prepper
        md_files_dir = output_dir / "md_prepper" / "md_files"
        
        # Validate VMD path (required for MDPrepper init)
        if not self.main_window._vmd_path:
            self.main_window._show_warning("Error", "Please select VMD executable path.")
            return

        try:
            # Initialize MDPrepper (using dummy paths for unused arguments)
            prepper = md_prepper.MDPrepper(
                receptor_pdb_path="",
                ligand_pdb_path="",
                flex_res_pdb_path="", 
                ligand_smiles="",
                ph=7.0,
                output_dir=str(output_dir),
                vmd_path=self.main_window._vmd_path
            )
            
            # Call copy_openmm_scripts
            prepper.copy_openmm_scripts(str(md_files_dir))
            
            self.main_window._show_info("Success", 
                f"MD configuration files generated successfully!\n\n"
                f"Files copied to:\n{md_files_dir}")
                
        except Exception as e:
            self.main_window._show_error("Error", f"Error generating MD configuration files:\n{str(e)}")
            traceback.print_exc()

    def _generate_enhanced_sampling_config_files(self):
        """Generate Enhanced Sampling configuration files (copy WTM-TABF scripts)."""
        # Validate output directory
        if not self.main_window._output_directory and not self.main_window._output_dir_input.text():
            self.main_window._show_warning("Error", "Please select an output directory first.")
            return

        # Use the directory from input if self._output_directory is not set
        output_dir_str = self.main_window._output_directory if self.main_window._output_directory else self.main_window._output_dir_input.text()
        output_dir = pathlib.Path(output_dir_str).resolve()
        
        # Define md_files directory inside md_prepper
        md_files_dir = output_dir / "md_prepper" / "md_files"
        
        # Validate VMD path (required for MDPrepper init)
        if not self.main_window._vmd_path:
            self.main_window._show_warning("Error", "Please select VMD executable path.")
            return

        try:
            # Initialize MDPrepper (using dummy paths for unused arguments)
            prepper = md_prepper.MDPrepper(
                receptor_pdb_path="",
                ligand_pdb_path="",
                flex_res_pdb_path="", 
                ligand_smiles="",
                ph=7.0,
                output_dir=str(output_dir),
                vmd_path=self.main_window._vmd_path
            )
            
            # Call copy_openmm_scripts with enhanced_sampling=True
            prepper.copy_openmm_scripts(str(md_files_dir), enhanced_sampling=True)
            
            self.main_window._show_info("Success", 
                f"Enhanced Sampling configuration files generated successfully!\n\n"
                f"Files copied to:\n{md_files_dir}")
                
        except Exception as e:
            self.main_window._show_error("Error", f"Error generating Enhanced Sampling configuration files:\n{str(e)}")
            traceback.print_exc()

    def _on_copy_and_rename_md_folder_clicked(self):
        """Handle the Copy&Paste and Rename MD File Folder button click."""
        # Show input dialog to get new folder name
        new_name, ok = QInputDialog.getText(
            self,
            "Copy&Paste and Rename MD File Folder",
            "Enter the new folder name for the copy:",
            QLineEdit.EchoMode.Normal,
            "md_prepper_copy"
        )
        
        if not ok:
            return  # User cancelled
        
        try:
            new_path = self.copy_and_rename_md_prepper_folder(new_name)
            self.main_window._show_info("Success", f"Folder copied successfully!\n\nNew path:\n{new_path}")
        except ValueError as e:
            self.main_window._show_warning("Invalid Input", str(e))
        except FileNotFoundError as e:
            self.main_window._show_warning("Folder Not Found", str(e))
        except FileExistsError as e:
            self.main_window._show_warning("Folder Already Exists", str(e))
        except Exception as e:
            self.main_window._show_error("Error", f"Error copying folder:\n{str(e)}")
            traceback.print_exc()
