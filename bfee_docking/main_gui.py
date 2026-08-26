# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Haohao Fu <fhh2626@nankai.edu.cn, fhh2626@gmail.com>

# Standard library imports
import pathlib
import sys
import traceback
import configparser
from platformdirs import user_config_dir

# PySide6 imports
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QFileDialog, QTabWidget,
    QMessageBox, QComboBox
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QIcon

# Project module imports
from . import __version__
from .gui_utils import get_app_icon
from .tabs.protein_tab import ProteinTab
from .tabs.ligand_tab import LigandTab
from .tabs.docking_tab import DockingTab
from .tabs.results_tab import ResultsTab
from .tabs.md_prepper_tab import MDPrepperTab


def _get_settings_path() -> pathlib.Path:
    """Return the per-user settings path."""
    settings_dir = pathlib.Path(user_config_dir("BFEE-Docking", "BFEE-Docking"))
    settings_dir.mkdir(parents=True, exist_ok=True)
    return settings_dir / "setting.ini"

class MainWindow(QMainWindow):
    """Main application window."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"BFEE-Docking {__version__}")
        
        # Set window icon
        self.setWindowIcon(get_app_icon())
        
        # State variables (protected)
        self._output_directory = None
        self._protein_file = None
        self._protein_parser = None
        self._ligand_files = []
        self._ligand_parser = None
        self._ligand_processing_worker = None
        self._ligand_smiles_list = []  # Store processed SMILES
        self._ligand_pdbqt_files = []  # Store processed PDBQT files
        self._docking_instance = None
        self._docking_worker = None
        self._detected_pockets = []  # Store detected pocket regions from pocketeer
        self._vmd_path = None
        self._saved_theme = None  # Store the saved theme name
        
        # Remember last used directory for better UX (shared across all browse buttons)
        self._last_dir = str(pathlib.Path.home())
        
        # Load settings
        self._load_settings()
        
        self._setup_ui()
        
        # Set the initial size to the minimum possible to fit the content
        self.adjustSize()
    
    def _show_info(self, title: str, message: str):
        """Show an information message box with the application icon."""
        msg_box = QMessageBox(self)
        msg_box.setWindowTitle(title)
        msg_box.setText(message)
        msg_box.setIcon(QMessageBox.Information)
        msg_box.setWindowIcon(get_app_icon())
        msg_box.exec()
    
    def _show_warning(self, title: str, message: str):
        """Show a warning message box with a single OK button."""
        msg_box = QMessageBox(self)
        msg_box.setWindowTitle(title)
        msg_box.setText(message)
        msg_box.setIcon(QMessageBox.Warning)
        msg_box.setWindowIcon(get_app_icon())
        msg_box.setStandardButtons(QMessageBox.Ok)
        msg_box.setDefaultButton(QMessageBox.Ok)
        msg_box.exec()
    
    def _show_error(self, title: str, message: str):
        """Show an error/critical message box with the application icon."""
        msg_box = QMessageBox(self)
        msg_box.setWindowTitle(title)
        msg_box.setText(message)
        msg_box.setIcon(QMessageBox.Critical)
        msg_box.setWindowIcon(get_app_icon())
        msg_box.exec()
        
    def _setup_ui(self):
        """Set up the main UI."""
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        
        # Output Directory Section
        output_dir_layout = QHBoxLayout()
        
        output_dir_label = QLabel("Output Directory:")
        output_dir_layout.addWidget(output_dir_label)
        
        self._output_dir_input = QLineEdit()
        self._output_dir_input.setPlaceholderText("Select output directory...")
        output_dir_layout.addWidget(self._output_dir_input)
        
        browse_output_dir_btn = QPushButton("Browse")
        browse_output_dir_btn.clicked.connect(self._browse_output_directory)
        output_dir_layout.addWidget(browse_output_dir_btn)
        
        main_layout.addLayout(output_dir_layout)
        
        # Tab Widget
        self.tabs = QTabWidget()
        
        # Instantiate Tab Classes
        self._protein_tab = ProteinTab(self)
        self._ligand_tab = LigandTab(self)
        self._docking_tab = DockingTab(self)
        self._results_tab = ResultsTab(self)
        self._md_prepper_tab = MDPrepperTab(self)
        
        # Add Tabs
        self.tabs.addTab(self._protein_tab, "Protein")
        self.tabs.addTab(self._ligand_tab, "Ligand")
        self.tabs.addTab(self._docking_tab, "Docking")
        self.tabs.addTab(self._results_tab, "Results and Visualization")
        self.tabs.addTab(self._md_prepper_tab, "MD Prepper")
        
        main_layout.addWidget(self.tabs)
        
        # Theme Selection Section
        theme_layout = QHBoxLayout()
        
        theme_label = QLabel("Theme:")
        theme_layout.addWidget(theme_label)
        
        self._theme_combo = QComboBox()
        self._load_themes()
        self._theme_combo.currentTextChanged.connect(self._apply_theme)
        theme_layout.addWidget(self._theme_combo)
        
        # Add stretch to push theme selector to the left
        theme_layout.addStretch()
        
        main_layout.addLayout(theme_layout)
        
        # Apply saved theme
        if self._saved_theme and self._theme_combo.findText(self._saved_theme) >= 0:
            self._theme_combo.setCurrentText(self._saved_theme)
        else:
            self._apply_theme(self._theme_combo.currentText())
        
    def _browse_output_directory(self):
        """Browse for output directory."""
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Output Directory",
            self._last_dir
        )
        
        if directory:
            self._output_directory = directory
            self._output_dir_input.setText(directory)
            self._last_dir = directory

    def _load_settings(self):
        """Load settings from setting.ini."""
        config = configparser.ConfigParser()
        settings_path = _get_settings_path()
        
        if settings_path.exists():
            try:
                config.read(settings_path)
                if "General" in config:
                    if "vmd_path" in config["General"]:
                        self._vmd_path = config["General"]["vmd_path"]
                    if "theme" in config["General"]:
                        self._saved_theme = config["General"]["theme"]
            except Exception as e:
                print(f"Error loading settings: {e}")

    def _save_settings(self):
        """Save settings to setting.ini."""
        config = configparser.ConfigParser()
        settings_path = _get_settings_path()
        
        # Read existing config to preserve other settings if any
        if settings_path.exists():
            config.read(settings_path)
        
        if "General" not in config:
            config["General"] = {}
            
        if self._vmd_path:
            config["General"]["vmd_path"] = self._vmd_path
            
        if hasattr(self, '_theme_combo') and self._theme_combo.currentText():
            config["General"]["theme"] = self._theme_combo.currentText()
            
        try:
            with open(settings_path, 'w') as configfile:
                config.write(configfile)
        except Exception as e:
            print(f"Error saving settings: {e}")
    
    def _load_themes(self):
        """Load available themes from the styles folder."""
        styles_dir = pathlib.Path(__file__).parent / "styles"
        self._theme_combo.addItem("Default")  # Add default (no theme) option
        
        if styles_dir.exists():
            for qss_file in sorted(styles_dir.glob("*.qss")):
                # Convert filename to display name (e.g., modern_dark.qss -> Modern Dark)
                theme_name = qss_file.stem.replace("_", " ").title()
                self._theme_combo.addItem(theme_name, str(qss_file))
    
    def _apply_theme(self, theme_name: str):
        """Apply the selected theme."""
        if theme_name == "Default":
            QApplication.instance().setStyleSheet("")
        else:
            # Get the file path from combo box data
            index = self._theme_combo.findText(theme_name)
            if index >= 0:
                qss_path = self._theme_combo.itemData(index)
                if qss_path and pathlib.Path(qss_path).exists():
                    try:
                        with open(qss_path, 'r', encoding='utf-8') as f:
                            stylesheet = f.read()
                        QApplication.instance().setStyleSheet(stylesheet)
                    except Exception as e:
                        print(f"Error loading theme: {e}")
                        QApplication.instance().setStyleSheet("")
        
        # Save the theme preference
        self._save_settings()

def main():
    """Main application entry point."""
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
