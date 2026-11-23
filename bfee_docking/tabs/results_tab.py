import math
import pathlib
import traceback

from PySide6.QtCore import QSignalBlocker, Qt
from PySide6.QtGui import QPalette
from PySide6.QtWidgets import (
    QComboBox, QDoubleSpinBox, QFileDialog, QGroupBox, QHBoxLayout, QLabel,
    QLineEdit, QProxyStyle, QPushButton, QStyle, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)


class _PoseComboStyle(QProxyStyle):
    def styleHint(self, hint, option=None, widget=None, returnData=None):
        # A menu-style popup ignores maxVisibleItems and uses scroll arrows.
        if hint == QStyle.SH_ComboBox_Popup:
            return 0
        return super().styleHint(hint, option, widget, returnData)


class ResultsTab(QWidget):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self._pose_affinities = {}
        self._result_rows = {}
        self._pose_molecule = None
        self._setup_ui()
        self.update_results_table()

    def _setup_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.addWidget(self._create_results_group(), 1)
        main_layout.addWidget(self._create_pose_group())
        main_layout.addWidget(self._create_redocking_group())

        self._molecule_input.textChanged.connect(self._refresh_poses)
        self._pose_combo.currentIndexChanged.connect(self._update_redocking_preview)
        self._box_expansion_spinbox.valueChanged.connect(self._update_redocking_preview)

    def _create_results_group(self):
        group = QGroupBox("Results")
        layout = QVBoxLayout(group)

        self._results_table = QTableWidget()
        self._results_table.setColumnCount(6)
        self._results_table.setHorizontalHeaderLabels([
            "Number", "SMILES", "Best Affinity", "Affinity", "Molecule Number", "Filename"
        ])
        self._results_table.setSelectionBehavior(QTableWidget.SelectRows)
        self._results_table.setSelectionMode(QTableWidget.SingleSelection)
        self._results_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self._results_table.itemSelectionChanged.connect(self._on_table_selection_changed)
        layout.addWidget(self._results_table)

        actions = QHBoxLayout()
        self._sort_results_btn = QPushButton("Sort by Best Affinity")
        self._sort_results_btn.clicked.connect(self._sort_results_table)
        actions.addWidget(self._sort_results_btn)

        actions.addWidget(QLabel("Cluster threshold"))
        self._similarity_threshold_input = QDoubleSpinBox()
        self._similarity_threshold_input.setRange(0.01, 1.0)
        self._similarity_threshold_input.setValue(0.70)
        self._similarity_threshold_input.setSingleStep(0.1)
        self._similarity_threshold_input.setDecimals(2)
        actions.addWidget(self._similarity_threshold_input)

        self._cluster_btn = QPushButton("Cluster and Save CSV...")
        self._cluster_btn.clicked.connect(self._cluster_results)
        actions.addWidget(self._cluster_btn)
        actions.addStretch()
        layout.addLayout(actions)
        return group

    def _create_pose_group(self):
        group = QGroupBox("Selected Pose")
        layout = QVBoxLayout(group)
        selection = QHBoxLayout()

        selection.addWidget(QLabel("Molecule"))
        self._molecule_input = QLineEdit("0")
        self._molecule_input.setMaximumWidth(90)
        selection.addWidget(self._molecule_input)

        selection.addWidget(QLabel("Pose"))
        self._pose_combo = QComboBox()
        self._pose_combo_style = _PoseComboStyle(self._pose_combo.style().objectName())
        self._pose_combo_style.setParent(self._pose_combo)
        self._pose_combo.setStyle(self._pose_combo_style)
        self._pose_combo.setMaxVisibleItems(10)
        self._pose_combo.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self._pose_combo.setMinimumContentsLength(22)
        self._pose_combo.view().setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        selection.addWidget(self._pose_combo)

        self._pose_hint = QLabel()
        self._pose_hint.setForegroundRole(QPalette.PlaceholderText)
        selection.addWidget(self._pose_hint)
        selection.addStretch()
        layout.addLayout(selection)

        actions = QHBoxLayout()
        self._show_in_pymol_btn = QPushButton("Show in PyMOL")
        self._show_in_pymol_btn.clicked.connect(self._show_in_pymol)
        actions.addWidget(self._show_in_pymol_btn)
        self._output_md_prepper_btn = QPushButton("Export to MD Prepper")
        self._output_md_prepper_btn.clicked.connect(self._output_md_prepper_files)
        actions.addWidget(self._output_md_prepper_btn)
        actions.addStretch()
        layout.addLayout(actions)
        return group

    def _create_redocking_group(self):
        group = QGroupBox("Redocking")
        layout = QVBoxLayout(group)

        expansion = QHBoxLayout()
        expansion.addWidget(QLabel("Box expansion (Å)"))
        self._box_expansion_spinbox = QDoubleSpinBox()
        self._box_expansion_spinbox.setRange(0.0, 100.0)
        self._box_expansion_spinbox.setDecimals(1)
        self._box_expansion_spinbox.setSingleStep(0.5)
        self._box_expansion_spinbox.setValue(6.0)
        self._box_expansion_spinbox.setToolTip(
            "Adds this amount to each box dimension. Each side expands by half."
        )
        expansion.addWidget(self._box_expansion_spinbox)
        expansion.addStretch()
        layout.addLayout(expansion)

        preview = QHBoxLayout()
        preview.addWidget(QLabel("Center (Å)"))
        self._redocking_center_input = QLineEdit()
        self._redocking_center_input.setReadOnly(True)
        self._redocking_center_input.setPlaceholderText("x, y, z")
        preview.addWidget(self._redocking_center_input)
        preview.addWidget(QLabel("Size (Å)"))
        self._redocking_size_input = QLineEdit()
        self._redocking_size_input.setReadOnly(True)
        self._redocking_size_input.setPlaceholderText("x, y, z")
        preview.addWidget(self._redocking_size_input)
        layout.addLayout(preview)

        actions = QHBoxLayout()
        self._apply_redocking_btn = QPushButton("Apply to Protein Tab")
        self._apply_redocking_btn.clicked.connect(self._apply_redocking_range)
        actions.addWidget(self._apply_redocking_btn)
        self._redocking_hint = QLabel()
        self._redocking_hint.setTextFormat(Qt.PlainText)
        self._redocking_hint.setWordWrap(True)
        actions.addWidget(self._redocking_hint, 1)
        layout.addLayout(actions)
        return group

    def _molecule_index(self):
        try:
            return int(self._molecule_input.text())
        except ValueError:
            return None

    def _on_table_selection_changed(self):
        selected = self._results_table.selectedRanges()
        if selected:
            item = self._results_table.item(selected[0].topRow(), 4)
            if item and int(item.text()) != self._molecule_index():
                self._molecule_input.setText(item.text())

    def update_results_table(self):
        """Load results and refresh the available poses for the selected molecule."""
        docking = self.main_window._docking_instance
        results = docking.get_combined_results() if docking else []
        self._pose_affinities = {row['molecule_number']: row['affinity'] for row in results}
        self._pose_molecule = None
        self._populate_results_table_rows(results)

        has_results = bool(results)
        self._molecule_input.setEnabled(has_results)
        self._sort_results_btn.setEnabled(has_results)
        self._similarity_threshold_input.setEnabled(has_results)
        self._cluster_btn.setEnabled(has_results)
        molecule = self._molecule_index()
        if results and not self._pose_affinities.get(molecule):
            first_available = next(
                (index for index, poses in self._pose_affinities.items() if poses),
                results[0]['molecule_number'],
            )
            with QSignalBlocker(self._molecule_input):
                self._molecule_input.setText(str(first_available))
        self._refresh_poses()

    def _populate_results_table_rows(self, results: list[dict]) -> None:
        self._result_rows = {result['molecule_number']: row for row, result in enumerate(results)}
        with QSignalBlocker(self._results_table):
            self._results_table.setRowCount(0)
            self._results_table.setRowCount(len(results))

        for row, result in enumerate(results):
            best = result['best_affinity']
            affinities = ", ".join(f"{affinity:.2f}" for affinity in result['affinity'])
            values = (
                str(result['number']), result['smiles'],
                f"{best:.2f}" if isinstance(best, (int, float)) else str(best),
                f"[{affinities}]", str(result['molecule_number']), result.get('filename', ''),
            )
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                if column in (0, 2, 4):
                    item.setTextAlignment(Qt.AlignCenter)
                self._results_table.setItem(row, column, item)

        self._results_table.resizeColumnsToContents()
        self._select_molecule_row(self._molecule_index())

    def _select_molecule_row(self, molecule):
        row = self._result_rows.get(molecule)
        with QSignalBlocker(self._results_table):
            if row is None:
                self._results_table.clearSelection()
            else:
                selected = self._results_table.selectedRanges()
                if not selected or selected[0].topRow() != row:
                    self._results_table.selectRow(row)

    def _refresh_poses(self):
        molecule = self._molecule_index()
        affinities = self._pose_affinities.get(molecule, [])
        selected_pose = self._pose_combo.currentData() if molecule == self._pose_molecule else None
        self._pose_molecule = molecule
        self._select_molecule_row(molecule)
        with QSignalBlocker(self._pose_combo):
            self._pose_combo.clear()
            for index, affinity in enumerate(affinities):
                self._pose_combo.addItem(f"{index} ({affinity:.2f} kcal/mol)", index)
            if affinities:
                index = self._pose_combo.findData(selected_pose) if selected_pose is not None else -1
                self._pose_combo.setCurrentIndex(max(0, index))

        if not self._pose_affinities:
            hint = "Run docking or load results to select a pose."
        elif molecule not in self._pose_affinities:
            hint = "Enter a molecule number from the Results table."
        elif not affinities:
            hint = "This molecule has no poses."
        else:
            hint = ""
        self._pose_hint.setText(hint)
        self._pose_hint.setVisible(bool(hint))
        for widget in (self._pose_combo, self._show_in_pymol_btn,
                       self._output_md_prepper_btn, self._box_expansion_spinbox):
            widget.setEnabled(bool(affinities))
        self._update_redocking_preview()

    def _selected_pose(self):
        molecule = self._molecule_index()
        pose = self._pose_combo.currentData()
        if self.main_window._docking_instance is not None and pose is not None:
            if 0 <= pose < len(self._pose_affinities.get(molecule, [])):
                return molecule, pose
        return None

    def _show_in_pymol(self):
        selection = self._selected_pose()
        if selection is None:
            return
        try:
            self.main_window._docking_instance.open_results(*selection)
        except Exception as e:
            self.main_window._show_error("Error", f"Error opening PyMOL:\n{e}")

    def _output_md_prepper_files(self):
        selection = self._selected_pose()
        if selection is None:
            return
        molecule, pose = selection
        try:
            output_text = self.main_window._output_dir_input.text().strip()
            output_dir = pathlib.Path(output_text)
            if not output_text or not output_dir.is_dir():
                self.main_window._show_warning("Error", "Please specify a valid output directory.")
                return

            md_prepper_dir = output_dir / "md_prepper"
            md_prepper_dir.mkdir(parents=True, exist_ok=True)
            result_files = self.main_window._docking_instance.output_MD_prepper_files(
                ligand_index=molecule, pose_index=pose, output_dir=md_prepper_dir
            )
            names = "\n".join(f"  - {path.name}" for path in result_files.values())
            self.main_window._show_info(
                "Success", f"MD prepper files saved to:\n{md_prepper_dir}\n\nGenerated files:\n{names}"
            )
        except Exception as e:
            self.main_window._show_error("Error", f"Error outputting MD prepper files:\n{e}")
            traceback.print_exc()

    def _calculate_redocking_box(self, selection):
        center, lengths = self.main_window._docking_instance.get_redocking_range(*selection)
        expansion = self._box_expansion_spinbox.value()
        size = [length + expansion for length in lengths]
        if not all(math.isfinite(value) for value in (*center, *size)):
            raise ValueError("The pose contains invalid coordinates.")
        if any(length <= 0.0 for length in size):
            raise ValueError("Increase Box expansion so that all box dimensions are positive.")
        return center, size

    @staticmethod
    def _format_coordinates(values):
        return ", ".join(f"{value:.4f}".rstrip("0").rstrip(".") for value in values)

    def _display_redocking_box(self, box=None, hint=""):
        self._apply_redocking_btn.setEnabled(box is not None)
        if box is not None:
            center, size = box
            self._redocking_center_input.setText(self._format_coordinates(center))
            self._redocking_size_input.setText(self._format_coordinates(size))
            self._redocking_hint.setText("Sets Specify Manually and switches to the Protein tab.")
        else:
            self._redocking_center_input.clear()
            self._redocking_size_input.clear()
            self._redocking_hint.setText(hint)

    def _update_redocking_preview(self):
        selection = self._selected_pose()
        if selection is None:
            self._display_redocking_box()
            return
        try:
            box = self._calculate_redocking_box(selection)
        except (OSError, ValueError, IndexError, RuntimeError) as e:
            self._display_redocking_box(hint=f"Cannot calculate the box: {e}")
        else:
            self._display_redocking_box(box)

    def _apply_redocking_range(self):
        selection = self._selected_pose()
        if selection is None:
            return
        protein_tab = self.main_window._protein_tab
        try:
            box = self._calculate_redocking_box(selection)
        except (OSError, ValueError, IndexError, RuntimeError) as e:
            self._display_redocking_box(hint=f"Cannot calculate the box: {e}")
            return
        try:
            protein_tab.apply_manual_docking_box(*box)
        except (ValueError, RuntimeError) as e:
            self.main_window._show_error("Error", f"Failed to apply redocking range:\n{e}")
            return
        self._display_redocking_box(box)
        self.main_window.tabs.setCurrentWidget(protein_tab)

    def _sort_results_table(self):
        docking = self.main_window._docking_instance
        if docking:
            self._populate_results_table_rows(docking.get_sorted_combined_results())

    def _cluster_results(self):
        """Cluster results and save to CSV."""
        docking = self.main_window._docking_instance
        if not docking or not docking.get_combined_results():
            self.main_window._show_warning("Error", "No results to cluster.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Clustered Results", self.main_window._last_dir, "CSV Files (*.csv)"
        )
        if not file_path:
            return

        try:
            docking.cluster_results(similarity_threshold=self._similarity_threshold_input.value())
            docking.save_results_to_csv(file_path)
            self.main_window._show_info("Success", f"Results clustered and saved to:\n{file_path}")
        except Exception as e:
            self.main_window._show_error("Error", f"Error during clustering/saving:\n{e}")
            print(f"Error: {e}")
