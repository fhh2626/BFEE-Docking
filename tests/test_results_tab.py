import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from PySide6.QtCore import Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog, QGroupBox, QLineEdit, QMessageBox, QTabWidget

from bfee_docking.docking import Docking
from bfee_docking.tabs.docking_tab import DockingTab
from bfee_docking.tabs.protein_tab import ProteinTab
from bfee_docking.tabs.results_tab import ResultsTab


class ResultsTabTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.rows = [
            {"number": 1, "molecule_number": 0, "smiles": "CC", "affinity": [-7.23, -6.73],
             "best_affinity": -7.23},
            {"number": 2, "molecule_number": 1, "smiles": "CCO", "affinity": [-9.21, -8.0, -7.1],
             "best_affinity": -9.21},
            {"number": 3, "molecule_number": 2, "smiles": "O", "affinity": [],
             "best_affinity": "N/A"},
        ]
        self.docking = Mock()
        self.docking.get_combined_results.return_value = self.rows
        self.docking.get_sorted_combined_results.return_value = [self.rows[1], self.rows[0], self.rows[2]]
        self.docking.get_redocking_range.side_effect = (
            lambda molecule, pose: ([molecule + 0.125, pose + 0.0005, 3.0], [4.0, 5.0, 6.0])
        )
        output = TemporaryDirectory()
        self.addCleanup(output.cleanup)
        self.output = Path(output.name)
        self.window = SimpleNamespace(
            _docking_instance=self.docking,
            _docking_worker=None,
            _protein_parser=Mock(),
            _output_directory=output.name,
            _ligand_pdbqt_files=[self.output / 'ligand.pdbqt', self.output / 'second.pdbqt'],
            _ligand_smiles_list=['CC', 'CCO'],
            _output_dir_input=QLineEdit(output.name),
            _last_dir=output.name,
            _show_warning=Mock(), _show_error=Mock(), _show_info=Mock(),
            tabs=QTabWidget(),
        )
        self.window._protein_parser.get_generated_pdbqt_path.return_value = self.output / 'protein.pdbqt'
        self.window._protein_parser.get_generated_docking_range.return_value = [[0, 0, 0], [20, 20, 20]]
        self.protein = self.window._protein_tab = ProteinTab(self.window)
        self.docking_tab = self.window._docking_tab = DockingTab(self.window)
        self.results = self.window._results_tab = ResultsTab(self.window)
        self.window.tabs.addTab(self.protein, "Protein")
        self.window.tabs.addTab(self.results, "Results and Visualization")
        self.window.tabs.setCurrentWidget(self.results)
        self.addCleanup(self.window.tabs.close)

    def test_results_controls_are_grouped_and_sort_keeps_molecule_id(self):
        group = self.results._results_table.parentWidget()
        self.assertEqual(group.title(), "Results")
        self.assertIs(self.results._sort_results_btn.parentWidget(), group)
        self.assertIs(self.results._cluster_btn.parentWidget(), group)
        self.assertEqual([box.title() for box in self.results.findChildren(QGroupBox)],
                         ["Results", "Selected Pose", "Redocking"])
        self.results._sort_results_btn.click()
        self.assertEqual(self.results._results_table.item(0, 4).text(), "1")
        self.assertEqual(self.results._molecule_input.text(), "0")
        self.results._results_table.selectRow(0)
        self.assertEqual(self.results._molecule_input.text(), "1")
        self.assertEqual(self.results._pose_combo.count(), 3)
        self.assertEqual(self.results._pose_combo.currentText(), "0 (-9.21 kcal/mol)")
        self.docking.get_redocking_range.assert_called_with(1, 0)
        with patch.object(QFileDialog, "getSaveFileName", return_value=(str(self.output / "results.csv"), "")):
            self.results._cluster_btn.click()
        self.docking.cluster_results.assert_called_once_with(similarity_threshold=0.70)
        self.docking.save_results_to_csv.assert_called_once_with(str(self.output / "results.csv"))

    def test_keyboard_row_selection_and_same_row_click_keep_actions_in_sync(self):
        table = self.results._results_table
        self.window.tabs.show()
        table.setFocus()
        QTest.keyClick(table, Qt.Key_Down)
        self.assertEqual(self.results._molecule_input.text(), "1")
        self.results._pose_combo.setCurrentIndex(2)
        cell = table.visualItemRect(table.item(1, 0)).center()
        QTest.mouseClick(table.viewport(), Qt.LeftButton, pos=cell)
        self.assertEqual(self.results._pose_combo.currentData(), 2)
        self.results._sort_results_btn.click()
        self.results._show_in_pymol_btn.click()
        self.docking.open_results.assert_called_once_with(1, 2)

    def test_typing_molecule_refreshes_poses_and_disables_invalid_selection(self):
        field = self.results._molecule_input
        field.selectAll()
        QTest.keyClicks(field, "1")
        self.assertEqual(self.results._pose_combo.count(), 3)
        self.assertEqual(self.results._results_table.currentRow(), 1)
        self.results._pose_combo.setCurrentIndex(2)
        self.assertEqual(self.results._pose_combo.currentData(), 2)
        self.docking.get_redocking_range.assert_called_with(1, 2)
        for invalid in ("", "abc", "-1", "100", "2"):
            with self.subTest(invalid=invalid):
                field.setText(invalid)
                self.assertEqual(self.results._pose_combo.count(), 0)
                self.assertFalse(self.results._show_in_pymol_btn.isEnabled())
                self.assertFalse(self.results._output_md_prepper_btn.isEnabled())
                self.assertFalse(self.results._apply_redocking_btn.isEnabled())
                self.assertEqual(self.results._redocking_center_input.text(), "")
                self.assertEqual(self.results._redocking_size_input.text(), "")
                self.assertTrue(self.results._pose_hint.text())
        field.setText("0")
        self.assertEqual(self.results._pose_combo.currentData(), 0)
        self.assertTrue(self.results._apply_redocking_btn.isEnabled())
        self.assertFalse(self.results._pose_hint.text())

    def test_long_pose_menu_can_scroll_and_actions_use_selected_pose(self):
        self.rows[0]["affinity"] = [-7.23] * 50
        self.results.update_results_table()
        combo = self.results._pose_combo
        self.window.tabs.show()
        combo.showPopup()
        self.app.processEvents()
        self.assertEqual(combo.count(), 50)
        self.assertGreater(combo.view().verticalScrollBar().maximum(), 0)
        combo.hidePopup()
        combo.setCurrentIndex(49)
        self.results._show_in_pymol_btn.click()
        self.docking.open_results.assert_called_once_with(0, 49)
        md_dir = self.output / "md_prepper"
        self.docking.output_MD_prepper_files.return_value = {"ligand": md_dir / "ligand.pdb"}
        self.results._output_md_prepper_btn.click()
        self.docking.output_MD_prepper_files.assert_called_once_with(ligand_index=0, pose_index=49, output_dir=md_dir)
        self.window._show_error.assert_not_called()

    def test_apply_recalculates_box_even_when_manual_mode_is_already_selected(self):
        self.results._molecule_input.setText("1")
        self.results._pose_combo.setCurrentIndex(2)
        self.results._box_expansion_spinbox.setValue(8)
        self.assertEqual(self.results._redocking_center_input.text(), "1.125, 2.0005, 3")
        self.assertEqual(self.results._redocking_size_input.text(), "12, 13, 14")
        self.protein._flex_result_input.setText("A:123")
        self.results._redocking_size_input.setText("stale display")
        self.results._apply_redocking_btn.click()
        parser = self.window._protein_parser
        parser.set_custom_docking_range.assert_called_once_with(center=[1.125, 2.0005, 3.0], side_length=[12, 13, 14])
        self.assertEqual(self.results._redocking_size_input.text(), "12, 13, 14")
        self.assertEqual(self.window.tabs.currentWidget(), self.protein)
        self.assertTrue(self.protein._radio_specify_manually.isChecked())
        self.assertTrue(self.protein._center_y_input.isEnabled())
        self.assertEqual(self.protein._center_y_input.text(), "2.0005")
        self.assertFalse(self.protein._flex_result_input.text())
        self.assertIn("Detect flexible residues again", self.protein._flex_result_input.placeholderText())
        self.window.tabs.setCurrentWidget(self.results)
        self.results._box_expansion_spinbox.setValue(10)
        self.results._apply_redocking_btn.click()
        parser.set_custom_docking_range.assert_called_with(center=[1.125, 2.0005, 3.0], side_length=[14, 15, 16])
        self.assertEqual(parser.set_custom_docking_range.call_count, 2)
        self.assertEqual(self.protein._length_x_input.text(), "14")
        self.assertEqual(self.window.tabs.currentWidget(), self.protein)
        self.window._show_info.assert_not_called()

    def test_pose_bounds_from_pdbqt_and_failed_preview(self):
        pdbqt = self.output / "ligand.pdbqt"
        def atom(i, x, y, z):
            return f"ATOM  {i:5d}  C   LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}\n"
        pdbqt.write_text(
            "MODEL 1\n" + atom(1, 1, 2, 3) + atom(2, 5, 2, 9) + "ENDMDL\n"
            "MODEL 2\n" + atom(1, 25, 25, 25) + "ENDMDL\n", encoding="ascii"
        )
        docking = object.__new__(Docking)
        docking._docking_results = [{"output": pdbqt}]
        self.docking.get_redocking_range.side_effect = docking.get_redocking_range
        self.results.update_results_table()
        self.assertEqual(self.results._redocking_center_input.text(), "3, 2, 6")
        self.assertEqual(self.results._redocking_size_input.text(), "10, 6, 12")
        self.results._box_expansion_spinbox.setValue(0)
        self.assertFalse(self.results._apply_redocking_btn.isEnabled())
        self.assertIn("positive", self.results._redocking_hint.text())
        self.results._box_expansion_spinbox.setValue(6)
        self.assertTrue(self.results._apply_redocking_btn.isEnabled())
        pdbqt.unlink()
        self.results._apply_redocking_btn.click()
        self.assertEqual(self.window.tabs.currentWidget(), self.results)
        self.assertFalse(self.results._apply_redocking_btn.isEnabled())
        self.assertEqual(self.results._redocking_center_input.text(), "")
        self.assertIn("PDBQT file not found", self.results._redocking_hint.text())
        self.window._protein_parser.set_custom_docking_range.assert_not_called()

    def test_old_flexible_residues_require_review_after_applying_box(self):
        tab = self.docking_tab
        tab._engine_combo.setCurrentText('gnina')
        tab._flexible_residues_input.setText('A:123')
        self.results._apply_redocking_btn.click()
        self.assertEqual(tab._flexible_residues_input.text(), 'A:123')
        with patch.object(QMessageBox, 'question', return_value=QMessageBox.No) as confirm, \
                patch('bfee_docking.tabs.docking_tab.docking.Docking') as create:
            tab.start_docking()
            confirm.assert_called_once()
            create.assert_not_called()
        tab._flexible_residues_input.setText('A:124')
        with patch.object(QMessageBox, 'question') as confirm, \
                patch('bfee_docking.tabs.docking_tab.docking.Docking') as create:
            tab.start_docking()
            confirm.assert_not_called()
            self.assertEqual(create.call_args.kwargs['flexible_residues'], 'A:124')

    def test_remote_job_generation_keeps_previous_results_even_on_failure(self):
        for engine, method in [('gnina', 'generate_gnina_job_package'), ('DSDP', 'generate_dsdp_job_package')]:
            with self.subTest(engine=engine):
                self.docking_tab._engine_combo.setCurrentText(engine)
                with patch('bfee_docking.tabs.docking_tab.docking.Docking') as create:
                    generate = getattr(create.return_value, method)
                    generate.return_value = self.output / 'job'
                    self.docking_tab.start_docking()
                    generate.assert_called_once()
                    self.assertIs(self.window._docking_instance, self.docking)
                    generate.side_effect = OSError('cannot write package')
                    with patch('bfee_docking.tabs.docking_tab.traceback.print_exc'):
                        self.docking_tab.start_docking()
                    self.assertIs(self.window._docking_instance, self.docking)
                    self.assertEqual(self.results._pose_combo.currentData(), 0)

    def test_import_keeps_previous_results_until_a_real_pose_is_available(self):
        index = self.output / 'index.csv'
        index.write_text(
            'Filename,OriginalPath,SMILES\n'
            'job.pdbqt,C:\\old\\ligand.pdbqt,CC\n'
            'other.pdbqt,/other/ligand.pdbqt,CCC\n',
            encoding='ascii',
        )
        (self.output / 'job_log.txt').write_text('job.pdbqt -9.21\n', encoding='ascii')
        outcome = self.output / 'job_out.pdbqt'
        outcome.write_text('MODEL 1\nENDMDL\n', encoding='ascii')
        flex = self.output / 'job_out_flex_res.pdbqt'
        atom = f'ATOM      1  C   LIG A   1    {1.0:8.3f}{2.0:8.3f}{3.0:8.3f}\n'
        flex.write_text(atom, encoding='ascii')
        candidate = object.__new__(Docking)
        candidate._ligand_pdbqts = self.window._ligand_pdbqt_files
        candidate._ligand_smiles = self.window._ligand_smiles_list
        self.docking_tab._engine_combo.setCurrentText('DSDP')
        with patch.object(QFileDialog, 'getExistingDirectory', return_value=str(self.output)), \
                patch('bfee_docking.tabs.docking_tab.docking.Docking', return_value=candidate):
            self.docking_tab.load_results()
            self.assertIs(self.window._docking_instance, self.docking)
            self.assertIn('No usable results', self.window._show_warning.call_args.args[1])
            self.assertNotIn('flexible_residue_output', candidate._docking_results[0])
            outcome.write_text(atom, encoding='ascii')
            self.docking_tab.load_results()
        self.assertIs(self.window._docking_instance, candidate)
        self.assertEqual(candidate.get_last_docking_type(), 'flexible')
        self.assertEqual(candidate._docking_results[0]['flexible_residue_output'], flex)
        self.assertEqual(self.results._pose_combo.currentText(), '0 (-9.21 kcal/mol)')
        warning = self.window._show_warning.call_args.args[1]
        self.assertIn('No matching ligand', warning)
        self.assertIn(f'Full logs: {self.output}', warning)

        # Different paths cannot disambiguate identical name/SMILES pairs after a move.
        with index.open('a', encoding='ascii') as manifest:
            manifest.write('duplicate.pdbqt,/another/ligand.pdbqt,CC\n')
        results = candidate.load_dsdp_results(self.output)
        self.assertFalse(results[0]['success'])
        self.assertIn('Multiple matching ligands', results[0]['error'])

    def test_pose_scores_do_not_include_numeric_diagnostics_or_shift_missing_modes(self):
        docking = object.__new__(Docking)
        log = ['mode | affinity | rmsd\n', '----+----------\n', '1 -9.21 0\n']
        self.assertEqual(docking._parse_vina_log(log + ['RUNTIME 12\n', '2 -8.7 2\n']), [-9.21])
        self.assertEqual(docking._parse_vina_log(log + ['2 nan 1\n', '3 -8.7 2\n']), [-9.21])
        self.assertEqual(docking._parse_vina_log(log + ['2 -8.7 2\n']), [-9.21, -8.7])

    def test_running_local_job_blocks_loading_even_after_engine_change(self):
        self.window._docking_worker = Mock()
        self.docking_tab._engine_combo.setCurrentText('gnina')
        self.assertFalse(self.docking_tab._load_results_btn.isEnabled())
        with patch.object(QFileDialog, 'getExistingDirectory') as choose:
            self.docking_tab.load_results()
            choose.assert_not_called()
        self.docking_tab._on_local_worker_stopped()
        self.assertTrue(self.docking_tab._load_results_btn.isEnabled())

    def test_empty_results_clear_previous_pose(self):
        self.window._docking_instance = None
        self.results.update_results_table()
        self.assertEqual(self.results._results_table.rowCount(), 0)
        self.assertEqual(self.results._pose_combo.count(), 0)
        self.assertFalse(self.results._sort_results_btn.isEnabled())
        self.assertFalse(self.results._apply_redocking_btn.isEnabled())
        self.assertEqual(self.results._pose_hint.text(), "")
        self.assertTrue(self.results._pose_hint.isHidden())
