import os
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from PySide6.QtWidgets import QApplication, QFileDialog
from rdkit import Chem
from rdkit.Chem import AllChem

from bfee_docking.md_prepper import MDPrepper, third_party_tools
from bfee_docking.tabs.md_prepper_tab import MDPrepperTab


class LigandPrepTests(unittest.TestCase):
    def setUp(self):
        temp = TemporaryDirectory(prefix="bfee ligand test ")
        self.addCleanup(temp.cleanup)
        self.output = Path(temp.name)
        self.source = self.output / "md_prepper" / "ligand.pdb"
        self.source.parent.mkdir()
        self.source.write_text("test pose\n", encoding="ascii")
        self.fixed = self.source.parent / "fixed_files"
        self.stream = self.fixed / "ligand.str"
        self.base = self.output / "force fields"
        self.base.mkdir()
        for filename in ("top_all36_cgenff.rtf", "par_all36_cgenff.prm", "par_all36_cgenff_increments.prm"):
            (self.base / filename).write_text(f"* CGenFF Force Field v. 5.0\n{filename}\n", encoding="ascii")
        self.prepper = object.__new__(MDPrepper)
        self.prepper._output_dir = str(self.output)
        self.prepper._ligand_pdb_path = str(self.source)
        self.prepper._ligand_smiles = "CCO"

    def stub_mol2(self, pdb_path, smiles, output):
        self.assertEqual(pdb_path.read_bytes(), self.source.read_bytes())
        self.assertEqual(smiles, "CCO")
        output.write_text("mol2\n", encoding="ascii")
        return self.molecule

    def stub_stream(self, molecule, force_field, output):
        self.assertIs(molecule, self.molecule)
        for file in self.base.iterdir():
            self.assertEqual((force_field / file.name).read_bytes(), file.read_bytes())
        output.write_text("* AutoCGENFF local-v5 ligand stream\nRESI UNL 0\n", encoding="ascii")

    @staticmethod
    def stub_ionize(**kwargs):
        for key in ("output_psf_path", "output_pdb_path"):
            Path(kwargs[key]).write_text("ionized\n", encoding="ascii")

    def test_process_ligand_and_preserve_matching_base_files(self):
        self.molecule = object()
        with patch.object(self.prepper, "_write_ligand_mol2", side_effect=self.stub_mol2), \
             patch.object(self.prepper, "_write_ligand_stream", side_effect=self.stub_stream), \
             patch.object(self.prepper, "_write_autocgenff_pdb", side_effect=lambda _, p: p.touch()):
            self.assertEqual(self.prepper.process_ligand(str(self.base)), str(self.stream))
            self.assertEqual(self.prepper.process_ligand(str(self.base)), str(self.stream))
        self.assertTrue(self.stream.is_file())
        self.assertEqual(self.stream.with_suffix(".pdb").read_bytes(), self.source.read_bytes())
        self.assertTrue(self.stream.with_suffix(".mol2").is_file())
        saved = self.prepper.get_ligand_force_field_dir(self.stream, "a different directory")
        self.assertEqual(saved, self.fixed / "cgenff_v5")
        self.assertEqual((saved / "par_all36_cgenff.prm").read_bytes(),
                         (self.base / "par_all36_cgenff.prm").read_bytes())

    def test_assignment_failure_keeps_pdb_and_mol2_without_reusing_old_stream(self):
        self.molecule = object()
        self.fixed.mkdir()
        self.stream.write_text("* AutoCGENFF local-v5 ligand stream\n", encoding="ascii")
        with patch.object(self.prepper, "_write_ligand_mol2", side_effect=self.stub_mol2), \
             patch.object(self.prepper, "_write_ligand_stream", side_effect=ValueError("unsupported ligand")):
            with self.assertRaisesRegex(ValueError, "unsupported ligand"):
                self.prepper.process_ligand(str(self.base))
        self.assertFalse(self.stream.exists())
        self.assertEqual(self.stream.with_suffix(".pdb").read_bytes(), self.source.read_bytes())
        self.assertTrue(self.stream.with_suffix(".mol2").is_file())

    def test_does_not_overwrite_a_manually_supplied_stream(self):
        self.molecule = object()
        self.fixed.mkdir()
        self.stream.write_text("custom topology\n", encoding="ascii")
        with patch.object(self.prepper, "_write_ligand_mol2", side_effect=self.stub_mol2):
            with self.assertRaisesRegex(ValueError, "Move or rename your existing file"):
                self.prepper.process_ligand(str(self.base))
        self.assertEqual(self.stream.read_text(encoding="ascii"), "custom topology\n")

    def test_psf_uses_current_stream_and_corresponding_base(self):
        manual = self.output / "custom.top"
        manual.write_text("* Custom\n*\n36 1\nRESI MYL 0.0\nEND\n", encoding="ascii")
        self.fixed.mkdir()
        self.stream.write_text("* AutoCGENFF local-v5 ligand stream\nRESI UNL 0\n", encoding="ascii")
        temp_dir = self.output / "temp_files"
        md_dir = self.output / "md_files"
        self.prepper.generate_psf = Mock(return_value=("protein.psf", "protein.pdb"))
        self.prepper.generate_ligand_psf = Mock(return_value=("ligand.psf", "ligand.pdb"))
        self.prepper.merge_psf_pdb = Mock(return_value=("complex.psf", "complex.pdb"))
        self.prepper.solvate_system = Mock(return_value=("solvated.psf", "solvated.pdb"))
        self.prepper.ionize_system = Mock(side_effect=self.stub_ionize)

        for selected, base in ((self.stream, self.base), (manual, third_party_tools.get_force_field_dir())):
            with self.subTest(selected=selected.name):
                self.prepper.generate_system_psf(
                    "protein.pdb", str(self.source), str(selected), 15.0, "NaCl", 0.1,
                    str(temp_dir), str(md_dir), cgenff_dir=str(self.base),
                )
                topologies = self.prepper.generate_ligand_psf.call_args.kwargs["topology_files"]
                self.assertEqual(Path(topologies[0]), base / "top_all36_cgenff.rtf")
                self.assertEqual((md_dir / "par_all36_cgenff.prm").read_bytes(),
                                 (base / "par_all36_cgenff.prm").read_bytes())
                exported = (md_dir / "ligand.str").read_text(encoding="ascii")
                self.assertIn("RESI UNL", exported)
                if selected == manual:
                    self.assertIn("read rtf card append", exported)
                    self.assertIn("RESI MYL", manual.read_text(encoding="ascii"))

    def test_manual_v5_stream_requires_matching_base(self):
        stream = self.output / "download.str"
        stream.write_text("* For use with CGenFF version 5.0\nRESI LIG 0\n", encoding="ascii")
        with self.assertRaisesRegex(FileNotFoundError, "Select the CGenFF 5.0 directory"):
            self.prepper.get_ligand_force_field_dir(stream)
        self.assertEqual(self.prepper.get_ligand_force_field_dir(stream, str(self.base)), self.base)
        (self.base / "top_all36_cgenff.rtf").write_text("* CGenFF Force Field v. 4.6\n", encoding="ascii")
        with self.assertRaisesRegex(ValueError, "Expected a CGenFF 5.0 data file"):
            self.prepper.get_ligand_force_field_dir(stream, str(self.base))

    def test_incompatible_pose_is_not_mapped_by_index(self):
        molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertEqual(AllChem.EmbedMolecule(molecule, randomSeed=19), 0)
        Chem.MolToPDBFile(molecule, str(self.source))
        template = Chem.AddHs(Chem.MolFromSmiles("COC"))
        with self.assertRaisesRegex(ValueError, "connectivity does not match"):
            self.prepper._assign_bond_orders_from_template(template, self.source)

    def test_vmd_rejects_uncaught_script_errors_and_stale_outputs(self):
        path = self.output / "ligand.psf"
        path.write_text("previous run", encoding="ascii")
        self.prepper._vmd_executable = self.output / "vmd"
        args = self.output / "test.tcl", self.output, path
        with patch.object(third_party_tools, "run_vmd", return_value="Tcl error"):
            with self.assertRaisesRegex(RuntimeError, "did not complete"):
                self.prepper._run_vmd_outputs(*args)
        with patch.object(third_party_tools, "run_vmd", return_value="BFEE_DOCKING_VMD_DONE\n"):
            with self.assertRaisesRegex(RuntimeError, "did not create a fresh"):
                self.prepper._run_vmd_outputs(*args)
        self.assertEqual(path.read_text(encoding="ascii"), "previous run")

    @unittest.skipUnless(os.environ.get("AUTOCGENFF_V5_DIR") and os.environ.get("BFEE_TEST_VMD"),
                         "set AUTOCGENFF_V5_DIR and BFEE_TEST_VMD for full PSF assembly")
    def test_full_system_retains_ligand_lone_pairs(self):
        protein = self.output / "protein [test].pdb"
        molecule = Chem.AddHs(Chem.MolFromSequence("AAA"))
        self.assertEqual(AllChem.EmbedMolecule(molecule, randomSeed=19), 0)
        Chem.MolToPDBFile(Chem.RemoveHs(molecule), str(protein))
        self.prepper = MDPrepper(str(protein), str(self.source), "", "Clc1ccccc1", 7.0,
                                 str(self.output), os.environ["BFEE_TEST_VMD"])
        molecule = Chem.AddHs(Chem.MolFromSmiles(self.prepper._ligand_smiles))
        self.assertEqual(AllChem.EmbedMolecule(molecule, randomSeed=19), 0)
        Chem.MolToPDBFile(Chem.RemoveHs(molecule), str(self.source))
        stream = Path(self.prepper.process_ligand(os.environ["AUTOCGENFF_V5_DIR"]))
        pose = stream.with_name("ligand_autocgenff.pdb")
        self.assertEqual(sum(line.startswith("HETATM") for line in pose.read_text().splitlines()), 12)
        psf, _ = self.prepper.generate_system_psf(
            str(protein), str(stream.with_suffix(".pdb")), str(stream), 12.0, "NaCl", 0.1,
            str(self.output / "temp [system]"), str(self.output / "md [system]"),
        )
        lone_pairs = next(line for line in Path(psf).read_text().splitlines() if "!NUMLP" in line)
        self.assertEqual(lone_pairs.split()[0], "1")
        if python := os.environ.get("BFEE_TEST_OPENMM_PYTHON"):
            subprocess.run([python, "-c", """
from openmm.app import CharmmParameterSet, CharmmPsfFile, NoCutoff
from openmm import Context, Platform, VerletIntegrator
from openmm.app import PDBFile
from openmm.unit import kilojoule_per_mole
import math
params = CharmmParameterSet('par_all36m_prot.prm', 'par_all36_cgenff.prm', 'toppar_water_ions.str', 'ligand.str')
psf = CharmmPsfFile('ionized.psf')
system = psf.createSystem(params, nonbondedMethod=NoCutoff)
assert sum(system.isVirtualSite(i) for i in range(system.getNumParticles())) == 1
context = Context(system, VerletIntegrator(0.001), Platform.getPlatformByName('Reference'))
context.setPositions(PDBFile('ionized.pdb').positions)
context.computeVirtualSites()
energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole)
assert math.isfinite(energy), energy
"""], cwd=Path(psf).parent, check=True, timeout=60)

    @unittest.skipUnless(os.environ.get("AUTOCGENFF_V5_DIR"), "set AUTOCGENFF_V5_DIR for real assignment")
    def test_real_assignment_and_optional_psfgen(self):
        molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertEqual(AllChem.EmbedMolecule(molecule, randomSeed=19), 0)
        # Exercise a non-UNL docked residue and a pose missing nonpolar hydrogens.
        pose = Chem.MolToPDBBlock(Chem.RemoveHs(molecule)).replace(" UNL ", " MOL ")
        self.source.write_text(pose, encoding="ascii")
        stream = Path(self.prepper.process_ligand(os.environ["AUTOCGENFF_V5_DIR"]))
        names = [line.split()[1] for line in stream.read_text().splitlines() if line.startswith("ATOM ")]
        self.assertEqual(len(names), 9)
        self.assertTrue({"C1", "C2", "O1"}.issubset(names))
        self.assertEqual(stream.with_suffix(".pdb").read_bytes(), self.source.read_bytes())
        if vmd := os.environ.get("BFEE_TEST_VMD"):
            self.prepper._vmd_executable = Path(vmd)
            base = self.prepper.get_ligand_force_field_dir(stream)
            psf, pdb = self.prepper.generate_ligand_psf(
                str(stream.with_suffix(".pdb")), [str(base / "top_all36_cgenff.rtf"), str(stream)]
            )
            self.assertIn("9 !NATOM", Path(psf).read_text(encoding="ascii"))
            self.assertTrue(Path(pdb).is_file())

            manual = self.output / "manual.top"
            topology = stream.read_text().split("read rtf card append\n", 1)[1].split("\nEND", 1)[0]
            manual.write_text(topology.replace("RESI UNL", "RESI MAN") + "\nEND\n", encoding="ascii")
            # Exercise the app's custom-topology handling and real ligand psfgen.
            # The rest of the system pipeline does not need to run for this test.
            self.prepper.generate_psf = Mock(return_value=("protein.psf", "protein.pdb"))
            self.prepper.merge_psf_pdb = Mock(return_value=("complex.psf", "complex.pdb"))
            self.prepper.solvate_system = Mock(return_value=("solvated.psf", "solvated.pdb"))
            self.prepper.ionize_system = Mock(side_effect=self.stub_ionize)
            temp_dir = self.output / "temp system"
            self.prepper.generate_system_psf(
                "protein.pdb", str(stream.with_suffix(".pdb")), str(manual), 15.0, "NaCl", 0.1,
                str(temp_dir), str(self.output / "md system"),
            )
            self.assertIn("9 !NATOM", (temp_dir / "ligand.psf").read_text(encoding="ascii"))
            self.assertIn("RESI MAN", manual.read_text(encoding="ascii"))


class LigandPrepGuiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        cls.app = QApplication.instance() or QApplication([])

    def test_autofill_override_and_failure(self):
        with TemporaryDirectory() as temp:
            output = Path(temp)
            working = output / "md_prepper"
            working.mkdir()
            (working / "ligand.pdb").touch()
            (working / "ligand.smi").write_text("CCO", encoding="ascii")
            fixed_dir = working / "fixed_files"
            fixed_dir.mkdir()
            (fixed_dir / "protein_final.pdb").touch()
            (fixed_dir / "ligand.pdb").touch()
            (fixed_dir / "ligand_input.smi").write_text("CCO", encoding="ascii")
            main = SimpleNamespace(
                _output_directory=temp, _vmd_path="test-vmd", _cgenff_dir=temp, _last_dir=temp,
                _show_info=Mock(), _show_warning=Mock(), _show_error=Mock(),
            )
            tab = MDPrepperTab(main)
            self.addCleanup(tab.close)
            auto = str(working / "fixed_files" / "ligand.str")
            manual = str(output / "manual.top")
            with patch("bfee_docking.tabs.md_prepper_tab.md_prepper.MDPrepper") as factory:
                backend = factory.return_value
                backend.process_ligand.return_value = auto
                backend.generate_system_psf.return_value = ("ionized.psf", "ionized.pdb")
                tab._md_process_ligand_btn.click()
                self.assertEqual(tab._ligand_str_edit.text(), auto)
                backend.process_ligand.assert_called_once_with(cgenff_dir=temp)

                with patch.object(QFileDialog, "getOpenFileName", return_value=(manual, "")) as browse:
                    tab._browse_ligand_str_file()
                self.assertIn("*.top", browse.call_args.args[-1])
                self.assertEqual(tab._ligand_str_edit.text(), manual)
                tab._ligand_str_edit.setText(str(output / "typed.str"))
                tab._generate_psf_btn.click()
                self.assertEqual(backend.generate_system_psf.call_args.kwargs["ligand_str_path"],
                                 str(output / "typed.str"))
                main._show_error.assert_not_called()

                (working / "ligand.smi").write_text("COC", encoding="ascii")
                tab._generate_psf_btn.click()
                backend.generate_system_psf.assert_called_once()
                self.assertIn("exported ligand has changed", main._show_warning.call_args.args[-1])

                backend.process_ligand.side_effect = ValueError("unsupported ligand")
                with patch("bfee_docking.tabs.md_prepper_tab.traceback.print_exc"):
                    tab._md_process_ligand_btn.click()
                main._show_error.assert_called_once()
                self.assertEqual(tab._ligand_str_edit.text(), str(output / "typed.str"))


if __name__ == "__main__":
    unittest.main()
