"""Installation orchestration without downloading or modifying real environments."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from test_multi_coverage import MODULE as cli

SOURCE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("ebp_install", SOURCE / "scripts/install.py")
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)


class InstallerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.prefix = Path(self.tmp.name).resolve() / "ebp"
        self.commands = []

    def fake_run(self, command):
        self.commands.append(command)
        if command[1] in {"create", "install"}:
            target = Path(command[command.index("--prefix") + 1])
            (target / "bin").mkdir(parents=True, exist_ok=True)
            (target / "bin/python").touch()
            (target / "conda-meta").mkdir(exist_ok=True)
            (target / "conda-meta/history").write_text("test")

    def test_install_configures_all_tools_without_shell_exports(self):
        with patch.dict(os.environ, {}, clear=False), patch.object(installer, "run", self.fake_run), \
                patch.object(installer.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "Usage", "")):
            installer.install(self.prefix, "/mock/conda")
        package = self.prefix / "share/eukbinpicker"
        config = json.loads((package / "installation.json").read_text())
        self.assertEqual(config["version"], "0.41")
        self.assertFalse((self.prefix.with_name("." + self.prefix.name + ".ebp-install") / "settings.json").exists())
        self.assertFalse((package / "install-cache").exists())
        self.assertTrue((self.prefix / "bin/eukbinpicker").is_file())
        for tool in ("concoct", "metabat2", "eukcc", "busco"):
            self.assertEqual(config["environments"][tool.upper() + "_ENV"],
                             str(package / "environments" / tool))
        flattened = " ".join(" ".join(x) for x in self.commands)
        for requirement in ("vamb==5.0.4", "remag==0.4.4", "samtools>=1.17", "busco=6"):
            self.assertIn(requirement, flattened)
        self.assertIn("https://download.pytorch.org/whl/cu124", flattened)
        with patch.object(cli, "ROOT", package), patch.dict(os.environ, {}, clear=True):
            self.assertEqual(cli.environment_value("EUKCC_ENV"), str(package / "environments/eukcc"))
            with patch.dict(os.environ, {"EUKBINPICKER_EUKCC_ENV": "explicit"}):
                self.assertEqual(cli.environment_value("EUKCC_ENV"), "explicit")

    def test_failed_probe_leaves_no_runnable_cli_and_can_resume(self):
        failed = subprocess.CompletedProcess([], 127, "", "missing shared library")
        with patch.dict(os.environ, {}, clear=False), patch.object(installer, "run", self.fake_run), \
                patch.object(installer.subprocess, "run", return_value=failed):
            with self.assertRaisesRegex(RuntimeError, "validation failed"):
                installer.install(self.prefix, "/mock/conda")
        self.assertFalse((self.prefix / "bin/eukbinpicker").exists())
        self.assertTrue((self.prefix.with_name("." + self.prefix.name + ".ebp-install") / "settings.json").is_file())
        self.commands.clear()
        with patch.dict(os.environ, {}, clear=False), patch.object(installer, "run", self.fake_run), \
                patch.object(installer.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "", "")):
            installer.install(self.prefix, "/mock/conda")
        self.assertFalse(any("--prefix" in command for command in self.commands))
        self.assertTrue((self.prefix / "bin/eukbinpicker").is_file())

    def test_existing_environment_is_not_overwritten(self):
        self.prefix.mkdir()
        sentinel = self.prefix / "keep"
        sentinel.write_text("unchanged")
        with self.assertRaisesRegex(ValueError, "new installation directory"):
            installer.install(self.prefix, "/mock/conda")
        self.assertEqual(sentinel.read_text(), "unchanged")

    def test_resume_repairs_package_records_left_without_python(self):
        failed = subprocess.CompletedProcess([], 127, "", "probe failed")
        with patch.dict(os.environ, {}, clear=False), patch.object(installer, "run", self.fake_run), \
                patch.object(installer.subprocess, "run", return_value=failed):
            with self.assertRaises(RuntimeError):
                installer.install(self.prefix, "/mock/conda")
        (self.prefix / "bin/python").unlink()
        (self.prefix / "conda-meta/python.json").write_text(json.dumps(
            {"name": "python", "version": "3.12.16", "build": "test_0"}))
        self.commands.clear()
        with patch.dict(os.environ, {}, clear=False), patch.object(installer, "run", self.fake_run), \
                patch.object(installer.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "", "")):
            installer.install(self.prefix, "/mock/conda")
        self.assertIn("--force-reinstall", self.commands[0])
        self.assertIn("python=3.12.16=test_0", self.commands[0])
        self.assertLess(self.commands[0].index("--force-reinstall"),
                        self.commands[0].index("python=3.12"))

    def test_readme_documents_every_public_option(self):
        readme = (SOURCE / "README.md").read_text()
        for action in cli.build_parser()._actions:
            for option in action.option_strings:
                self.assertIn(f"`{option}`", readme)
