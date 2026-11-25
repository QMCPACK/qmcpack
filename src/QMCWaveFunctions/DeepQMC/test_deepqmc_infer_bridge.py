"""Focused tests for DeepQMC Hydra configuration discovery."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from deepqmc_infer_bridge import DeepQMCInferBridge, _find_hydra_config


class HydraConfigDiscoveryTest(unittest.TestCase):
    def test_missing_adjacent_config_fails_before_checkpoint_loading(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            checkpoint = Path(temp_dir) / "missing-checkpoint.pt"
            with patch("deepqmc_infer_bridge._load_checkpoint_without_h5py") as load_checkpoint:
                with self.assertRaisesRegex(FileNotFoundError, "beside the checkpoint"):
                    DeepQMCInferBridge(checkpoint)
                load_checkpoint.assert_not_called()

    def test_explicit_valid_config_is_used(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            explicit = root / "explicit.yaml"
            explicit.touch()
            associated = root / ".hydra" / "config.yaml"
            associated.parent.mkdir()
            associated.touch()
            self.assertEqual(_find_hydra_config(root / "checkpoint.pt", explicit), explicit)

    def test_explicit_missing_config_does_not_fall_back(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            associated = root / ".hydra" / "config.yaml"
            associated.parent.mkdir()
            associated.touch()
            with self.assertRaisesRegex(FileNotFoundError, "config path does not exist"):
                _find_hydra_config(root / "checkpoint.pt", root / "missing.yaml")

    def test_explicit_relative_path_uses_runtime_working_directory(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            checkpoint_dir = root / "training"
            checkpoint_dir.mkdir()
            cwd_dir = root / "runtime"
            cwd_dir.mkdir()
            (cwd_dir / "config.yaml").touch()
            original_cwd = Path.cwd()
            try:
                os.chdir(cwd_dir)
                self.assertEqual(
                    _find_hydra_config(checkpoint_dir / "checkpoint.pt", "config.yaml"),
                    Path("config.yaml"),
                )
            finally:
                os.chdir(original_cwd)

    def test_associated_config_is_discovered_next_to_checkpoint(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            checkpoint_dir = Path(temp_dir) / "training"
            associated = checkpoint_dir / ".hydra" / "config.yaml"
            associated.parent.mkdir(parents=True)
            associated.touch()
            self.assertEqual(_find_hydra_config(checkpoint_dir / "checkpoint.pt"), associated)

    def test_ancestor_config_is_not_discovered(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            ancestor_config = root / ".hydra" / "config.yaml"
            ancestor_config.parent.mkdir()
            ancestor_config.touch()
            checkpoint = root / "training" / "checkpoint.pt"
            with self.assertRaisesRegex(FileNotFoundError, "beside the checkpoint"):
                _find_hydra_config(checkpoint)

    def test_legacy_environment_variable_is_ignored(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            explicit = root / "legacy.yaml"
            explicit.touch()
            with patch.dict(os.environ, {"DEEPQMC_HYDRA_CONFIG": str(explicit)}):
                with self.assertRaisesRegex(FileNotFoundError, "beside the checkpoint"):
                    _find_hydra_config(root / "checkpoint.pt")

    def test_empty_explicit_path_uses_adjacent_default(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            checkpoint_dir = Path(temp_dir)
            associated = checkpoint_dir / ".hydra" / "config.yaml"
            associated.parent.mkdir()
            associated.touch()
            self.assertEqual(_find_hydra_config(checkpoint_dir / "checkpoint.pt", ""), associated)


if __name__ == "__main__":
    unittest.main()
