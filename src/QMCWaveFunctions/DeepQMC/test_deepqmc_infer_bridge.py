"""Focused tests for DeepQMC Hydra configuration discovery."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from deepqmc_infer_bridge import DeepQMCInferBridge, _find_hydra_config


class HydraConfigDiscoveryTest(unittest.TestCase):
    def test_missing_config_fails_before_checkpoint_loading(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            checkpoint = Path(temp_dir) / "missing-checkpoint.pt"
            with patch.dict(os.environ, {}, clear=True), patch.object(Path, "exists", autospec=True, return_value=False):
                with self.assertRaisesRegex(FileNotFoundError, "DEEPQMC_HYDRA_CONFIG"):
                    DeepQMCInferBridge(checkpoint)

    def test_explicit_missing_config_fails_without_falling_back(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            associated = root / ".hydra" / "config.yaml"
            explicit = root / "explicit-missing.yaml"
            with patch.dict(os.environ, {"DEEPQMC_HYDRA_CONFIG": str(explicit)}), patch.object(
                Path, "exists", autospec=True, side_effect=lambda self: self == associated
            ):
                with self.assertRaisesRegex(FileNotFoundError, "DEEPQMC_HYDRA_CONFIG does not exist"):
                    _find_hydra_config(root / "checkpoint.pt")

    def test_explicit_config_wins_over_associated_config(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            explicit = root / "explicit.yaml"
            associated = root / ".hydra" / "config.yaml"
            with patch.dict(os.environ, {"DEEPQMC_HYDRA_CONFIG": str(explicit)}), patch.object(
                Path, "exists", autospec=True, side_effect=lambda self: self in {explicit, associated}
            ):
                self.assertEqual(_find_hydra_config(root / "checkpoint.pt"), explicit)

    def test_associated_config_is_discovered_in_ancestor(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            associated = root / ".hydra" / "config.yaml"
            checkpoint = root / "training" / "checkpoint.pt"
            with patch.dict(os.environ, {}, clear=True), patch.object(
                Path, "exists", autospec=True, side_effect=lambda self: self == associated
            ):
                self.assertEqual(_find_hydra_config(checkpoint), associated)


if __name__ == "__main__":
    unittest.main()
