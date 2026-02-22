"""R Script Runner - executes R scripts as subprocesses for bioinformatics computation."""

import json
import logging
import os
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Directory where R scripts live (relative to this file: ../../r_scripts)
R_SCRIPTS_DIR = Path(__file__).parent.parent / "r_scripts"


@dataclass
class RScriptResult:
    returncode: int
    stdout: str
    stderr: str
    elapsed_seconds: float

    @property
    def success(self) -> bool:
        return self.returncode == 0


class RScriptError(Exception):
    """Raised when an R script exits with non-zero return code."""

    def __init__(self, script: str, result: RScriptResult):
        self.script = script
        self.result = result
        super().__init__(
            f"R script '{script}' failed with return code {result.returncode}.\n"
            f"STDERR:\n{result.stderr[-3000:] if len(result.stderr) > 3000 else result.stderr}"
        )


class RScriptRunner:
    """Manages R script execution and Python-R communication via config JSON files."""

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.config_path = self.output_dir / "config.json"

    def check_r_available(self) -> bool:
        """Check if Rscript is available on PATH."""
        try:
            result = subprocess.run(
                ["Rscript", "--version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    def write_config(self, config: dict) -> None:
        """Atomically write config JSON for R scripts to read.

        Uses a temp file + os.replace for atomic write to avoid partial reads.
        """
        tmp_fd, tmp_path = tempfile.mkstemp(
            dir=self.output_dir, suffix=".json.tmp", prefix="config_"
        )
        try:
            with os.fdopen(tmp_fd, "w", encoding="utf-8") as f:
                json.dump(config, f, indent=2, ensure_ascii=False)
            os.replace(tmp_path, self.config_path)
            logger.debug(f"Config written to {self.config_path}")
        except Exception:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
            raise

    def read_config(self) -> dict:
        """Read current config JSON."""
        with open(self.config_path, "r", encoding="utf-8") as f:
            return json.load(f)

    def run(
        self,
        script_name: str,
        timeout: int = 3600,
        extra_args: Optional[list] = None,
    ) -> RScriptResult:
        """Run an R script from r_scripts/ directory.

        Args:
            script_name: Filename of the R script (e.g., '01_data_prep.R')
            timeout: Maximum seconds to wait (default 3600 = 1 hour)
            extra_args: Optional additional CLI arguments passed to Rscript

        Returns:
            RScriptResult with returncode, stdout, stderr, elapsed_seconds

        Raises:
            FileNotFoundError: If the R script does not exist
            RScriptError: If the R script exits with non-zero return code
        """
        script_path = R_SCRIPTS_DIR / script_name
        if not script_path.exists():
            raise FileNotFoundError(f"R script not found: {script_path}")

        cmd = ["Rscript", "--vanilla", str(script_path), str(self.config_path)]
        if extra_args:
            cmd.extend(extra_args)

        logger.info(f"Running R script: {script_name}")
        start_time = time.time()

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            raise RScriptError(
                script_name,
                RScriptResult(
                    returncode=-1,
                    stdout="",
                    stderr=f"R script timed out after {timeout} seconds",
                    elapsed_seconds=time.time() - start_time,
                ),
            )

        elapsed = time.time() - start_time
        result = RScriptResult(
            returncode=proc.returncode,
            stdout=proc.stdout,
            stderr=proc.stderr,
            elapsed_seconds=elapsed,
        )

        if proc.stdout:
            logger.debug(f"R stdout:\n{proc.stdout[-2000:]}")
        if proc.stderr:
            logger.debug(f"R stderr:\n{proc.stderr[-2000:]}")

        logger.info(f"R script {script_name} finished in {elapsed:.1f}s (rc={proc.returncode})")

        if not result.success:
            raise RScriptError(script_name, result)

        return result

    def validate_outputs(self, expected_files: list[str]) -> list[str]:
        """Check that expected output files exist.

        Args:
            expected_files: List of file paths relative to output_dir, or absolute paths

        Returns:
            List of missing file paths (empty if all present)
        """
        missing = []
        for f in expected_files:
            path = Path(f) if Path(f).is_absolute() else self.output_dir / f
            if not path.exists():
                missing.append(str(path))
        return missing
