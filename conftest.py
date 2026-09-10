"""Test configuration: add src/ to sys.path so tests can import methylation modules."""
import sys
from pathlib import Path

# Add the parent src/ directory to sys.path
_src_path = Path(__file__).resolve().parent.parent / "src"
if str(_src_path) not in sys.path:
    sys.path.insert(0, str(_src_path))
