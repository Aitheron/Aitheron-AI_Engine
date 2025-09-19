from pathlib import Path

def is_ready_file(path) -> bool:
    p = Path(path)
    return p.is_file() and p.stat().st_size > 0
