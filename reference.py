"""Pinned poseidon-tools dependency. Install with: python3 reference.py"""

from hashlib import sha256
from io import BytesIO
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from urllib.request import urlopen
from zipfile import ZipFile


REVISION = "60075da7c0521d9493749a035b1f30d4eda37138"
ARCHIVE_SHA256 = "05321e138871ce63eeb247583ac827905d5ff570eb198cae79d5f593fe09ebed"
URL = f"https://codeload.github.com/khovratovich/poseidon-tools/zip/{REVISION}"
DEPENDENCY = Path(__file__).resolve().parent / ".deps" / f"poseidon-tools-{REVISION}"


def install():
    """Download the unmodified source once; subsequent runs work offline."""
    if DEPENDENCY.is_dir():
        print(f"poseidon-tools {REVISION} is already installed.")
        return
    with urlopen(URL, timeout=60) as response:
        archive = response.read()
    if sha256(archive).hexdigest() != ARCHIVE_SHA256:
        raise RuntimeError("poseidon-tools archive checksum mismatch.")
    DEPENDENCY.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=DEPENDENCY.parent) as temporary:
        with ZipFile(BytesIO(archive)) as source:
            source.extractall(temporary)
        (Path(temporary) / DEPENDENCY.name).rename(DEPENDENCY)
    print(f"Installed poseidon-tools {REVISION} in {DEPENDENCY}")


if __name__ == "__main__":
    install()
else:
    if not (DEPENDENCY / "poseidon" / "poseidon.py").is_file():
        raise ImportError("Missing poseidon-tools dependency. "
                          "Run 'python3 reference.py' from the Slipway directory.")
    sys.path.insert(0, str(DEPENDENCY))
    from poseidon.poseidon import Poseidon
