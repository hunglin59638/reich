#!/usr/bin/env python3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from modules.server import app, API_PORT

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=API_PORT, debug=True)
