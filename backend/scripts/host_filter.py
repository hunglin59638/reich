#!/usr/bin/env python

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from modules.host_filter import cli

if __name__ == "__main__":
    cli()
