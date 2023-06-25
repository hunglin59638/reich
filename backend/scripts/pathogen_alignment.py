#!/usr/env/bin python
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from backend.alignment import cli

if __name__ == "__main__":
    import click
    from pathlib import Path

    @click.command()
    # @click.argument("input-file", type=Path, required=True, nargs=-1)
    @click.argument("reference", type=Path, required=True, nargs=1)
    @click.argument("query", type=Path, required=True, nargs=-1)
    def main(query):
        print(query)

    main()
