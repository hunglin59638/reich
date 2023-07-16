#!/usr/bin/env python3

import os
import gzip
import click
import json
from pathlib import Path
from logging import getLogger, StreamHandler, FileHandler, DEBUG, WARN, INFO, Formatter


logger = getLogger(__file__)
APP_ROOT = Path(__file__).parent.parent
MAX_THREADS = os.cpu_count()

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


def set_threads(func):
    def wrapper(*args, **kwargs):
        threads = kwargs.get("threads", 0)
        kwargs["threads"] = str(MAX_THREADS) if int(threads) == 0 else str(threads)
        return func(*args, **kwargs)

    return wrapper


def set_out_dir(func):
    def wrapper(*args, **kwargs):
        out_dir = kwargs.get("out_dir", None)
        out_dir = Path().cwd() if out_dir is None or out_dir == "." else Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        kwargs["out_dir"] = out_dir
        return func(*args, **kwargs)

    return wrapper


def get_basename(path):
    path = Path(path)
    return (
        path.with_suffix("").with_suffix("")
        if path.name.endswith(".gz")
        else path.with_suffix("")
    )


def set_binaries_path():
    bin_dir = os.path.join(APP_ROOT, "Reich/bin")
    logger.info(f"Adding {bin_dir} to PATH.")
    os.environ["PATH"] = bin_dir + ":" + os.environ["PATH"]


def open_gz(file, mode=None):
    GZIP_MAGIC_NUMBER = b"\x1f\x8b"

    def _is_gz_file():
        if mode is None or "r" in mode:
            with open(file, "rb") as f:
                return f.read(2) == GZIP_MAGIC_NUMBER
        else:
            True

    mode = mode if mode else "rt"
    if _is_gz_file():
        return gzip.open(file, mode)
    else:
        return open(file, mode)


def check_seq_format(file):
    with open_gz(file) as f:
        header = f.read(1)
        if header == "@":
            return "fastq"
        elif header == ">":
            return "fasta"
        else:
            return "unknown"


def generate_batch(file, out_dir=None):
    # deprecated
    file = Path(file)
    if not out_dir:
        out_dir = Path().cwd()
    else:
        out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bps = 0
    n = 0
    BATCH_COUNT = 1000000000
    seq_fmt = check_seq_format(file)
    batch_seq = ""
    if seq_fmt == "unknown":
        raise TypeError(f"{file} is not fasta or fastq")
    for rec in SeqIO.parse(open_gz(file), seq_fmt):
        bps += len(rec)
        batch_seq += rec.format(seq_fmt)
        if bps < BATCH_COUNT:
            pass
        else:
            batch_file = f"{file}_batch_{str(n)}"
            with open(batch_file, "w") as f:
                f.write(batch_seq)
                batch_seq = ""
            bps = 0
            n += 1
            yield batch_file
    # finally
    if batch_seq:
        n += 1
        batch_file = f"{file}_batch_{str(n)}"
        with open(batch_file, "w") as f:
            f.write(batch_seq)
        return batch_file


class AdvancedJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, Path):
            return str(obj.resolve())
        return json.JSONEncoder.default(self, obj)


class BaseNameType(click.ParamType):
    name = "basename"

    def convert(self, value, param, ctx):
        try:
            if isinstance(value, str):
                value = Path(value)
                value = value.parent.resolve() / value.name
                if len([i for i in value.parent.glob(f"{value.name}.*")]) > 0:
                    return value
                else:
                    self.fail(f"{value} is not a basename path", param, ctx)
            elif value is None:
                return value
            else:
                self.fail(f"{value} is not a basename path", param, ctx)
        except ValueError:
            self.fail(f"{value} is not a basename path", param, ctx)
