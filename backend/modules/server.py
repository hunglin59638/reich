#!/usr/bin/env python3
from pathlib import Path
import subprocess
from flask import Flask, request, jsonify

from modules.common import read_config

API_PORT = read_config()["api_port"]
app = Flask(__name__)


@app.route("/version", methods=["GET"])
def version():
    return jsonify(read_config()["version"])


@app.route("/run", methods=["GET"])
def run():
    project_id = request.args.get("project_id")
    return jsonify({"project_id": project_id})
