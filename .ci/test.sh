#!/bin/bash

set -e

python3 -bb -m pytest -v -x -W always --cov PIL --cov Tests --cov-report term Tests
