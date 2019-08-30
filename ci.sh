#!/bin/bash

black -l 88 --check .
flake8 --max-line-length=88 .
pytest -s
