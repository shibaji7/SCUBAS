#!/usr/bin/env python

import argparse
import os
import shutil


def clean():
    if os.path.exists("dist/"):
        shutil.rmtree("dist/")
    if os.path.exists("build/"):
        shutil.rmtree("build/")
    if os.path.exists("scubas.egg-info/"):
        shutil.rmtree("scubas.egg-info/")
    os.system("find . -type d -name '.ipynb_checkpoints' -exec rm -rf {} +")
    os.system("find . -type d -name '__pycache__' -exec rm -rf {} +")
    return


def build():
    os.system("isort -rc -sl .")
    os.system("autoflake --remove-all-unused-imports -i -r .")
    os.system("isort -rc -m 3 .")
    os.system("black .")
    os.system("python setup.py sdist bdist_wheel")
    return


def uplaod_pip():
    clean()
    build()
    os.system("twine upload dist/* --verbose")
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-b", "--build", action="store_true", help="Build the project from scratch"
    )
    parser.add_argument(
        "-rm", "--clean", action="store_true", help="Cleaning pre-existing files"
    )
    parser.add_argument(
        "-rb",
        "--rebuild",
        action="store_true",
        help="Clean and rebuild the project from scratch",
    )
    parser.add_argument(
        "-upip",
        "--upip",
        action="store_true",
        help="Upload code to PIP repository",
    )
    args = parser.parse_args()
    if args.clean:
        clean()
    if args.build:
        build()
    if args.rebuild:
        clean()
        build()
    if args.upip:
        uplaod_pip()
        
