import hashlib
import io
import logging
import subprocess
import sys
from typing import List
import urllib
from urllib.request import Request, urlopen


def run(args: List[str], show_output: bool = True) -> subprocess.CompletedProcess:
    """Run a command line.

    Parameters
    ----------
    args: List[str]
        A list of argument
    show_output: bool (default: True)
        Output command line

    Return
    ------
    subprocess.CompletedProcess
        Return result obtained with subprocess
    """
    ret = subprocess.run(args, capture_output=True, encoding="utf8")
    if show_output and ret.stdout is not None:
        logging.info(ret.stdout)
    if show_output and ret.stderr is not None:
        logging.warning(ret.stderr)
    return ret


def md5(path: str) -> str:
    hash_md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def url_download(url: str, path: str) -> None:
    try:
        with urlopen(Request(url)) as fod:
            with open(path, "wb") as dst:
                while True:
                    chunk = fod.read(2**10)
                    if chunk:
                        dst.write(chunk)
                    else:
                        break
    except Exception as e:
        print(str(e))

def url_download_to_memory(url: str) -> io.BytesIO:
    count = 0
    error_code = 0
    while count < 3:
        try:
            memory_buffer = io.BytesIO()
            with urlopen(Request(url)) as fod:
                while True:
                    chunk = fod.read(2**10)
                    if chunk:
                        memory_buffer.write(chunk)
                    else:
                        break
            memory_buffer.seek(0)  # Move the cursor to the beginning of the buffer
            return memory_buffer, error_code
        except urllib.error.HTTPError as e:
            if e.code > 501:
                count += 1
            else:
                return None, e.code
            error_code = e.code
    return None, error_code
