import os
import sys
import subprocess
from pathlib import Path


def main():
    args = sys.argv
    cmd = args

    exec_dir = None
    cmd[0] = '/usr/bin/git'
    for i, opt in enumerate(cmd):
        if opt == '-C':
            exec_dir = cmd[i+1]
            del cmd[i:i+2]
            break

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=exec_dir)
    while True:
        line = proc.stdout.readline()
        sys.stdout.write(line.decode(encoding='utf-8'))

        if not line and proc.poll() is not None:
            break


if __name__ == '__main__':
    main()
