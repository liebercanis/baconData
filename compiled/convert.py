#!/usr/bin/env python3
import sys
import os
from subprocess import Popen, PIPE
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import pprint


def main(args):
    """ run convertRaw """
    dirs = []
    p = os.listdir('data')
    print(p)
    for i in p:
        if os.path.isdir('data/'+i):
            dirs.append(i)
            print("\n is dir ", i)
        else:
            print("\n not dir ", i)

    n = len(dirs)
    if (n < 1):
        return

    print(" files %i ", len(p), " dirs %i ", len(dirs))
    if (len(sys.argv) > 1):
        n = int(args[0])

    print(" args ", args, " number of directories to convert  ", n)
    for i in range(0, n):
        print(" run job %i %d", i, dirs[i])
        process = Popen(['test', dirs[i]], stdout=PIPE, stderr=PIPE)
        process.wait()
        stdout, stderr = process.communicate()
        print(stdout)


if __name__ == '__main__':
    main(sys.argv[1:])
