# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import sys

from monty.serialization import loadfn


def main():
    lines = [loadfn(filename).__str__() for filename in sys.argv[1:]]
    print("\n\n".join(lines))


if __name__ == "__main__":
    main()
