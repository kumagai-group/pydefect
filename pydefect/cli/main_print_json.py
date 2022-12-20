# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import sys

from monty.serialization import loadfn


def main():
    for filename in sys.argv[1:]:
        print("-"*80)
        print(f"file: {filename}")
        print(loadfn(filename))


if __name__ == "__main__":
    main()
