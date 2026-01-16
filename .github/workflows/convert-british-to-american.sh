#!/usr/bin/env bash

set -e

find . -type f -print0 -iname "*.r" -or -iname "*.[qr]md" -or -iname "*.[r]markdown" |
    xargs -0 sed -i \
        -e 's/visualize/visualize/' \
        -e 's/summarize/summarize/' \
        -e 's/analyzed/analyzed/' \
        -e 's/\<analyze\>/analyze/' \
        -e 's/favor/favor/' \
        -e 's/visualisation/visualization/' \
        -e 's/visualisation/visualization/'
