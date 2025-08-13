#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Usage: $0 ./cesar"
  exit 1
fi

set -e

HERE=$(dirname $0)
cesar=$1

# check if cesar successfully returns the test samples.

echo "Checking for correct behaviour."
for f in $(find $HERE/extra/samples -name '*.in'); do
  diff -q <($cesar $f) $HERE/extra/samples/$(basename $f .in).out
  if [ $? -ne 0 ]; then
    echo "$f... NOK"
    exit 1
  else
    echo "$f... OK"
  fi
done
exit 0
