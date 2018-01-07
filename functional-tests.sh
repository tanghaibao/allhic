#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest
set -uo pipefail

go build -o main_test main.go

run test_optimize ./main_test optimize tests/test.clm
assert_in_stderr "Success"

run test_optimize_skipga ./main_test optimize tests/test.clm --skipGA
assert_in_stderr "Success"
