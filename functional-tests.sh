#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest
set -uo pipefail

go build -o allhic_test cmd/main.go

run test_optimize ./allhic_test optimize tests/simulation/test.ids tests/simulation/test.clm
assert_in_stderr "Success"

run test_optimize_skipga ./allhic_test optimize tests/simulation/test.ids tests/simulation/test.clm --skipGA
assert_in_stderr "Success"
