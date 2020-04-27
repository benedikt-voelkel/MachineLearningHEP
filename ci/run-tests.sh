#!/bin/bash -e

TEST=$1
set -o pipefail
cd "$(dirname "$0")"/..

#############################################################################
##  © Copyright CERN 2018. All rights not expressly granted are reserved.  ##
##                 Author: Gian.Michele.Innocenti@cern.ch                  ##
## This program is free software: you can redistribute it and/or modify it ##
##  under the terms of the GNU General Public License as published by the  ##
## Free Software Foundation, either version 3 of the License, or (at your  ##
## option) any later version. This program is distributed in the hope that ##
##  it will be useful, but WITHOUT ANY WARRANTY; without even the implied  ##
##     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    ##
##           See the GNU General Public License for more details.          ##
##    You should have received a copy of the GNU General Public License    ##
##   along with this program. if not, see <https://www.gnu.org/licenses/>. ##
#############################################################################




function test-pylint()
{
    local test_files=$@
    echo "run test: pylint"
    type pylint
    for tf in $test_files; do
        echo "File $tf "
        pylint $tf
    done
}


function test-case()
{
    case "$1" in
        pylint)
            shift
            test-pylint $@
            ;;
        *)
            echo "Unknown test case $1"
            ;;
    esac
}


function test-all()
{
    echo "Run all tests"
    test-pylint $@
}


function print-help()
{
    echo
    echo "run_tests.sh usage to run CI tests"
    echo ""
    echo "run_tests.sh [<testcase>|all]  # defaults to all"
    echo ""
    echo "Possible test cases are:"
    echo "  pylint                        # run style tests for copyright and pylint"
    echo ""
    echo "--help|-h                         # Show this message and exit"
}


[[ $# == 0 ]] && { echo "ERROR: Arguments required" ; print-help ; exit 1; }

FILES=""
TESTS=""
CURRENT_ARG=""

while [[ $# -gt 0 ]]; do
    case "$1" in

        --tests)
            CURRENT_ARG="tests"
            ;;
        --files)
            CURRENT_ARG="files"
            ;;
        --help|-h)
            print-help
            exit 1
            ;;

        *)
            case "$CURRENT_ARG" in
                files)
                    FILES+=" $1"
                    ;;
                tests)
                    TESTS+=" $1"
                    ;;
                *)
                    echo "Unknown option $1"
                    print-help
                    exit 2
                    ;;
            esac
            ;;
    esac
    shift
done


[[ "$FILES" == "" ]] && { echo "ERROR: No files to test. Exit..."; exit 1; }




if [[ "$TESTS" == "" ]]
then
    test-all $FILES
else
    for t in $TESTS
    do
        echo "Do test for $t"
        test-case $t $FILES
    done
fi

