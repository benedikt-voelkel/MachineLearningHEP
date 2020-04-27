#!/bin/bash -e

TEST=$1
#set -o pipefail
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

ERR=0


function swallow() {
    local ERR=0
    local TMPF=$(mktemp /tmp/swallow.XXXX)
    local MSG=$1
    shift
    printf "[    ] $MSG" >&2
    "$@" &> $TMPF || ERR=$?
    if [[ $ERR != 0 ]]; then
        printf "\r[\033[31mFAIL\033[m] $MSG (log follows)\n" >&2
        cat $TMPF
        printf "\n" >&2
    else
        printf "\r[ \033[32mOK\033[m ] $MSG\n" >&2
    fi
    rm -f $TMPF
    return $ERR
}


function test-pylint()
{
    local err=0
    local test_files=$@
    echo "run test: pylint"
    type pylint
    for tf in $test_files; do
        echo "File $tf "
        swallow "linting $tf" pylint $tf || ERR=1
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

exit $ERR
