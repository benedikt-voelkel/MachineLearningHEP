#!/bin/bash -e

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

cd "$(dirname "$0")"/..


function test-dryrun()
{
    echo "run test: dryrun"
    local databases=$(find data -maxdepth 1 -name "database_ml_parameters*")
    for db in $databases
    do
        local has_MBvspt_ntrkl=$(grep "MBvspt_ntrkl:" $db)
        if [[ "$has_MBvspt_ntrkl" != "" ]]
        then
            continue
        fi
        python ci/dryrun.py -r submisison/default_complete.yaml -d $db -a MBvspt_ntrkl
    done
}


function test-pylint()
{
    local test_files="$@"
    echo "run test: pylint"
    type pylint
    for tf in $test_files; do
        pylint $tf
    done
}


function test-all()
{
    echo "Run all tests"
    test-pylint "$1" "$2" "$3"
    test-dryrun
}


function print-help()
{
    echo
    echo "run_tests.sh usage to run CI tests"
    echo ""
    echo "run_tests.sh [<testcase>|all]  # defaults to all"
    echo ""
    echo "Possible test cases are:"
    echo "  style                        # run style tests for copyright and pylint"
    echo "  dryrun                       # run a dryrun test"
    echo ""
    echo "--help|-h                         # Show this message and exit"
}


[[ $# == 0 ]] && test-all


FILES_CREATED=""
FILES_CHANGED=""
FILES_DELETED=""

DO_PYLINT=""
DO_DRYRUN=""


while [[ $# -gt 0 ]]; do
    case "$1" in

        all) test-all ;;
        pylint) DO_PYLINT="1" ;;
        dryrun) DO_DRYRUN="1" ;;

        --files-created)
            shift
            FILES_CREATED="$1"
            ;;
        --files-changed)
            shift
            FILES_CHANGED="$1"
            ;;
        --files-deleted)
            shift
            FILES_DELETED="$1"
            ;;
        --help|-h)
            print-help
            exit 1
            ;;

        *)
            echo "Unknown option $1"
            print-help
            exit 2
            ;;
    esac
    shift
done
