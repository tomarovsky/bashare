#!/bin/bash

set -e -u

while getopts 'i:' flag; do
        case "${flag}" in
                i) username="${OPTARG}" ;;
                *) exit 1 ;;
        esac
done

# script
getent passwd | grep ${username}


