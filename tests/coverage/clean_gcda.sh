#!/bin/sh

find . -name \*.gcda | xargs -i rm {}
