#!/bin/bash
find . -size +50M | sed 's|^./||g' | grep -v ".git" | cat >> .gitignore


