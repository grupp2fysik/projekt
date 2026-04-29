#!/bin/bash

# Remove built file from project
find . -type f -name "*.out" -delete
find . -type f -name "*.xml" -delete
find . -type f -name "slurm-*.out" -delete
find . -type d -name "*.save" -exec rm -rf {} +

echo "Clean complete."
