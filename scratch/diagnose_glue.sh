#!/bin/bash
# scratch/diagnose_glue.sh
# Run this on the server 'tartarus' inside the true_runner directory

set -euo pipefail

echo "=== Diagnostic Run for GLUE on Server ==="
echo "Current directory: $(pwd)"
echo "Physical path: $(realpath .)"
echo "User: $(whoami)"
echo "Docker group check: $(groups)"

# Check permissions of path components
echo "=== Path Permissions ==="
path_acc=""
for part in $(echo "$(realpath .)" | tr '/' ' '); do
    path_acc="${path_acc}/${part}"
    ls -ld "$path_acc"
done

# Create test directory and files
echo "=== Creating test mount files ==="
mkdir -p results/hbv
TEST_FILE="results/hbv/test_mount.txt"
echo "hello from host physical path" > "$TEST_FILE"
echo "Test file exists on host: $(ls -l "$TEST_FILE")"

# Test Docker volume mount accessibility with Alpine
echo "=== Testing Docker Volume Mount (Alpine) ==="
echo "Running Docker alpine to cat test file..."
set +e
docker run --rm -v "$(realpath results/hbv):/work" alpine cat /work/test_mount.txt
DOCKER_EXIT=$?
set -e

if [ $DOCKER_EXIT -eq 0 ]; then
    echo "✅ Success: Docker (Alpine) can read files mounted from $(realpath results/hbv)"
else
    echo "❌ Error: Docker (Alpine) failed to read mount. Exit code: $DOCKER_EXIT"
fi

# Test Docker volume mount accessibility with GLUE tools image
echo "=== Testing Docker Volume Mount (GLUE tools image) ==="
echo "Running GLUE tools image to list files in /work..."
set +e
docker run --rm -v "$(realpath results/hbv):/work" cvrbioinformatics/gluetools:latest ls -la /work
GLUE_EXIT=$?
set -e

if [ $GLUE_EXIT -eq 0 ]; then
    echo "✅ Success: GLUE tools container can list /work"
else
    echo "❌ Error: GLUE tools container failed. Exit code: $GLUE_EXIT"
fi

# Clean up test file
rm -f "$TEST_FILE"

echo "=== Diagnostic Completed ==="
