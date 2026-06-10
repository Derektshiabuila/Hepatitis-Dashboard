#!/bin/bash
# scratch/restore_mysql_databases.sh
# Run this on the server 'tartarus' to restore the GLUE databases

set -euo pipefail

# Ensure we are in the project folder
cd "$(dirname "$0")/.."

echo "=== Restoring GLUE Databases on Server ==="

for virus in hbv hcv hev; do
    container="gluetools-mysql-${virus}"
    dump_file="${virus}_db.sql.gz"
    
    if [ ! -f "$dump_file" ]; then
        echo "Error: Dump file ${dump_file} not found in root directory!"
        exit 1
    fi
    
    echo "Restoring ${virus} database inside ${container}..."
    
    # Check if container is running
    if ! docker ps --format '{{.Names}}' | grep -q "^${container}$"; then
        echo "Starting ${container}..."
        docker start "${container}"
        sleep 5
    fi
    
    # Ensure GLUE_TOOLS database is clean and recreated
    docker exec "${container}" mysql -uroot -e "DROP DATABASE IF EXISTS GLUE_TOOLS; CREATE DATABASE GLUE_TOOLS;"
    
    # Import dump
    gunzip -c "$dump_file" | docker exec -i "${container}" mysql -uroot GLUE_TOOLS
    
    # Verify tables
    echo "Verification for ${container}:"
    docker exec "${container}" mysql -uroot -e "SHOW TABLES IN GLUE_TOOLS;"
    echo "------------------------------------------------"
done

echo "=== Restore Completed Successfully ==="
