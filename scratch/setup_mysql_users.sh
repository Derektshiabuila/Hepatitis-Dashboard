#!/bin/bash
# scratch/setup_mysql_users.sh
# Run this on the server 'tartarus' to configure the MySQL containers

set -euo pipefail

for virus in hbv hcv hev; do
    container="gluetools-mysql-${virus}"
    echo "=== Configuring database for ${container} ==="
    
    # Check if container exists
    if ! docker ps -a --format '{{.Names}}' | grep -q "^${container}$"; then
        echo "Warning: Container ${container} does not exist. Skipping."
        continue
    fi

    # Check if container is running
    if ! docker ps --format '{{.Names}}' | grep -q "^${container}$"; then
        echo "Container ${container} is not running. Starting it..."
        docker start "${container}"
        sleep 5
    fi
    
    # Run MySQL commands to create GLUE_TOOLS database and gluetools user
    echo "Creating database GLUE_TOOLS and user 'gluetools'..."
    docker exec "${container}" mysql -uroot -e "
        CREATE DATABASE IF NOT EXISTS GLUE_TOOLS;
        CREATE USER IF NOT EXISTS 'gluetools'@'%' IDENTIFIED BY 'glue12345';
        ALTER USER 'gluetools'@'%' IDENTIFIED BY 'glue12345';
        GRANT ALL PRIVILEGES ON GLUE_TOOLS.* TO 'gluetools'@'%';
        FLUSH PRIVILEGES;
    " || {
        echo "Failed to configure database for ${container}."
        continue
    }
    
    # Verify configuration
    echo "Verifying user list for ${container}:"
    docker exec "${container}" mysql -uroot -e "SELECT user, host FROM mysql.user;"
    echo "Verifying databases for ${container}:"
    docker exec "${container}" mysql -uroot -e "SHOW DATABASES;"
    echo "------------------------------------------------"
done

echo "=== Setup Completed ==="
