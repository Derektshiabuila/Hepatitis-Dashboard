import subprocess
import os
import sys
from pathlib import Path

def run_diagnostics():
    print("=== Running Headless Display Diagnostics inside Docker ===")
    
    uid = os.getuid() if hasattr(os, "getuid") else 0
    gid = os.getgid() if hasattr(os, "getgid") else 0
    print(f"Host User: {uid}:{gid}")

    image_name = "local-wine-vb6-xvfb:v3"
    
    # We will execute a diagnostic shell script inside the container
    shell_script = """set -x
echo "=== Container Environment ==="
id
env

echo "=== Directory Permissions ==="
ls -la /tmp
ls -la /tmp/.X11-unix || echo "No .X11-unix"

echo "=== Test 1: Start Xvfb directly (Unix socket) ==="
Xvfb :99 -screen 0 640x480x8 > /tmp/xvfb_unix.log 2>&1 &
XVFB_PID=$!
sleep 2

echo "Check if Xvfb process is running:"
ps -p $XVFB_PID || echo "Xvfb is NOT running!"
cat /tmp/xvfb_unix.log

echo "Check socket directory:"
ls -la /tmp/.X11-unix

echo "Try connecting with xsec/xauth/xlsclients:"
DISPLAY=:99 xlsclients || echo "Failed to connect to :99 via Unix socket"

kill $XVFB_PID
wait $XVFB_PID 2>/dev/null

echo "=== Test 2: Start Xvfb with TCP listening ==="
# Modern Xvfb disables TCP by default. We enable it with -listen tcp
Xvfb :99 -screen 0 640x480x8 -listen tcp > /tmp/xvfb_tcp.log 2>&1 &
XVFB_PID=$!
sleep 2

echo "Check if Xvfb process is running:"
ps -p $XVFB_PID || echo "Xvfb is NOT running!"
cat /tmp/xvfb_tcp.log

echo "Try connecting via TCP (127.0.0.1:99):"
DISPLAY=127.0.0.1:99 xlsclients || echo "Failed to connect to 127.0.0.1:99 via TCP"

kill $XVFB_PID
wait $XVFB_PID 2>/dev/null
"""
    
    cmd = [
        "docker", "run", "--rm",
        "-e", "HOME=/tmp",
        "-u", f"{uid}:{gid}",
        image_name,
        "sh", "-c", shell_script
    ]
    
    print("Running docker command...")
    res = subprocess.run(cmd, capture_output=False)
    print(f"Docker finished with code: {res.returncode}")

if __name__ == "__main__":
    run_diagnostics()
