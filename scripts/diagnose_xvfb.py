import subprocess
import os
import sys
from pathlib import Path

def run_diagnostics():
    print("=== Running Headless Display Diagnostics inside Docker ===")
    
    uid = os.getuid() if hasattr(os, "getuid") else 0
    gid = os.getgid() if hasattr(os, "getgid") else 0
    print(f"Host User: {uid}:{gid}")

    # Find scripts dir relative to this file
    scripts_dir = Path(__file__).parent.resolve()
    if scripts_dir.name != "scripts":
        scripts_dir = scripts_dir / "scripts"
    
    # In true_runner context, we want true_runner/scripts
    if not (scripts_dir / "RDP5CL.exe").exists():
        scripts_dir = scripts_dir.parent / "true_runner" / "scripts"
        
    print(f"Scripts Directory: {scripts_dir}")

    image_name = "local-wine-vb6-xvfb:v3"
    
    # We will execute a diagnostic shell script inside the container
    shell_script = """set -x
echo "=== Container Environment ==="
id

# 1. Initialize Wine prefix
mkdir -p /tmp/wineprefix
WINEARCH=win32 WINEPREFIX=/tmp/wineprefix wineboot --init

# 2. Copy RDP5 files
mkdir -p /tmp/wineprefix/drive_c/Program Files/RDP5
cp RDP.ini PairsScores BinProbs /tmp/wineprefix/drive_c/Program Files/RDP5/
cp /opt/dlls/MSVBVM60.DLL /tmp/wineprefix/drive_c/windows/system32/

# 3. Start Xvfb directly in background
Xvfb :99 -screen 0 640x480x8 > /tmp/xvfb.log 2>&1 &
XVFB_PID=$!
sleep 2

# Verify Xvfb started
ps -p $XVFB_PID || (echo "Xvfb failed to start" && cat /tmp/xvfb.log && exit 1)

# 4. Run RDP5CL.exe using wine pointing to display :99
export DISPLAY=:99
export WINEPREFIX=/tmp/wineprefix
wine RDP5CL.exe -fAL_7.fasta -nor
EXIT_CODE=$?
echo "Wine exit code: $EXIT_CODE"

kill $XVFB_PID
wait $XVFB_PID 2>/dev/null
exit $EXIT_CODE
"""
    
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{scripts_dir}:/work",
        "-w", "/work",
        "-e", "WINEDEBUG=+err",
        "-e", "WINEDLLOVERRIDES=mscoree,mshtml=d",
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
