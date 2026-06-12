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

    # Copy hev_RDP.ini from root/parent directory to scripts_dir/RDP.ini
    root_dir = scripts_dir.parent
    src_ini = root_dir / "hev_RDP.ini"
    dest_ini = scripts_dir / "RDP.ini"
    if src_ini.exists():
        import shutil
        print(f"Copying config: {src_ini} -> {dest_ini}")
        shutil.copy(src_ini, dest_ini)
    else:
        print(f"Warning: {src_ini} not found")

    image_name = "local-wine-vb6-xvfb:v3"

    
    # We will execute a diagnostic shell script inside the container
    shell_script = """set -x
echo "=== Container Environment ==="
id

# Verify files exist in work dir
ls -la

# Run the entire sequence inside xvfb-run (without -nolisten unix)
xvfb-run -a --server-args="-screen 0 640x480x8" sh -c "
  set -x
  mkdir -p /tmp/wineprefix
  WINEARCH=win32 WINEPREFIX=/tmp/wineprefix wineboot --init
  mkdir -p '/tmp/wineprefix/drive_c/Program Files/RDP5'
  cp RDP.ini PairsScores BinProbs '/tmp/wineprefix/drive_c/Program Files/RDP5/'
  cp /opt/dlls/MSVBVM60.DLL '/tmp/wineprefix/drive_c/windows/system32/'
  WINEPREFIX=/tmp/wineprefix wine RDP5CL.exe -fAL_7.fasta -nor
"
EXIT_CODE=$?
echo "Execution finished with exit code: $EXIT_CODE"
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
    
    # Clean up RDP.ini
    if dest_ini.exists():
        print(f"Removing temporary config: {dest_ini}")
        os.remove(dest_ini)

if __name__ == "__main__":
    run_diagnostics()
