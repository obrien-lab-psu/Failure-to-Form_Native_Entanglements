import os
import pymol

# Manually set the path to your PyMOL license file
os.environ["PYMOL_LICENSE_FILE"] = "/storage/home/ims86/pymol_license_inv58630.lic"
print("License file path set to:", os.environ["PYMOL_LICENSE_FILE"])

# Launch PyMOL
pymol.finish_launching(["pymol", "-qc"])

# Test PyMOL functionality (e.g., a licensed-only feature)
try:
    pymol.cmd.ray(800, 600)  # Ray tracing requires a valid license
    print("PyMOL ray tracing succeeded: License is valid.")
except pymol.CmdException as e:
    print("PyMOL ray tracing failed:", str(e))

