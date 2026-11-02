import os, glob, zipfile, tarfile
import pandas as pd

print("=== 03_Troposphere Files ===")
for f in glob.glob("03_Troposphere/*"):
    print(" ", f, os.path.getsize(f), "bytes")
    if f.endswith(".zip"):
        with zipfile.ZipFile(f, "r") as z:
            names = z.namelist()
            print(f"    ZIP contains {len(names)} files. First 3:", names[:3])
            with z.open(names[0]) as fp:
                lines = [fp.readline().decode("utf-8", errors="ignore").strip() for _ in range(25)]
                for ln in lines:
                    if ln.startswith("+TROP/SOLUTION") or "TROTOT" in ln or not ln.startswith("%"):
                        print("    Sample:", ln[:100])

print("\n=== 04_Gas_S5P/multi_gas_s5p.csv Summary ===")
mg = pd.read_csv("04_Gas_S5P/multi_gas_s5p.csv")
print(mg.info())
print(mg.tail(12))