"""Builds B1v5_collect_zip.ipynb: mounts Drive, packs the results folder skala-b1-v5 into one zip and starts the download.
    python build_collect_zip.py [--out B1v5_collect_zip.ipynb]"""
import argparse, nbformat as nbf

CELL1 = '''
# ---- Cell 1: Google Drive (click Allow in the popup) and find the results folder ----
import os, sys, time
ROOT_NAME = "skala-b1-v5"
ROOT = os.environ.get("B1V5_ROOT")            # (the environment variable is only for tests)
if not ROOT:
    from google.colab import drive
    print("=" * 80); print("CLICK 'Connect to Google Drive' and then Allow in the popup NOW."); print("=" * 80)
    for attempt in range(3):
        try:
            drive.mount("/content/drive", force_remount=(attempt > 0))
            if os.path.isdir("/content/drive/MyDrive"):
                break
        except Exception as ex:
            print(f"Drive mount attempt {attempt + 1} of 3 failed ({repr(ex)[:80]}); trying again in 20 s: please click Allow.")
            time.sleep(20)
    ROOT = "/content/drive/MyDrive/" + ROOT_NAME
print("Results folder:", ROOT, "| exists:", os.path.isdir(ROOT))
if not os.path.isdir(ROOT):
    print("NOT FOUND. Check the folder name in Google Drive (My Drive > skala-b1-v5) and that the right Google account is connected.")
'''

CELL2 = '''
# ---- Cell 2: what is in the folder (nothing is changed) ----
SKIP_DIRS = {"cache", "__pycache__", ".ipynb_checkpoints"}      # cache = Skala weights and Hugging Face files (large, can be downloaded again)
SKIP_SUFFIX = (".tmp",)
files, skipped, total = [], [], 0
for dp, dn, fn in os.walk(ROOT):
    rel = os.path.relpath(dp, ROOT)
    if any(part in SKIP_DIRS for part in rel.split(os.sep)):
        skipped.append(rel); dn[:] = []; continue
    for f in sorted(fn):
        p = os.path.join(dp, f)
        if f.endswith(SKIP_SUFFIX):
            continue
        files.append(p); total += os.path.getsize(p)
print(f"{len(files)} files, {total / 1e6:.1f} MB (left out: {', '.join(skipped) or 'nothing'})")
top = {}
for p in files:
    k = os.path.relpath(p, ROOT).split(os.sep)[0]
    top[k] = top.get(k, 0) + os.path.getsize(p)
for k, v in sorted(top.items()):
    print(f"  {k:45s} {v / 1e6:9.2f} MB")
key = ["RUN_SUMMARY.md", "check_report.md", "check_report.csv", "phase11_dft.csv", "phase12_basis.csv", "phase13_scout.csv", "phase_errors.csv",
       "phase03_ccsdt.csv", "phase08_mp2_basis.csv", "queue_not_run.csv", "findings_outcome.csv", "numbers_v5.csv", "summary_tables.xlsx"]
print("\\nKey files present:")
for k in key:
    print(("  yes  " if os.path.exists(os.path.join(ROOT, k)) else "  NO   ") + k)
sm = os.path.join(ROOT, "RUN_SUMMARY.md")
if os.path.exists(sm):
    print("\\n--- RUN_SUMMARY.md (first lines) ---")
    print("".join(open(sm).readlines()[:40]))
'''

CELL3 = '''
# ---- Cell 3: make the zip and start the download ----
import zipfile, datetime
stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
ZIP = os.path.join("/content" if os.path.isdir("/content") else os.getcwd(), f"skala-b1-v5_{stamp}.zip")
with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in files:
        z.write(p, os.path.join(ROOT_NAME, os.path.relpath(p, ROOT)))
print(f"Zip written: {ZIP} ({os.path.getsize(ZIP) / 1e6:.1f} MB, {len(files)} files)")
try:
    from google.colab import files as colab_files
    colab_files.download(ZIP)
    print("The download should start now (allow downloads if your browser asks). If it does not, open the Files panel on the left, right-click the zip and choose Download.")
except Exception as ex:
    print("Automatic download did not start (", repr(ex)[:80], "): open the Files panel on the left (folder icon), right-click", os.path.basename(ZIP), "and choose Download.")
print("\\nWhen the file is on your computer: Runtime > Disconnect and delete runtime (this notebook does not do it, because that could cancel the download).")
'''

def build(out):
    nb = nbf.v4.new_notebook()
    nb.metadata["colab"] = {"provenance": []}
    nb.cells = [
        nbf.v4.new_markdown_cell("# B1 v5: pack the results folder into one zip and download it\n\nRun on any runtime (a CPU runtime is enough). **Run all**, click **Allow** in the Google Drive popup, wait about a minute, and the zip downloads. "
                                 "It packs `My Drive/skala-b1-v5` (all CSV files, reports, figures, geometries) and leaves out the large `cache` folder (Skala weights). Nothing in Drive is changed or deleted."),
        nbf.v4.new_code_cell(CELL1.strip()),
        nbf.v4.new_code_cell(CELL2.strip()),
        nbf.v4.new_code_cell(CELL3.strip()),
    ]
    nbf.write(nb, out)
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="B1v5_collect_zip.ipynb"); a = ap.parse_args()
    print("wrote", build(a.out))
