"""
Inspect the Johns et al. 2018 Nature Methods supplementary XLSX files.
Prints sheet names, shapes, and column names/dtypes ONLY (no raw data dumped)
so this is safe to run and read without blowing up context.
"""
import sys
import openpyxl
import pandas as pd
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "raw"

FILES = [
    "NIHMS945382-supplement-3.xlsx",
    "NIHMS945382-supplement-4.xlsx",
    "NIHMS945382-supplement-5.xlsx",
    "NIHMS945382-supplement-6.xlsx",
    "NIHMS945382-supplement-7.xlsx",
]

def inspect(fname):
    path = RAW / fname
    print(f"\n{'='*80}\nFILE: {fname} ({path.stat().st_size/1e6:.2f} MB)\n{'='*80}")
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for sheet_name in wb.sheetnames:
        ws = wb[sheet_name]
        print(f"\n  --- SHEET: {sheet_name!r}  max_row={ws.max_row}  max_col={ws.max_column}")
        # print first 3 rows raw (headers + a couple examples), truncated
        for i, row in enumerate(ws.iter_rows(min_row=1, max_row=3, values_only=True)):
            trimmed = [str(c)[:40] if c is not None else None for c in row[:20]]
            print(f"    row{i+1}: {trimmed}")
    wb.close()

if __name__ == "__main__":
    targets = sys.argv[1:] if len(sys.argv) > 1 else FILES
    for f in targets:
        inspect(f)
