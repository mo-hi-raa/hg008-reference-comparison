"""Print the capstone result at the end of pixi run reproduce."""
import csv

with open("results/capstone/comparison/summary.tsv") as handle:
    row = next(csv.DictReader(handle, delimiter="\t"))

print()
print("CBIO 7501 capstone result")
print(f"Known: {row['known']}  Calls: {row['calls']}  Found: {row['found']}  "
      f"Missed: {row['missed']}  Extra: {row['extra']}")
print(f"Recall: {row['recall']}  Precision: {row['precision']}")
print("Report: results/capstone/report/capstone-report.html")
