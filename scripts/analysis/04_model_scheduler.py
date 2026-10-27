"""Run independent R workers within a fourteen-core CPU budget.

Inputs: ENDO_WORK_DIR and the ENDO_MANIFEST, ENDO_SCRIPT and ENDO_BATCH settings.
Outputs: worker logs and a batch status record under the private workspace.

Four workers use at most 4 + 4 + 4 + 2 simultaneous Stan chains. Each Bayesian
model retains four chains, including models run by the two-core worker.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import csv
import json
import os
import subprocess

ROOT = Path(os.environ["ENDO_WORK_DIR"]).expanduser().resolve()
CODE_DIR = Path(__file__).resolve().parent
manifest = os.environ.get("ENDO_MANIFEST", "calendar_new_models_manifest.csv")
script = os.environ.get("ENDO_SCRIPT", "04_run_calendar_models.R")
batch = os.environ.get("ENDO_BATCH", "calendar")

with (ROOT / "08_qa" / manifest).open() as handle:
    names = [row["id"] for row in csv.DictReader(handle)]


def run_group(index):
    group = names[index::4]
    if not group:
        return dict(worker=index + 1, endpoints=[], returncode=0)

    env = dict(
        os.environ,
        ENDO_MODEL_WORKER="1",
        ENDO_CHAIN_CORES=str([4, 4, 4, 2][index]),
    )
    with (ROOT / f"10_logs/{batch}_worker_{index + 1}.log").open("w") as log:
        result = subprocess.run(
            ["Rscript", str(CODE_DIR / script), *group],
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return dict(worker=index + 1, endpoints=group, returncode=result.returncode)


with ThreadPoolExecutor(max_workers=4) as pool:
    results = list(pool.map(run_group, range(4)))

(ROOT / f"10_logs/{batch}_batch_results.json").write_text(
    json.dumps(results, indent=2)
)
raise SystemExit(1 if any(x["returncode"] != 0 for x in results) else 0)
