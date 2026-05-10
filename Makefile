PYTHON ?= python3
CONFIG ?= config/local_paths.yml
FILTERS ?= config/filters.yml
PDBS ?= config/pdbs.yml

.PHONY: setup inventory build phase1 phase2-smoke phase2-mvp phase3 phase4-build phase4-figures phase4 phase5 phase6 phase7 phase8 phase9 phase7-9 app-export app-install app-dev app-build metrics tokenization structure test clean

setup:
	pip install -e .

inventory:
	$(PYTHON) scripts/01_inventory_data.py --config $(CONFIG) --out results/tables/data_inventory.csv

build:
	$(PYTHON) scripts/02_build_dataset.py --config $(CONFIG) --filters $(FILTERS)

phase1:
	PYTHONPATH=src $(PYTHON) scripts/06_build_phase1_dataset.py --config $(CONFIG) --filters $(FILTERS)

phase2-smoke:
	PYTHONPATH=src $(PYTHON) scripts/07_compute_phase2_metrics.py --panel smoke

phase2-mvp:
	PYTHONPATH=src $(PYTHON) scripts/07_compute_phase2_metrics.py --panel mvp

phase3:
	PYTHONPATH=src $(PYTHON) scripts/08_refine_cds_phase3.py --filters $(FILTERS)

phase4-build:
	PYTHONPATH=src $(PYTHON) scripts/09_build_phase4_representations.py

phase4-figures:
	PYTHONPATH=src $(PYTHON) scripts/10_make_phase4_figures.py

phase4: phase4-build phase4-figures

phase5:
	PYTHONPATH=src $(PYTHON) scripts/11_run_phase5_tokenization_audit.py

phase6:
	PYTHONPATH=src $(PYTHON) scripts/12_run_phase6_tokenization_stability.py

phase7:
	PYTHONPATH=src $(PYTHON) scripts/13_build_phase7_antigenlm_bridge.py

phase8:
	PYTHONPATH=src $(PYTHON) scripts/14_build_phase8_latent_atlas.py

phase9:
	PYTHONPATH=src $(PYTHON) scripts/15_build_phase9_structure_mapping.py

phase7-9: phase7 phase8 phase9 app-export

app-export:
	$(PYTHON) data_export/export_vercel_safe_bundle.py

app-install:
	cd app && npm install

app-dev:
	cd app && npm run dev

app-build:
	cd app && npm run typecheck && npm run build

metrics:
	$(PYTHON) scripts/03_compute_sequence_metrics.py --fasta data/processed/ha_na_sequences.fasta --metadata data/processed/ha_na_metadata.parquet

tokenization:
	$(PYTHON) scripts/04_run_tokenization_audit.py --fasta data/processed/ha_na_sequences.fasta --out results/tables/tokenization_metrics.csv

structure:
	$(PYTHON) scripts/05_make_structure_demo.py --pdb-config $(PDBS) --metrics results/tables/sequence_metrics.csv

test:
	PYTHONPATH=src pytest -q

clean:
	rm -rf __pycache__ .pytest_cache
