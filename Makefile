PYTHON ?= python3
VENV ?= .venv
PIP := $(VENV)/bin/python -m pip
RUN := $(VENV)/bin/python

.PHONY: venv install smoke smoke-perm scan-smoke compare-smoke baseline mps clean

venv:
	$(PYTHON) -m venv $(VENV)
	$(RUN) -m ensurepip --upgrade

install: venv
	$(PIP) install -U pip
	$(PIP) install -r requirements.txt

smoke:
	$(RUN) qa_adiabatic_steps_bench.py -n 6 --instances 5 --t-max 3 --shots 32 --aer-method statevector --opt-ref exact

smoke-perm:
	$(RUN) qa_adiabatic_steps_bench.py -n 6 --instances 5 --t-max 3 --shots 32 \
		--aer-method statevector --opt-ref exact --stats-method perm --perm-iterations 1000 --no-plots --outdir qa_out_perm

scan-smoke:
	$(RUN) qa_adiabatic_steps_bench.py --n-list 4,5 --instances 3 --t-max 2 --shots 32 \
		--aer-method statevector --opt-ref exact --stats-method mw --no-plots --outdir qa_scan_smoke

compare-smoke:
	$(RUN) compare_qa_sa_prr.py -n 4 --instances 1 --t-max 1 --delta-t 0.5 --shots 16 \
		--aer-method statevector \
		--sa-reads 16 --sa-sweep-checkpoints 8,16 \
		--prr-total-time 1 --prr-segments 4 --prr-maxiter-bfgs 4 --prr-maxiter-nm 8 \
		--no-plots --outdir diagnostics_local/compare_smoke

baseline:
	$(RUN) qa_adiabatic_steps_bench.py -n 6 --instances 100 --t-max 10 --shots 128 --aer-method statevector --opt-ref exact

mps:
	$(RUN) qa_adiabatic_steps_bench.py -n 30 --instances 50 --t-max 10 --shots 64 \
		--aer-method matrix_product_state \
		--mps-max-bond-dimension 128 \
		--mps-truncation-threshold 1e-10

clean:
	rm -rf qa_out qa_out_perm qa_scan_smoke
