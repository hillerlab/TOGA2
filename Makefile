## TOGA2 supports bash only; make's default recipe shell is /bin/sh, which is dash on
## Debian/Ubuntu. Set it explicitly so recipes and check_shell agree.
SHELL := /bin/bash

ARCH=$(shell uname -m)
CHECK_DEPS=check_dependencies.py
DELIM="=============================="
MAIN=toga2.py
EXEC_SCRIPTS=cesar_exec.py cesar_preprocess.py classify_chains.py feature_extractor.py fine_orthology_resolver.py get_contig_sizes.py predict_with_spliceai.py train_model.py
EXEC_MODULES=chain_bst_index.py get_names_from_bed.py
.PHONY: all build build_cesar build_c build_cython build_rust pull_submodules check check_essentials check_managers check_python check_shell check_third_party chmod install_python_packages install_binaries install_postoga train_models

all: build

check: check_shell check_essentials check_managers

build: chmod install check build_c build_cesar build_cython build_rust train_models

install: install_binaries install_python install_third_party install_postoga

build_c:
	if [ ARCH = "arm64" ]; then \
		CFLAGS="-Wall -Wextra -O2 -g -std=c99 -arch arm64"; \
	else \
		CFLAGS="-Wall -Wextra -O2 -g -std=c99"; \
	fi; \
	mkdir -p bin && gcc $$CFLAGS -o bin/chain_filter_by_id src/c/chain_filter_by_id.c
	gcc $$CFLAGS -fPIC -shared -o src/python/modules/chain_coords_converter_slib.so src/c/chain_coords_converter_slib.c
	gcc $$CFLAGS -fPIC -shared -o src/python/modules/chain_bst_lib.so src/c/chain_bst_lib.c
	echo ${DELIM}

build_cesar:
	cd CESAR2.0 && make
	echo ${DELIM}

build_cython:
	mkdir -p bin && python3 setup.py build_ext --build-lib=src
	echo ${DELIM}

build_rust:
	cd src/rust && cargo build --release
	echo ${DELIM}
	cd bed2gtf && cargo build --release

pull_submodules:
	git submodule update --init --recursive

check_essentials:
	python3 ${CHECK_DEPS} essentials
	echo ${DELIM}

check_managers:
	python3 ${CHECK_DEPS} managers
	echo ${DELIM}

check_python:
	python3 ${CHECK_DEPS} python --installation_mode

check_shell:
	@if [ -n "$${BASH_VERSION:-}" ]; then \
		echo "bash detected ($$BASH_VERSION); shell type check successfully passed"; \
	else \
		echo "ERROR: TOGA2 currently only supports bash as operating shell. Please change your default shell or create a separate environment with bash as default" >&2; \
		exit 1; \
	fi

check_third_party:
	./${CHECK_DEPS} third_party

chmod:
	chmod +x ${CHECK_DEPS} ; \
	chmod +x ${MAIN}; \
	for exec in ${EXEC_SCRIPTS}; do \
		chmod +x src/python/$${exec}; \
	done ; \
	for exec in ${EXEC_MODULES}; do \
		chmod +x src/python/modules/$${exec}; \
	done

install_binaries:
	mkdir -p bin && \
	wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/twoBitToFa && \
	wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/faToTwoBit && \
	wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/wigToBigWig && \
#	wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bigWigToWig && \
	wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bigBedToBed && \
	wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bedToBigBed && \
	wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/ixIxx && \
		chmod +x bin/*
	echo ${DELIM}

install_python_packages:
	python3 -m pip install -r assets/requirements/requirements.txt; \
	echo ${DELIM}

install_python:
	python3 -m pip install -r assets/requirements/requirements.txt

install_third_party:
	python3 ${CHECK_DEPS} install_third_party

install_postoga:
	cd postoga/rustools && \
	maturin develop --release

train_models:
	python3 src/python/train_model.py
	echo ${DELIM}
