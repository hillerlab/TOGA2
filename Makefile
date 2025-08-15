ARCH=$(shell uname -m)
CHECK_DEPS=check_dependencies.py
DELIM="=============================="
MAIN=toga2.py
EXEC_SCRIPTS=cesar_exec.py cesar_preprocess.py classify_chains.py feature_extractor.py train_model.py
EXEC_MODULES=chain_bst_index.py get_contig_sizes.py get_names_from_bed.py
VENV ?= false
VENV_NAME ?= toga2
.PHONY: all build build_cesar build_c build_cython build_rust check check_essentials check_managers check_python check_shell check_third_party chmod install_python_packages install_binaries train_models

all: build

check: check_shell check_essentials check_python check_managers check_third_party

build: chmod check install build_c build_cesar build_cython build_rust train_models

install: install_binaries install_python_packages install_third_party

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

check_essentials:
	./${CHECK_DEPS} essentials
	echo ${DELIM}

check_managers:
	./${CHECK_DEPS} managers
	echo ${DELIM}

check_python:
	./${CHECK_DEPS} python --installation_mode

check_shell:
	if [ $(echo $0) != "-bash" ]; then \
		echo "ERROR: TOGA2 currently only supports bash as operating shell. Please change you default shell or create a separate environment with bash as default"; \
	else \
		echo "bash has been found to be a default current shell; shell type check successfully passed" ; \
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
	wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bigWigToWig && \
	wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bigBedToBed && \
	wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/bedToBigBed && \
	wget -P bin/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v479/ixIxx && \
        chmod +x bin/*
	echo ${DELIM}

install_python_packages:
	if [[ ! -f missing_packages.txt ]] && [[ ${VENV} == false ]]; then \
		echo "No missing Python packages in the current environment; exiting"; \
	else \
		if [[ ${VENV} == true ]]; then \
			echo "Creating a separate Python environment called ${VENV_NAME}"; \
			python3 -m venv ${VENV_NAME} && \
			source ${VENV_NAME}/bin/activate && \
			python3 -m pip install -r requirements.txt; \
		else \
			echo "Installing missing packages globally"; \
                        cat missing_packages.txt ; \
			python3 -m pip install -r missing_packages.txt --root-user-action=ignore; \
		fi; \
	fi ; \
	echo ${DELIM}

install_third_party:
	if [[ ${VENV} == true ]]; then \
		source ${VENV_NAME}/bin/activate ; \
	fi; \
	./${CHECK_DEPS} install_third_party

train_models:
	src/python/train_model.py
	echo ${DELIM}
