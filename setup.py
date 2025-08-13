#!/usr/bin/env python3

from setuptools import Extension, setup
from Cython.Build import cythonize
from typing import List

import os
import platform

C_FLAGS: List[str] = ['-Wall', '-Wextra', '-O2', '-g', '-std=c99']
if platform.machine() == 'arm64':
    C_FLAGS.extend(['-arch arm64'])
C_FLAGS_EXTENDED: List[str] = [*C_FLAGS, '-fPIC', '-shared']

setup(
    name='TOGA2',
    version='2.0.2',
    python_requires=">=3.6",
    ext_modules = [
        *cythonize(
        [
            Extension(
                # os.path.join('python', 'modules', 'extract_chain_projections'), 
                'python.modules.extract_chain_projections',
                [os.path.join('src', 'cython', 'extract_chain_projections.pyx')]
            ),
            Extension(
                # os.path.join('python', 'modules', 'chain_bed_intersect'),
                'python.modules.chain_bed_intersect',
                [os.path.join('src', 'cython', 'chain_bed_intersect.pyx')]
            ),
            Extension(
                # os.path.join('python', 'modules', 'string_splitter'), 
                'python.modules.string_splitter',
                [os.path.join('src', 'cython', 'string_splitter.pyx')]
            )
        ],
        annotate=True,
        language_level=3
    )
    ],
    setup_requires=['cython']
)
