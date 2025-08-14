#!/usr/bin/env python3

"""
Checks for third-party dependencies necessary for the full TOGA2 experience
"""

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from importlib.util import find_spec
from shutil import which
from typing import (
    Any, Dict, Iterable, List, Optional, Tuple, Union
)

import click
import os
import subprocess
import sys

__author__ = 'Yury V. Malovichko'
__year__ = '2025'

CONTEXT_SETTINGS: Dict[str, Any] = {
    'help_option_names': [
        '-h', '-help', '--help'
    ],
    'ignore_unknown_options': True,
    'allow_extra_args': True,
    'max_content_width': 150
}

## minimal Python version accepted; so far the lowest tested with T2 was 3.9
MIN_MINOR_VERSION: int = 9

## requirement list file
REQUIREMENTS: str = 'requirements.txt'

## utilities essential for TOGA2 installation and minimal functionality
ESSENTIALS: Tuple[str] = (
    'awk', 'cargo', 'gcc'
)

NEXTFLOW: str = 'nextflow'
PARASOL: str = 'para'
NEXTFLOW_MANAGERS: Dict[str, Tuple[str]] = {
    'awsbatch' : ('aws',), 
    'azurebatch': ('az',), 
    'bridge': ('NONSENSICAL_BRIDGE_PLACEHOLDER',), ## TODO: Cannot find documentation on Bridge invocation from cc 
    'flux': ('flux',),
    'google-batch': ('gcloud',), 
    'condor': ('condor_submit',), 
    'hq': ('hq',), 
    'k8s': ('kubectl',),
    'lsf': ('bsub',), 
    'moab': ('msub',), 
    'nqsii': ('qsub',), 
    'oar': ('oarsub',),
    'pbs': ('qsub',), 
    'pbspro': ('qsub',), 
    'sge': ('qsub',), 
    'slurm': (
        'sacct', 'salloc', 'sattach', 
        'sbatch', 'sbcast', 'scancel', 
        'scontrol', 'sinfo', 'sprio',
        'squeue', 'srun', 'sshare',
        'sstat', 'strigger', 'sview'
    )
}

## Python packages to check; mirror the `requirements.txt` contents, with importing aliases considered
PACKAGES: Dict[str, Tuple[str, str]] = {
    'xgboost': ('xgboost', '3.0.0'), 
    'joblib': ('joblib', '1.4.2'), 
    'numpy': ('numpy', '2.2.5'), 
    'pandas': ('pandas', '2.2.3'),
    'numexpr': ('numexpr', '2.8.1'),
    'scikit-learn': ('sklearn', '1.7.1'), 
    'scipy': ('scipy', '1.15.2'),
    'Cython': ('cython', '3.0.12'), 
    'h5py': ('h5py', '3.12.1'),
    'click': ('click', '8.1.8'),
    'click-option-group': ('click_option_group', '0.5.6'),
    'contourpy': ('contourpy', '1.2.1'),
    'tensorflow': ('tensorflow', '2.8.0'),
}

class ThirdPartyChecker:

    __slots__ = ('name', 'min_version', 'path', 'cmd', 'version')

    def check(self) -> bool:
        if not self.check_presence():
            return False
        self.get_version()
        if not self.check_version():
            return False
        return True

    def check_presence(self) -> bool:
        src_in_path: Union[str, None] = which(self.name)
        if src_in_path is None:
            click.echo('WARNING: Program %s is missing from $PATH' % self.name)
        self.path = src_in_path
        return src_in_path is not None

    def get_version(self) -> str:
        pr: subprocess.Pipe = subprocess.Popen(
            self.cmd.format(self.path), 
            shell=True, 
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, _ = pr.communicate()
        rc: int = pr.returncode
        if rc != 0:
            return None
        stdout: str = stdout.decode('utf8')
        return stdout

    def check_version(self) -> bool:
        self.version: str = self.get_version()
        comps: Iterable[Tuple[str, str]] = zip(
            self.min_version.split('.'), self.version.split('.')
        )
        is_greater: bool = True
        for m, v in comps:
            m: int = int(m)
            v: int = int(v)
            if v > m:
                break
            if m > v:
                is_greater = False
        if not is_greater:
            click.echo(
                (
                    'WARNING: Version %s for program %s does not satisfy '
                    'the minimal version requirement %s' % (self.version, self.name, self.min_version)
                )
            )
        return is_greater


class PrankChecker(ThirdPartyChecker):

    def __init__(self) -> None:
        self.name: str = 'prank'
        self.min_version: str = '250331'
        self.path: Union[str, None] = None
        self.cmd: str = '{} -version'

    def get_version(self) -> str:
        version: Union[str, None] = None
        for line in super().get_version().split('\n'):
            version_point: int = line.find('PRANK')
            if version_point != -1:
                version =  line.rstrip()[version_point+8:].rstrip('.')
                return version
        raise Exception(
            'Unexpected "-version" formatting for PRANK instance at %s' % self.path
        )

class IqTreeChecker(ThirdPartyChecker):

    def __init__(self) -> None:
        self.name: str = 'iqtree2'
        self.min_version: str = '2.4.0'
        self.path: Union[str, None] = None
        self.cmd: str = '{} --version'

    def get_version(self) -> str:
        comps: List[str] = super().get_version().split(' ')
        if len(comps) < 4 or comps[2] != 'version':
            raise Exception(
                'Unexpected "--version" formatting for IqTree2 instance at %s' % self.path
            )
        return comps[3]


class IntronIcChecker(ThirdPartyChecker):

    def __init__(self) -> None:
        self.name: str = 'intronIC'
        self.min_version: str = '1.5.2'
        self.path: Union[str, None] = None
        self.cmd: str = '{} --version'

    def get_version(self) -> str:
        comps: List[str] = super().get_version().split(' ')
        if len(comps) < 2 or comps[0] != 'intronIC':
            raise Exception(
                'Unexpected "--version" formatting for intronIC instance at %s' % self.path
            )
        return comps[1].lstrip('v').rstrip()


class SpliceAiChecker(ThirdPartyChecker):

    def __init__(self) -> None:
        self.name: str = 'spliceai'
        self.min_version: str = '1.3.1'
        self.path: Union[str, None] = None
        self.cmd: str = '{} -h'

    def get_version(self) -> str:
        for line in super().get_version().split('\n'):
            if line.startswith('Version'):
                comps: List[str] = line.rstrip().split()
                if len(comps) < 2:
                    break
                version: str = comps[1]
                return version
        raise Exception(
            'Unexpected "--help" formatting for SpliceAI instance at %s' % self.path
        )


INSTALLABLES: List[ThirdPartyChecker] = [ 
    IntronIcChecker(), SpliceAiChecker(),
    PrankChecker(), IqTreeChecker(),
]

## missing dependency files
MISSING_PYTHON_PACKAGES: str = 'missing_packages.txt'
MISSING_THIRD_PARTY: str = 'missing_third_party.txt'

## installation commands for third-party software
# INTRONIC_INSTALL_CMD: str = 'git clone https://github.com/glarue/intronIC.git bin/intronIC'
INTRONIC_INSTALL_CMD: str = """
wget -P bin https://github.com/glarue/intronIC/archive/refs/tags/v1.5.2.tar.gz && \
tar -xzvf bin/v1.5.2.tar.gz -C bin/ && rm -rf bin/v1.5.2.tar.gz && mv bin/intronIC-1.5.2 bin/intronIC
"""
SPLICEAI_INSTALL_CMD: str = 'git clone https://github.com/Illumina/SpliceAI.git bin/SpliceAI && cd bin/SpliceAI && python setup.py install'
PRANK_INSTALL_CMD: str = """
git clone https://github.com/ariloytynoja/prank-msa.git bin/prank && cd bin/prank && make
"""
IQTREE2_INSTALL_CMD: str = """
wget -P bin/iqtree2 https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-intel.tar.gz && \
    tar -xzvf bin/iqtree2/iqtree-2.4.0-Linux-intel.tar.gz -C bin/iqtree2/
"""
INSTALL_CMDS: Dict[str, str] = {
    'intronIC': INTRONIC_INSTALL_CMD,
    'spliceai': SPLICEAI_INSTALL_CMD,
    'prank': PRANK_INSTALL_CMD,
    'iqtree2': IQTREE2_INSTALL_CMD
}

@click.group(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
def check_deps() -> None:
    """
    Non-Python dependencies check and installation module
    """
    pass

@check_deps.command(
    'essentials', 
    context_settings=CONTEXT_SETTINGS, 
    short_help='Checks for awk, gcc, and cargo'
)
def check_essentials() -> None:
    """
    Checks for essential dependencies necessary for TOGA2 installation (awk, gcc, and cargo).
    Absence of any of those makes TOGA2
    """
    for essential in ESSENTIALS:
        ess_in_path: Union[str, None] = which(essential)
        if ess_in_path is None:
            click.echo(
                'ERROR: Utility %s required for TOGA2 installation has not been found in $PATH' % essential
            )
            sys.exit(1)
        click.echo('Found %s at %s' % (essential, ess_in_path))
    click.echo('All essential utilities are present and available in $PATH')

@check_deps.command(
    'managers', 
    context_settings=CONTEXT_SETTINGS, 
    short_help='Checks for Nextflow, Parasol, and paralle process executors'
)
def check_managers() -> None:
    """
    Checks for parallel process managers. 
    By default, TOGA2 relies on Nextflow to invoke a third-party executor; the command will first check for 'nextflow' availability in $PATH
    and, once it is found, will further seek for any of those managers. If none of the Nextflow-supported managers are found, you can 
    still use Nextflow with a local executor.
    
    \b
    Alternatively, TOGA2 can use Parasol (https://genecats.gi.ucsc.edu/eng/parasol.html) as a standalone process manager. The command will seek for 'para'
    TOGA2 does NOT install parallel process manager on its own; if neither Nextflow nor Parasol are found, the command will throw an error and exit.

    \b
    NOTE: TOGA2 relies on default invoke command names (nextflow, sbatch, etc.) available in $PATH; currently it does not support alternative executor paths 
    and offers limited supported for custom parallel process management strategies. If you rely on either Nextflow- or Parasol-based strategy 
    and have the respective program already installed, please make sure it is available in $PATH under its default name.
    """
    nf_in_path: Union[str, None] = which(NEXTFLOW)
    nf_found: bool = nf_in_path is not None
    if nf_found:
        click.echo('%s found in path at %s' % (NEXTFLOW, nf_in_path))
        click.echo('Searching for available parallel process managers')
        managers_found: List[str] = []
        for manager, execs in NEXTFLOW_MANAGERS.items():
            click.echo('Checking for executor "%s" availability in $PATH' % manager)
            missing_execs: List[str] = []
            for ex in execs:
                ex_in_path: Union[str, None] = which(ex)
                if ex_in_path is None: 
                    missing_execs.append(ex)
            if missing_execs:
                click.echo(
                    'WARNING: The following executables for executor "%s" were not found: %s' % (
                        manager, ', '.join(missing_execs)
                    ) 
                )
                continue
            click.echo('Found all the necessary executables for executor "%s"' % manager)
            managers_found.append(manager)
        if not managers_found:
            click.echo(
                (
                    'WARNING: no parallel process managers compatible with Nextflow were found '
                    'in $PATH. If you intend to use TOGA2 with Nextflow, please install '
                    'an appropriate manager or default to local execution'
                )
            )
    para_in_path: Union[str, None] = which(PARASOL)
    para_found: bool = para_in_path is not None
    if para_found:
        click.echo('Found Parasol at %s' % para_in_path)
    if not nf_found and not para_found:
        click.echo(
            'ERROR: Neither Nextflow nor Parasol are available in $PATH. Please check '
            'if you have either manager installed and available in the environment intended '
            'to be used with TOGA2'
        )
        sys.exit(1)
    click.echo('=' * 30)
    click.echo('Parallel process manager check summary:')
    if nf_found:
        click.echo('* Nextflow is available in $PATH')
        click.echo('* The following executors are available for Nextflow: %s' % (
            'local, ' + ', '.join(managers_found))
        )
    else:
        click.echo(
            '* Nextflow was not found in $PATH. Use Parasol instead '
            'or install Nextflow on your machine'
        )
    if para_found:
        click.echo('* Para is available in $PATH')
    else:
        click.echo('* Parasol was not found in $PATH')
    

@check_deps.command(
    'python', 
    context_settings=CONTEXT_SETTINGS,
    short_help='Checks Python version and third-party package availability'
)
@click.option(
    '--installation_mode',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, suppresses missing/outdated package errors '
        'and lists them for further installation'
    )
)
def check_python(installation_mode: Optional[bool]) -> None:
    """
    Checks Python version and necessary packages availability
    """
    ## check the Python version
    minor_v: int = sys.version_info.minor
    if minor_v < MIN_MINOR_VERSION:
        click.echo(
            'ERROR: Used Python version is 3.%i'
        )

    packages_to_install: Dict[str, str] = {}
    for pack_name, (package, version_) in PACKAGES.items():
        if find_spec(package) is None:
            if not installation_mode:
                click.echo(
                    (
                        'ERROR: Module %s is not installed for the current Python version. '
                        'Please check your Python version, execute Makefile in "install" mode '
                        'or install the module manually'
                    ) % package
                )
                sys.exit(1)
            else:
                click.echo(
                    'WARNING: Module %s is not installed for the current Python version' % package
                )
                packages_to_install[pack_name] = version_
                continue
        try:
            ver: str = version(package)
        except PackageNotFoundError:
            pkg = import_module(package)
            ver: str = pkg.__version__
        click.echo('Found package %s, version %s' % (package, ver))
        for req, found in zip(version_.split('.'), ver.split('.')):
            req = int(req)
            found = int(found)
            if found > req:
                break
            if found < req:
                if installation_mode:
                    click.echo(
                        'WARNING: Package %s does not comply with the minimal required version %s' % (
                            package, version_
                        )
                    )
                    packages_to_install[pack_name] = version_
                    break
                else:
                    click.echo(
                        'ERROR: Found version %s for package %s with minimal required version %s' % (
                            ver, package, version_
                        )
                    )
                    sys.exit(1)
        if pack_name not in packages_to_install:
            click.echo(
                'Version %s for package %s satisfies the minimal requirement (%s)' % (
                    ver, package, version_
                )
            )
    if packages_to_install:
        click.echo(
            'Writing missing Python package names to %s ' % MISSING_PYTHON_PACKAGES
        )
        with open(MISSING_PYTHON_PACKAGES, 'w') as h:
            for package, version_ in packages_to_install.items():
                h.write(f'{package}>={version_}\n')
    
@check_deps.command(
    'third_party', 
    context_settings=CONTEXT_SETTINGS,
    short_help=(
        'Check for third-party tools required for auxiliary TOGA2 features '
        '(intronIC, SpliceAI, PRANK, and IqTree2)'
    )
)
def check_third_party() -> None:
    """
    Checks for the third-party tools used for TOGA2 input preparation and 
    annotation improvement. These include:\n
    \b
    \t* intronIC for reference intron classification;
    \t* SpliceAI for exon/intron boundary correction and query-specific intron annotation;
    \t* PRANK and IqTree2 for gene tree-based orthology resolution.


    The missing/outdated software will be saved to missing_third_party.txt file in the installation directory.
    The listed programs will be installed to bin/ once you run `make install` command.


    While having access to those programs is recommended for the full TOGA2 experience, 
    none of them are absolutely indispensable. If any of them are missing from your machine with no 
    chance to install them, consider the following workarounds:

    \t* Instead of installing intronIC and SpliceAI, you can use reference intron classification and query 
    SpliceAI annotations obtained elsewhere. For references/queries appearing in the TOGA2 companion dataset, 
    check https://hgdownload.soe.ucsc.edu/downloads.html ;

    \t* If you do not have PRANK>=v.250331 and/or IqTree2 and are not interested in gene tree-based 
    orthology refinement, disable this step by setting the `--skip_gene_trees` flag when runnning TOGA2 ;

    \t* Alternatively, you specify raxmlHPC-PTHREADS-AVX as phylogeny inference tool of choice when running TOGA2 
    by specifying `--use_raxml` flag. Note that currently raxmlHPC-PTHREADS-AVX is the only substitute for 
    IqTree2 supported by TOGA2, and its installation and configuration is relegated to the end user.
    """
    missing_soft: List[str] = []
    for checker in INSTALLABLES:
        found: bool = checker.check()
        if not found:
            missing_soft.append(checker.name)
            continue
        click.echo(
            'Found program "%s", version %s, in $PATH at %s' % (
                checker.name, checker.version, checker.path
            )
        )
    if missing_soft:
        click.echo('Writing names for missing programs to %s' % MISSING_THIRD_PARTY)
        with open(MISSING_THIRD_PARTY, 'w') as h:
            for missing in missing_soft:
                h.write(missing + '\n')


@check_deps.command(
    'install_third_party',
    context_settings=CONTEXT_SETTINGS,
    help='Installs missing/outdated third-party programs (intronIC, PRANK, IqTree2, SpliceAI)'
)
def install_third_party() -> None:
    """
    """
    if not os.path.exists(MISSING_THIRD_PARTY):
        click.echo(
            'All third-party programs have been successfully found in $PATH'
        )
        sys.exit(0)
    with open(MISSING_THIRD_PARTY, 'r') as h:
        for line in h:
            line = line.rstrip()
            if line not in INSTALL_CMDS:
                click.echo('ERROR: No installation recipe for third-party software %s' % line)
            cmd: str = INSTALL_CMDS[line]
            subprocess.run(cmd, shell=True)
            click.echo('Successfully installed %s' % line)

if __name__ == '__main__':
    check_deps()