#!/usr/bin/env python3

"""
Nucleotide/codon alignment module
"""

import os
import shutil
from collections import defaultdict
from shutil import which
from typing import Dict, List, Optional, TextIO, Tuple, Union

import click
from ete3 import Tree

from .cesar_wrapper_constants import CLASS_TO_NUM, LOSS_STATUSES
from .constants import Headers
from .shared import CommandLineManager, get_proj2trans, hex_code

__author__ = "Yury V. Malovichko"
__credits__ = "Bernhard Bein"
__year__ = "2025"

ALIGNERS_TO_USE: List[str] = ["macse2", "prank", "magus", "muscle"]
PRANK: str = "prank"
MACSE: str = "macse2"
MAGUS: str = "magus"
MUSCLE: str = "muscle"
TWOBIT2FA: str = "twoBitToFa"
ONE_TO_ONE: str = "one2one"
MANY_TO_ONE: str = "many2one"
SINGLE_COPY_CLASSES: Tuple[str] = (ONE_TO_ONE, MANY_TO_ONE)
ORTH_FILE: str = "orthology_classification.tsv"
LOSS_FILE: str = "loss_summary.tsv"
EXON_GZIP: str = "exon_aln.fa.gz"
EXON_GZIP_INDEX: str = ".exon_aln.fa.gz.ix"
TWOBIT_FILE: str = "exon_seqs.2bit"
EXTRACTION_ERROR: str = "Sequence extraction for {} from 2bit file {} failed"
EXON_SEP_DUMMY: str = "n"
PRANK_BEST_PRACTICE: str = "{} -d={} -F -DNA -o={} -seed={} "
PRANK_FOR_REGULAR_ALN: str = " -iterate=10"
PRANK_FOR_ANCESTRAL: str = " -once -showanc"
MACSE_BEST_PRACTICE: str = "{} -prog alignSequences -seq {} -out_NT {} -out_AA {}"
MAGUS_BEST_PRACTICE: str = "{} -i {} --datatype dna -o {}"
MUSCLE_ALIGN_CMD: str = "muscle -super5 {} -output {}"
MUSCLE_ENSEMBLE_CMD: str = "muscle -align {} -diversified -output {} -threads {}"
MUSCLE_MAXCC_CMD: str = "muscle -maxcc {} -output {}"
MUSCLE_LETTERCONF_CMD: str = "muscle -letterconf {} -ref {} -output {}"
# MUSCLE_ADDCONF_CMD: str = 'muscle -addconfseq {} -output {}'
DEV_NULL: str = "/dev/null"
ALN_ERROR: str = "Alignment with {} for transcript {}, exon {}, failed"
ONE2ZERO: str = "one2zero"
PROJECTION: str = "PROJECTION"


def filter_by_posterior(
    seqs: Dict[str, str], confidence: Dict[str, str], threshold: int
) -> str:
    """
    Given a MUSCLE alignment string, a letter confidence string,
    and a numeric threshold, return the alignment with positions having
    confidence lower than the threshold replaced with gap symbols
    """
    output: Dict[str, str] = {}
    for species, seq in seqs.items():
        scores: str = confidence.get(species, None)
        if scores is None:
            raise KeyError("Missing confidence scores for species %s" % species)
        if len(seq) != len(scores):
            raise ValueError(
                (
                    "Inconsistent score vector length for species %s: sequence consists of %i nucleotides, "
                    "score vector has %i positions"
                )
                % (species, len(seq), len(scores))
            )
        new_seq: str = ""
        for i, nuc in enumerate(seq):
            score: str = scores[i]
            if score.isdigit():
                val: int = int(score)
                if val < threshold:
                    nuc = "-"
            new_seq += nuc
        output[species] = new_seq
    return output


class CodonAligner(CommandLineManager):
    """ """

    __slots__: Tuple[str] = (
        "transcript",
        "exon_numbers",
        "aligner",
        "ref_exon_path",
        "ref_name",
        "loss_statuses",
        "use_one2zeros",
        "confidence_threshold",
        "muscle_threads",
        "tree",
        "aa_file",
        "seed",
        "show_ancestors",
        "ancestral_seq_dir",
        "aligner_exe",
        "twobit2fa",
        "tmp_dir",
        "keep_tmp",
        "output",
        "phase_split_codons",
        "add_projection_names",
        "confidence_score_file",
        "exon2phase",
        "exon2length",
        "exon_seqs",
        "concatenated_fasta",
        "confidence_scores",
        "exon2missing",
        "no_sequence_queries",
        "species2name",
        "clean_tmp",
    )

    @staticmethod
    def parse_exon_list(exon_numbers: str) -> str:
        """
        Splits comma-separated exon list
        """
        return [int(x) for x in exon_numbers.split(",") if x]

    @staticmethod
    def parse_loss_status(stat_list: str) -> Union[Tuple[str], None]:
        """Parses a comma-separated list of accepted loss statuses. Case-insensitive"""
        if stat_list is None:
            return None
        stats: List[str] = [x for x in stat_list.split(",") if x]
        unsupported_stats: List[str] = [x for x in stats if x not in LOSS_STATUSES]
        if unsupported_stats:
            raise ValueError(
                "The following loss statuses are not supported in TOGA2: %s"
                % ",".join(unsupported_stats)
            )
        return tuple(stats)

    def __init__(
        self,
        input_dirs: click.File,
        transcript_id: str,
        output: Optional[click.File],
        exon_numbers: Optional[Union[str, None]],
        reference_exons: Optional[Union[click.Path, None]],
        reference_name: Optional[Union[str, None]],
        accepted_loss_status: Optional[Union[str, None]],
        use_one2zeros: Optional[bool],
        confidence_threshold: Optional[int],
        aligner: Optional[str],
        aligner_exe: Optional[Union[click.Path, None]],
        phase_split_codons: Optional[bool],
        tree: Optional[Union[click.Path, None]],
        amino_acids_output: Optional[Union[click.Path, None]],
        show_ancestors: Optional[bool],
        path_to_ancestor_files: Optional[Union[click.Path, None]],
        seed: Optional[str],
        confidence_scores: Optional[Union[click.File, None]],
        muscle_threads: Optional[int],
        add_projection_names: Optional[bool],
        twobit2fa: Optional[Union[click.Path, None]],
        tmp_dir: Optional[Union[click.Path, None]],
        keep_tmp: Optional[bool],
        verbose: Optional[bool],
    ) -> None:
        self.v: bool = verbose
        self.set_logging()

        self.transcript: str = transcript_id
        self.exon_numbers: Union[List[str], None] = (
            None if exon_numbers is None else self.parse_exon_list(exon_numbers)
        )
        self.loss_statuses: Union[Tuple[str], None] = self.parse_loss_status(
            accepted_loss_status
        )
        self.use_one2zeros: bool = use_one2zeros
        self.confidence_threshold: int = confidence_threshold
        self.muscle_threads: int = muscle_threads
        self.aligner: str = aligner
        self.set_aligner_exe(aligner, aligner_exe)
        self.set_twobit_path(twobit2fa)
        self.tree: Union[str, None] = tree
        self.aa_file: Union[str, None] = amino_acids_output
        if amino_acids_output is not None and aligner != MACSE:
            self._to_log(
                "Amino acid output is not supported for aligners other than MACSE2",
                "warning",
            )
        self.show_ancestors: bool = show_ancestors
        self.ancestral_seq_dir: Union[click.Path, None] = path_to_ancestor_files
        self.seed: str = seed
        self.add_projection_names: bool = add_projection_names
        self.ref_exon_path: Union[str, None] = reference_exons
        self.ref_name: Union[str, None] = reference_name
        if self.ref_exon_path is not None and self.ref_name is None:
            self.ref_name = "ref"
        self.exon2phase: Tuple[int] = tuple()
        self.exon2length: Dict[int, int] = {}
        self.exon_seqs: Dict[int, List[str]] = defaultdict(list)
        self.concatenated_fasta: Dict[str, str] = defaultdict(str)
        self.confidence_scores: Dict[str, str] = defaultdict(str)
        self.exon2missing: Dict[int, List[str]] = defaultdict(list)
        self.no_sequence_queries: List[str] = []
        self.species2name: Dict[str, str] = {}
        self.tmp_dir: str = tmp_dir
        if not os.path.exists(tmp_dir):
            self._mkdir(self.tmp_dir)
            self.clean_tmp: bool = True
        else:
            self.clean_tmp: bool = False
        self.keep_tmp: bool = keep_tmp
        self.output: TextIO = output
        self.phase_split_codons: bool = phase_split_codons or aligner == MACSE
        self.confidence_score_file: TextIO = confidence_scores
        self.run(input_dirs)

    def run(self, input_dirs: TextIO) -> None:
        """Main method"""
        self._to_log("Input transcript is: %s" % self.transcript)
        if self.ancestral_seq_dir is not None:
            self._to_log(
                "Creating a storage directory for ancestral sequence reconstructions"
            )
            self._mkdir(self.ancestral_seq_dir)
        if self.ref_exon_path:
            self._to_log("Extracting reference sequences from %s" % self.ref_exon_path)
            ref_sequences: List[str] = self.extract_sequences(
                self.ref_exon_path, self.transcript, is_ref=True, infer_phases=True
            )
            for i, exon in enumerate(ref_sequences):
                if self.exon_numbers is not None:
                    exon_num: int = self.exon_numbers[i]
                else:
                    exon_num: int = i + 1
                header: str = f">{self.ref_name}_|_{self.transcript}_|_{exon_num}"
                entry: str = f"{header}\n{exon}"
                self.exon_seqs[exon_num].append(entry)
            self._to_log(
                "Successfully extracted %i exons for reference %s"
                % (len(ref_sequences), self.ref_name)
            )
        non_missing_queries: int = 0
        for i, line in enumerate(input_dirs):
            skip_loss_check: bool = False
            path: str = line.rstrip()
            if not path:
                continue
            species: str = (
                path.rstrip(os.sep).split(os.sep)[-1].replace("vs_", "")
            )  # s.path.basename(path).replace('vs_', '')
            projection_names, is_one2zero = self.get_orthologs(path)
            num_proj_found: int = len(projection_names)
            if not num_proj_found and self.use_one2zeros and is_one2zero:
                projection_names = self.best_candidate_in_loss_file(path)
                num_proj_found: int = len(projection_names)
                if num_proj_found > 1:
                    self._to_log(
                        (
                            "Multiple one2zero projections for query %s; "
                            "excluding it from the analysis"
                        ) % species,
                        "warning"
                    )
                    num_proj_found = 0
            if not num_proj_found:
                self._to_log("No projections found for species %s" % species, "warning")
                for exon in self.exon2length:
                    self.exon2missing[exon].append(species)
                self.no_sequence_queries.append(species)
                continue
            self._to_log(
                "Found %i projection names for species %s" % (num_proj_found, species)
            )
            non_missing_queries += 1
            if num_proj_found > 1:
                self._to_log(
                    "Duplicated entry for projection %s in species %s"
                    % (self.transcript, species),
                    "warning",
                )
                # self.no_sequence_queries.append(species)
                # continue
            ref_found: bool = (
                species == self.ref_name or i == 0
            ) and self.ref_exon_path is None
            if ref_found:
                self.ref_name = species
            rejected: bool = False
            if self.loss_statuses and not skip_loss_check:
                # accepted_loss_status: bool = self.check_loss_status(path, proj)
                # if not accepted_loss_status:
                #     self._to_log(
                #         (
                #             "Projection %s does not comply with the accepted lost statuses "
                #             "provided in species %s; accepted statuses are: %s"
                #         )
                #         % (proj, species, ",".join(self.loss_statuses)),
                #         "warning",
                #     )
                #     self.no_sequence_queries.append(species)
                #     rejected = True
                status2proj: Dict[str, str] = self.check_loss_status(path, projection_names)
                highest_status: str = max(status2proj.keys(), key=lambda x: CLASS_TO_NUM[x])
                if highest_status not in self.loss_statuses:
                    self._to_log(
                        (
                            "Neither of the found projections comply with the accepted lost statuses"
                            "in species %s; accepted statuses are: %s"
                        )
                        % (species, ",".join(self.loss_statuses)),
                        "warning",
                    )
                    self.no_sequence_queries.append(species)
                    rejected = True
                projection_names: List[str] = status2proj[highest_status]
            proj: str = min(projection_names, key=lambda x: int(x.split("#")[-1].split(",")[0]))
            self.species2name[species] = proj
            exons: Dict[str, int] = self.extract_sequences(
                path, proj, infer_phases=ref_found
            )
            for i, exon in enumerate(exons):
                if self.exon_numbers is not None:
                    exon_num: int = self.exon_numbers[i]
                else:
                    exon_num: int = i + 1
                if not exon or rejected:
                    self._to_log(
                        "Exon %i for species %s is fully deleted" % (exon_num, species)
                    )
                    self.exon2missing[exon_num].append(species)
                    continue
                header: str = f">{species}_|_{self.transcript}_|_{exon_num}"
                entry: str = f"{header}\n{exon}"
                self.exon_seqs[exon_num].append(entry)
            self._to_log(
                "Extracted %i exon sequences for species %s" % (len(exons), species)
            )
        if non_missing_queries == 0:
            self._to_log(
                "No 1:1 orthology found in any query for transcript %s" % self.transcript,
                "warning"
            )
            self._exit()
        self._to_log("Proceeding to alignment step")
        self.run_alignment()
        self._to_log("Writing output fasta to %s" % self.output.name)
        self.write_output()

    def set_aligner_exe(self, aligner: str, exe: Union[str, None]) -> None:
        """
        Sets aligner executable; if not provided, will search for executable in PATH
        """
        if exe is not None:
            self.aligner_exe: str = exe
            return
        exe_in_path: str = which(aligner)
        if exe_in_path is not None:
            self._to_log("Found aligner %s in PATH at %s" % (aligner, exe_in_path))
            self.aligner_exe: str = exe_in_path
            return
        self._die(
            "Aligner %s was not found in PATH, with no default path provided" % aligner
        )

    def set_twobit_path(self, twobit_exe: str) -> None:
        """
        Sets UCSC twoBitToFa executable; if not provided, will search for executable in PATH
        """
        if twobit_exe is not None:
            self.twobit2fa: str = twobit_exe
            return
        exe_in_path: str = which(TWOBIT2FA)
        if exe_in_path is not None:
            self._to_log("Found twoBitToFa in PATH at %s" % exe_in_path)
            self.twobit2fa: str = exe_in_path
            return
        self._die(
            "twoBitToFa executable was not found in PATH, with no default path provided"
        )

    def get_orthologs(self, path: str) -> Tuple[List[str], bool]:
        """
        Gets the one:one and many:one projection names for the focal transcript and a given query
        """
        table_file: str = os.path.join(path, ORTH_FILE)
        if not os.path.exists(table_file):
            self._die(
                "Orthology classification file does  not exist for input directory %s"
                % path
            )
        out_projs: List[str] = []
        is_one2zero: bool = False
        with open(table_file, "r") as h:
            for i, line in enumerate(h, start=1):
                if line == Headers.ORTHOLOGY_TABLE_HEADER:
                    continue
                data: List[str] = line.rstrip().split("\t")
                if len(data) != 5:
                    self._die(
                        "Numbers of fields at line %i in file %s differs from the expected (5)"
                        % (i, table_file)
                    )
                ref_tr: str = data[1]
                if ref_tr != self.transcript:
                    continue
                orth_status: str = data[4]
                if orth_status not in SINGLE_COPY_CLASSES:
                    if orth_status == ONE2ZERO:
                        is_one2zero = True
                    continue
                proj_name: str = data[3]
                out_projs.append(proj_name)
        return out_projs, is_one2zero

    def best_candidate_in_loss_file(self, path: str) -> List[str]:
        """
        Picks the best non-functional isoform from the loss file.

        Args:
            path: a path to the input directory

        Returns:
            A list containing a single projection name, or an empty list if either no 
            or multiple projections were found
        """
        loss_file: str = os.path.join(path, LOSS_FILE)
        if not os.path.exists(loss_file):
            self._die(
                "File %s is missing from the input directory %s" % (LOSS_FILE, path)
            )
        out_list: List[str] = []
        curr_status: str = ""
        with open(loss_file, "r") as h:
            for i, line in enumerate(h, start=1):
                line = line.rstrip()
                if line == Headers.LOSS_FILE_HEADER:
                    continue
                data: List[str] = line.split("\t")
                if not data or not data[0]:
                    continue
                if len(data) < 3:
                    self._die(
                        "Improperly formatting found in file %s at line %i"
                        % (loss_file, i)
                    )
                if data[0] != PROJECTION:
                    continue
                tr: str = get_proj2trans(data[1])[0]
                if tr != self.transcript:
                    continue
                if "#paralog" in data[1] or "#retro" in data[1]:
                    continue
                if not out_list:
                    out_list.append(data[1])
                    curr_status = data[2]
                    continue
                status: int = CLASS_TO_NUM[data[2]]
                prev_status: int = CLASS_TO_NUM[curr_status]
                if status == prev_status:
                    out_list.append(data[1])
                if status > prev_status:
                    out_list = [data[1]]
                    curr_status = data[2]
        return out_list


    # def check_loss_status(self, path: str, proj: str) -> bool:
    def check_loss_status(self, path: str, projections: List[str]) -> Dict[str, str]:
        """
        Checks whether an orthologous sequence has an accepted loss status
        defined by the user
        """
        loss_file: str = os.path.join(path, LOSS_FILE)
        if not os.path.exists(loss_file):
            self._die(
                "File %s is missing from the input directory %s" % (LOSS_FILE, path)
            )
        status2proj: Dict[str, List[str]] = defaultdict(list)
        with open(loss_file, "r") as h:
            for i, line in enumerate(h, start=1):
                line = line.rstrip()
                if line == Headers.LOSS_FILE_HEADER:
                    continue
                data: List[str] = line.split("\t")
                if not data or not data[0]:
                    continue
                if len(data) < 3:
                    self._die(
                        "Improperly formatting found in file %s at line %i"
                        % (loss_file, i)
                    )
                # if data[1] != proj:
                #     continue
                # return data[2] in self.loss_statuses
                if data[1] not in projections:
                    continue
                status2proj[data[2]].append(data[1])
        return status2proj

    def extract_sequences(
        self,
        path: str,
        projection: str,
        is_ref: bool = False,
        infer_phases: bool = False,
    ) -> List[str]:
        """
        Extracts exon sequences from the TwoBit sequence storage
        """
        if not is_ref:
            twobit_path: str = os.path.join(path, TWOBIT_FILE)
            if not os.path.exists(twobit_path):
                self._die(
                    "TwoBit file %s is missing for input directory %s"
                    % (TWOBIT_FILE, path)
                )
        else:
            twobit_path: str = path
        projection = projection.replace(",", ".")
        cmd: str = f"{self.twobit2fa} -seq={projection} {twobit_path} stdout"
        err: str = EXTRACTION_ERROR.format(projection, twobit_path)
        entry: str = self._exec(cmd, err)
        out_seqs: List[str] = []
        curr_seq: str = ""
        ex_num: int = 1
        if infer_phases:
            phases: List[int] = []
            prev_phase: int = 0
        for line in entry.split("\n"):
            if not line:
                continue
            if line[0] == ">":
                continue
            ex_seqs: List[str] = line.split(EXON_SEP_DUMMY)
            for i, ex_seq in enumerate(ex_seqs):
                ## if the nex exon's sequence is encountered,
                ## add the previous exon to the list unless it was excluded by the user
                if i:
                    if self.exon_numbers is None or ex_num in self.exon_numbers:
                        out_seqs.append(curr_seq)
                        if infer_phases:
                            self.exon2length[ex_num] = len(curr_seq)
                            phase: int = (len(curr_seq) - prev_phase) % 3
                            phases.append(phase)
                            prev_phase = phase
                    curr_seq = ""
                    ex_num += 1
                ## if the current exon was excluded by the user, proceed further
                if self.exon_numbers is not None and ex_num not in self.exon_numbers:
                    continue
                ## otherwise, add the subsequence to the current exon's sequence
                curr_seq += ex_seq
        if self.exon_numbers is None or ex_num in self.exon_numbers:
            out_seqs.append(curr_seq)
            if infer_phases:
                self.exon2length[ex_num] = len(curr_seq)
        if infer_phases:
            self.exon2phase = tuple(phases)
        out_seqs = self.phase_codons(out_seqs)
        return out_seqs

    def phase_codons(self, exons: List[str]) -> List[str]:
        """
        Given a list of exon sequences, restores split codons to zero split phase
        """
        if not self.exon2phase:
            return exons
        if not self.exon2length:
            return exons
        if len(exons) == 1:
            return exons
        ## a pair of adjacent exons is a minimal unit of codon restoration,
        ## therefore iteration starts from the second exon
        for e in range(1, len(exons)):
            phase: int = self.exon2phase[e - 1]
            ## intron separating the two exons is in zero phase;
            ## nothing to restore, thus proceed further
            if not phase:
                continue
            prev_exon: str = exons[e - 1]
            next_exon: str = exons[e]
            ## extract the reference exon lengths; remember that the numeration is one-based
            prev_len: int = self.exon2length[e]
            next_len: int = self.exon2length[e + 1]
            phase_remainder: int = 3 - phase
            ## if previous exon is deleted, crop the split portion from the next one
            if not len(prev_exon):
                next_exon = next_exon[phase_remainder:]
                exons[e] = next_exon
                continue
            ## if next exon is deleted, crop the split portion from the previous one
            if not len(next_exon):
                prev_exon = prev_exon[:-phase]
                exons[e - 1] = prev_exon
                continue
            ## otherwise, resolve the split codons in favour of the shorter of two exons
            if self.phase_split_codons:
                if prev_len < next_len:
                    if phase_remainder >= len(next_exon):
                        continue
                    prev_exon = prev_exon + next_exon[:phase_remainder]
                    next_exon = next_exon[phase_remainder:]
                    exons[e - 1] = prev_exon
                    exons[e] = next_exon
                else:
                    if phase >= len(next_exon):
                        continue
                    # _prev_exon = prev_exon[:-phase]
                    # next_exon = prev_exon[-phase:] + next_exon
                    exons[e - 1] = prev_exon[:-phase]  ##_prev_exon
                    exons[e] = prev_exon[-phase:] + next_exon  ##next_exon
        return exons

    def run_alignment(self) -> None:
        """
        Runs alignment command for each exon consecutively
        """
        for exon, seqs in self.exon_seqs.items():
            tmp_fasta_in_name: str = (
                f"{self.transcript}_{exon}_input_{hex_code()}.fasta"
            )
            tmp_fasta_in_path: str = os.path.join(self.tmp_dir, tmp_fasta_in_name)
            with open(tmp_fasta_in_path, "w") as h:
                for seq in seqs:
                    h.write(seq + "\n")
            tmp_fasta_out_name: str = (
                f"{self.transcript}_{exon}_output_{hex_code()}.fasta"
            )
            tmp_fasta_out_path: str = os.path.join(self.tmp_dir, tmp_fasta_out_name)
            anc_fa: str = ""
            anc_dnd: str = ""
            aln_format: Tuple[str, str, str] = (
                self.aligner_exe,
                tmp_fasta_in_path,
                tmp_fasta_out_path,
            )
            tmp_tree_path: str = ""
            if self.tree is not None:
                self._to_log("Pruning the guiding tree for alignment")
                tmp_tree_file: str = f"tree_{self.transcript}_{exon}_{hex_code()}"
                tmp_tree_path: str = os.path.join(self.tmp_dir, tmp_tree_file)
                specs_to_remove: List[str] = self.exon2missing.get(exon, [])
                postfix: str = f"_|_{self.transcript}_|_{exon}"
                self.prune_tree(specs_to_remove, postfix, tmp_tree_path)
            aa_file: str = ""
            confidence_file: str = ""
            if self.aligner == MUSCLE:
                ## MUSCLE was introduced primarily for the sake of column confidence scores report,
                ## therefore its use requires special handling
                self._to_log("Running the MUSCLE routine for exon %s" % exon)
                self.muscle_alignment(tmp_fasta_in_path, tmp_fasta_out_path)
            else:
                ## MACSE and PRANK behave similarly
                if self.aligner == MACSE:
                    if self.aa_file is None:
                        aa_file: str = os.path.join(
                            self.tmp_dir, "macse_aa_aln" + hex_code()
                        )
                    else:
                        aa_file: str = (
                            self.aa_file
                        )  ## TODO: Must be created for each exon
                    cmd: str = MACSE_BEST_PRACTICE.format(*aln_format, aa_file)
                elif self.aligner == PRANK:
                    cmd: str = PRANK_BEST_PRACTICE.format(*aln_format, self.seed)
                    if self.tree is not None:
                        cmd += f" -t={tmp_tree_path} -prunetree"
                    if self.show_ancestors:
                        cmd += PRANK_FOR_ANCESTRAL
                        anc_fa = tmp_fasta_out_path + ".best.anc.fas"
                        anc_dnd = tmp_fasta_out_path + ".best.anc.dnd"
                    else:
                        cmd += PRANK_FOR_REGULAR_ALN
                    tmp_fasta_out_path += ".best.fas"
                elif self.aligner == MAGUS:
                    cmd: str = MAGUS_BEST_PRACTICE.format(*aln_format)
                self._to_log("Running alignment for exon %s" % exon)
                self._echo(f"Alignment command: {cmd}")
                self.regular_alignment(
                    cmd, 
                    exon, 
                    tmp_fasta_in_path,
                    tmp_fasta_out_path,
                )
                self._exec(cmd, ALN_ERROR.format(self.aligner, self.transcript, exon))
                if self.aligner == MACSE and self.aa_file is None:
                    self._rm(aa_file)
            max_len: int = 0

            ## parse the output alignment file, record the aligned lines per species
            with open(tmp_fasta_out_path, "r") as h:
                species: str = ""
                seq: str = ""
                for line in h:
                    line = line.rstrip()
                    if not line:
                        continue
                    if line[0] == ">":
                        ## header found
                        if species:
                            self.concatenated_fasta[species] += seq
                            max_len: int = max(max_len, len(seq))
                            species = ""
                            seq = ""
                        species: str = line[1:].split("_|_")[0]
                        continue
                    seq += line
                if seq:
                    self.concatenated_fasta[species] += seq
                    max_len: int = max(max_len, len(seq))
            ## if MUSCLE was used for sequence alignment, record the column confidence scores in the same fashion
            if self.aligner == MUSCLE:
                confidence_file: str = tmp_fasta_in_path + ".letterconf.afa"
                with open(confidence_file, "r") as h:
                    species: str = ""
                    confidence: str = ""
                    for line in h:
                        line = line.rstrip()
                        if not line:
                            continue
                        if line[0] == ">":
                            ## header found
                            if species:
                                self.confidence_scores[species] += confidence
                                species = ""
                                confidence = ""
                            species: str = line[1:].split("_|_")[0]
                            continue
                        confidence += line
                    if confidence:
                        self.confidence_scores[species] += confidence
            if exon in self.exon2missing:
                for species in self.exon2missing[exon]:
                    if species in self.no_sequence_queries:
                        continue
                    self.concatenated_fasta[species] += "-" * max_len
                    if self.aligner == MUSCLE:
                        self.confidence_scores[species] += "0" * max_len
            if self.aligner == MUSCLE:
                self.concatenated_fasta = filter_by_posterior(
                    self.concatenated_fasta,
                    self.confidence_scores,
                    self.confidence_threshold,
                )  ## TODO: Rewrite as
            if self.show_ancestors:
                if self.ancestral_seq_dir is not None:
                    anc_fa_out: str = os.path.join(
                        self.ancestral_seq_dir,
                        f"{self.transcript}_exon{exon}.ancestral.fasta",
                    )
                    anc_dnd_out: str = os.path.join(
                        self.ancestral_seq_dir,
                        f"{self.transcript}_exon{exon}.ancestral.dnd",
                    )
                    shutil.move(anc_fa, anc_fa_out)
                    shutil.move(anc_dnd, anc_dnd_out)
            if not self.keep_tmp:
                self._rm(tmp_fasta_in_path)
                self._rm(tmp_fasta_out_path)
                if self.tree is not None:
                    self._rm(tmp_tree_path)

    def regular_alignment(
        self, 
        cmd: str, 
        exon: int,
        in_fasta: str,
        out_fasta: str
    ) -> None:
        """
        Staple sequence alignment method for PRANK/MACSE2/MAGUS aligners.
        If only a single sequence is found in the input file, the whole directory

        Args:
            cmd: command to execute
            exon: exon number to align; used in error message

        Returns:
            None

        Raises:
            Dies with self._die() if alignment command dies
        """
        present_seq_num: int = len(self.exon_seqs[exon])
        if present_seq_num < 2:
            self._to_log(
                "One or less input sequence for exon %s; reporting the alignment as"
            )
            copy_cmd: str = f"cp {in_fasta} {out_fasta}"
            _ = self._exec(copy_cmd, "Copying %s to %s failed" % (in_fasta, out_fasta))
            return
        self._exec(cmd, ALN_ERROR.format(self.aligner, self.transcript, exon))

    def muscle_alignment(self, seq_input: str, output: str) -> None:
        """
        MUSCLE supports per-column alignment confidence score extraction,
        which requires running a series of MUSCLE modes
        """
        ## step 1: acquire the diversified alignment ensemble (i.e., align the sequences)
        ensemble_file: str = seq_input + ".efa"
        aln_cmd: str = MUSCLE_ENSEMBLE_CMD.format(
            seq_input, ensemble_file, self.muscle_threads
        )
        self._to_log("Running the alignment procedure with MUSCLE")
        _ = self._exec(aln_cmd, "MUSCLE alignment for %s failed" % seq_input)

        ## step 2: extract the maximum confidence alignment (a reference)
        maxcc_cmd: str = MUSCLE_MAXCC_CMD.format(ensemble_file, output)
        self._to_log(
            "Extracting the most probable alignment from the MUSCLE alignment ensemble"
        )
        _ = self._exec(
            maxcc_cmd,
            "Failed to extract the most probable alignment from %s" % ensemble_file,
        )

        ## step 3: calculate column confidence scores
        confidence_file: str = seq_input + ".letterconf.afa"
        letterconf_cmd: str = MUSCLE_LETTERCONF_CMD.format(
            ensemble_file, output, confidence_file
        )
        _ = self._exec(
            letterconf_cmd,
            "Extracting column confidence intervals from %s failed" % ensemble_file,
        )
        self._rm(ensemble_file)

    def prune_tree(self, species: List[str], postfix: str, dest: str) -> None:
        """Removes species listed in the `species` argument from the tree"""
        in_tree = Tree(self.tree, format=1)
        leaves_to_leave: List[str] = [
            x.name for x in in_tree.iter_leaves() if x.name not in species
        ]
        in_tree.prune(leaves_to_leave, preserve_branch_length=True)
        for leaf in in_tree.iter_leaves():
            leaf.name = f"{leaf.name}{postfix}"
        in_tree.write(format=1, outfile=dest)

    def write_output(self) -> None:
        """Writes concatenated exonwise alignments into an output Fasta file"""
        for species, seq in self.concatenated_fasta.items():
            if species == self.ref_name:
                species = "REFERENCE"
            header: str = f">{species}"
            if self.add_projection_names:
                proj_name: Union[str] = self.species2name.get(species)
                if proj_name:
                    header = f"{header}.{proj_name}"
            self.output.write(header + "\n" + seq + "\n")
        if self.aligner == MUSCLE:
            for species, score in self.confidence_scores.items():
                if species == self.ref_name:
                    species = "REFERENCE"
                header: str = f">{species}"
                self.confidence_score_file.write(header + "\n" + score + "\n")

    def set_logging(self) -> None:
        """Overriding logging setup module"""
        super().set_logging()


if __name__ == "__main__":
    CodonAligner()
