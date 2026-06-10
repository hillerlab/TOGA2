#!/usr/bin/env python3

"""
SpliceAI prediction wrapper
"""

import os
from collections import defaultdict
from heapq import heappop, heappush

# from parallel_jobs_manager import (
#     CustomStrategy, NextflowStrategy,
#     ParaStrategy, ParallelJobsManager
# )
from typing import Dict, List, Optional, Tuple, Union

import click

from .constants import Constants
from .parallel_jobs_manager import (
    CustomStrategy,
    NextflowStrategy,
    ParallelJobsManager,
    ParaStrategy,
)
from .shared import CommandLineManager, dir_name_by_date, get_upper_dir


## GLOBAL TODO: Move all the constants to a separate storage class in constants.py
PYTHON_DIR: str = get_upper_dir(__file__, 2)
EXEC_SCRIPT: str = os.path.join(PYTHON_DIR, "predict_with_spliceai.py")
CONTIG_SIZE_SCRIPT: str = os.path.join(PYTHON_DIR, "get_contig_sizes.py")

DEFAULT_CHUNK_SIZE: int = 6_000_000
DEFAULT_FLANK_SIZE: int = 50_000
DEFAULT_MIN_CONTIG_SIZE: int = 500
DEFAULT_MEMORY_LIMIT: int = 5

STRANDS: Tuple[str, str] = ("+", "-")

FILE_NAME_TEMPLATES: Tuple[str] = (
    "AcceptorPlus",
    "AcceptorMinus",
    "DonorPlus",
    "DonorMinus",
)

DONOR_PLUS: str = "spliceAiDonorPlus.bw"
DONOR_MINUS: str = "spliceAiDonorMinus"
ACC_PLUS: str = "spliceAiAcceptorPlus"
ACC_MINUS: str = "spliceAiAcceptorMinus"

## TODO: Add --resume_from and --halt_at functionality here and to toga2.py
PIPELINE_STEPS: Tuple[str] = ("all", "prepare", "schedule", "run", "aggregate")
RESUME_ORDER: Dict[str, str] = {x: i for i, x in enumerate(PIPELINE_STEPS)}

def batch_num(file: str, template: str) -> int:
    """Returns the batch number from the batch file name

    Args:
        file: a name of the batch output file
        template: SpliceAi file template to remove

    Returns:
        Batch number as integer
    """
    return file.split("@")[-1].replace(template, "").replace("batch", "").replace(".wig", "")


class SpliceAiManager(CommandLineManager):
    """ """

    __slots__ = (
        "output",
        "tmp_dir",
        "nextflow_dir",
        "twobit",
        "chunk_size",
        "flank_size",
        "min_contig_size",
        "round_to",
        "min_prob",
        "project_name",
        "job_num",
        "parallel_strategy",
        "nextflow_exec_script",
        "max_number_of_retries",
        "nextflow_config_file",
        "max_parallel_time",
        "cluster_queue_name",
        "memory_limit",
        "bed_dir",
        "job_file",
        "job_list",
        "chunk_num",
        "tmp_fa",
        "unmasked_twobit",
        "chrom_sizes",
        "resume_from",
        "halt_at",
        "twobittofa_binary",
        "fatotwobit_binary",
        "wigtobigwig_binary",
        "v",
        "log_file",
        "keep_tmp",
    )

    def __init__(
        self,
        query_2bit: click.Path,
        output: Optional[click.Path],
        chunk_size: Optional[int],
        flank_size: Optional[int],
        min_contig_size: Optional[int],
        round_to: Optional[int],
        min_prob: Optional[float],
        job_number: Optional[int],
        parallel_strategy: Optional[str],
        nextflow_exec_script: Optional[click.Path],
        max_number_of_retries: Optional[int],
        nextflow_config_file: Optional[click.Path],
        max_parallel_time: Optional[int],
        cluster_queue_name: Optional[str],
        memory_limit: Optional[int],
        resume_from: Optional[str],
        halt_at: Optional[str],
        twobittofa_binary: Optional[click.Path],
        fatotwobit_binary: Optional[click.Path],
        wigtobigwig_binary: Optional[click.Path],
        project_name: Optional[Union[str, None]],
        keep_temporary_files: Optional[bool],
        verbose: Optional[bool],
    ) -> None:
        self.v: bool = verbose
        self.project_name: str = (
            project_name if project_name is not None else dir_name_by_date("spliceai")
        )
        self.output: str = (
            self._abspath(output)
            if output is not None
            else self._abspath(self.project_name)
        )
        self._mkdir(self.output)
        # self.tmp_dir: str = os.path.join(self.output, f"tmp_{self.project_name}")
        self.tmp_dir: str = os.path.join(self.output, "tmp")
        self.nextflow_dir: str = os.path.join(self.tmp_dir, "nextflow")
        self.log_file: str = os.path.join(self.output, f"{self.project_name}.txt")
        self.set_logging()
        self._to_log("Initializing SpliceAI annotation module")
        self.logger.propagate = False

        self.resume_from: str = resume_from
        self.halt_at: str = halt_at

        self.twobit: click.Path = query_2bit

        self.chunk_size: int = chunk_size
        self.flank_size: int = flank_size
        if self.flank_size > self.chunk_size:
            self._to_log(
                (
                    "Chunk flank cannot be greater than chunk size; "
                    "setting chunk flank to %i"
                )
                % self.chunk_size,
                "warning",
            )
            self.flank_size = self.chunk_size
        self.min_contig_size: int = min_contig_size

        self.round_to: int = round_to
        self.min_prob: float = min_prob

        self.job_num: int = job_number
        self.parallel_strategy: str = parallel_strategy
        self.nextflow_exec_script: str = nextflow_exec_script
        self.max_number_of_retries: int = max_number_of_retries
        self.nextflow_config_file: str = nextflow_config_file
        self.max_parallel_time: int = max_parallel_time
        self.cluster_queue_name: str = cluster_queue_name
        self.memory_limit: int = memory_limit
        self.chunk_num: int = 0

        self.bed_dir: str = os.path.join(self.tmp_dir, "bed_input")
        self.job_list: str = os.path.join(self.tmp_dir, "joblist")
        self.tmp_fa: str = os.path.join(self.tmp_dir, "unmasked.fa")
        self.unmasked_twobit: str = os.path.join(self.tmp_dir, "unmasked.2bit")
        self.chrom_sizes: str = os.path.join(self.tmp_dir, "chrom_sizes.txt")

        self.twobittofa_binary: Union[str, None] = twobittofa_binary
        self.fatotwobit_binary: Union[str, None] = fatotwobit_binary
        self.wigtobigwig_binary: Union[str, None] = wigtobigwig_binary
        self.keep_tmp: bool = keep_temporary_files

        self._mkdir(self.output)
        self._mkdir(self.tmp_dir)
        self._mkdir(self.bed_dir)
        self._mkdir(self.nextflow_dir)

        self.run()

        if not keep_temporary_files:
            self._rmdir(self.tmp_dir)

    def run(self) -> None:
        """Entry point"""
        if self._execute_step("prepare"):
            self._to_log("Preparing input data for SpliceAI")
            self.prepare_ref_genome()
        else:
            if self.halt_at == "prepare":
                self._to_log(
                    "Finishing SpliceAI annotation before the data preparation step as suggested"
                )
                self._exit()
            self._to_log("Skipping the data preparation step as suggested")
        if self._execute_step("schedule"):
            self._to_log("Scheduling parallel jobs")
            self.schedule_jobs()
        else:
            if self.halt_at == "schedule":
                self._to_log(
                    "Finishing SpliceAI annotation before the job scheduling step as suggested"
                )
                self._exit()
            self._to_log("Skipping the job scheduling step as suggested")
        if self._execute_step("run"):
            self._to_log("Preparing to run the parallel SpliceAI jobs")
            self.run_jobs()
        else:
            if self.halt_at == "run":
                self._to_log(
                    "Finishing SpliceAI annotation before the execution step as suggested"
                )
                self._exit()
            self._to_log("Skipping the execution step as suggested")
        if self._execute_step("aggregate"):
            self._to_log("Aggregating the results")
            self.aggregate_jobs()
        else:
            if self.halt_at == "aggregate":
                self._to_log(
                    "Finishing SpliceAI annotation before the results aggregation step as suggested"
                )
                self._exit()
            self._to_log("Skipping the results aggregation step as suggested")

    def prepare_ref_genome(self) -> None:
        """
        Replaces assembly gaps with polyA in the input genome
        since SpliceAi cannot process symbols other than A/C/G/T

        As in-place 2bit file modification is impossible, the original file
        is first decompressed, with non-standard nucleotide symbols
        in the sequence part then being removed

        In the same run, the code extract contig sizes from the resulting Fasta file
        """
        ## unmask the genome
        self._to_log("Replacing ambiguous symbols with adenosines in the genome file")
        unmask_cmd: str = (
            "set -eu; set -o pipefail; "
            f"{self.twobittofa_binary} {self.twobit} stdout | "
            f"sed '/^>/!s/[BD-FH-SU-Z]/A/Ig' | tee {self.tmp_fa} |"
            f"{CONTIG_SIZE_SCRIPT} - -o {self.chrom_sizes}"
        )
        _ = self._exec(unmask_cmd, "Genome unmasking & contig size retrieval failed:")
        ## convert the resulting fasta back into 2bit
        self._to_log("Compressing the modified genome sequence")
        compress_cmd: str = (
            f"{self.fatotwobit_binary} {self.tmp_fa} {self.unmasked_twobit}"
        )
        _ = self._exec(compress_cmd, "faToTwoBit compression failed:")

    def schedule_jobs(self) -> None:
        """
        Schedules SpliceAI jobs for cluster execution.
        Individual commands are grouped in {self.n_jobs} job bins. Jobs are equilibrated
        using longest-processing-time-first (LPT) algorithm with sequence lengths
        as proxies for execution time per command
        """
        bin2chunks: Dict[str, List[Tuple[str, int, int, str]]] = defaultdict(list)
        job_heap: List[Tuple[int, int]] = [(0, i) for i in range(self.job_num)]
        chunk_num: int = 0
        with open(self.chrom_sizes, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data:
                    continue
                if len(data) != 2:
                    self._die(
                        (
                            "Improper formatting in the chromosome size file at line %s; "
                            "expected 2 columns, got % i"
                        )
                        % (i, len(data))
                    )
                contig, length = data
                length = int(length)
                ## ignore contigs shorter than the minimal threshold
                if length < self.min_contig_size:
                    continue
                ## each contig is split into chunks of self.chunk_size length each + 2 * self.flank_size
                if length >= self.chunk_size:
                    ## longer sequences are split into individual chunks
                    chunk_num: int = length // self.chunk_size
                    if chunk_num * self.chunk_size < length:
                        chunk_num += 1
                    for s in range(0, chunk_num):
                        # for s in range(0, length, self.chunk_size):
                        start: int = s * self.chunk_size
                        flanked_start: int = max(0, start - self.flank_size)
                        end: int = min(length, start + self.chunk_size)
                        flanked_end: int = min(length, end + self.flank_size)
                        chunk_size: str = end - start
                        lightest_bin: Tuple[int, int] = heappop(job_heap)
                        total_bin_size, bin_id = lightest_bin
                        # name: str = f'{contig}:{start}-{end}'
                        name: str = f"chunk{chunk_num}"
                        bin2chunks[bin_id].append(
                            "\t".join(
                                map(
                                    str,
                                    (
                                        contig,
                                        flanked_start,
                                        flanked_end,
                                        name + "{}",
                                        0,
                                        "{}",
                                        start,
                                        end,
                                    ),
                                )
                            )
                        )
                        total_bin_size += chunk_size
                        heappush(job_heap, (total_bin_size, bin_id))
                        chunk_num += 1
                    ## if there is a trailing terminal portion,
                    ## flank it from the left side and push it as well
                    if self.chunk_size * chunk_num < length:
                        start: int = chunk_num * self.chunk_size
                        flanked_start: int = max(0, start - self.flank_size)
                        end: int = length
                        chunk_size = end - start
                        lightest_bin: Tuple[int, int] = heappop(job_heap)
                        total_bin_size, bin_id = lightest_bin
                        # name: str = f'{contig}:{flanked_start}-{end}'
                        # bin2chunks[bin_id].append((contig, start, end))
                        name: str = f"chunk{chunk_num}"
                        bin2chunks[bin_id].append(
                            "\t".join(
                                map(
                                    str,
                                    (
                                        contig,
                                        flanked_start,
                                        flanked_end,
                                        name + "{}",
                                        0,
                                        "{}",
                                        start,
                                        end,
                                    ),
                                )
                            )
                        )
                        total_bin_size += chunk_size
                        heappush(job_heap, (total_bin_size, bin_id))
                        chunk_num += 1
                else:
                    ## just push the whole chunk
                    lightest_bin: Tuple[int, int] = heappop(job_heap)
                    total_bin_size, bin_id = lightest_bin
                    # name: str = f'{contig}:{0}-{length}'
                    # bin2chunks[bin_id].append(
                    #     (contig, 0, length, name)
                    # )
                    name: str = f"chunk{chunk_num}"
                    bin2chunks[bin_id].append(
                        "\t".join(
                            map(
                                str,
                                (contig, 0, length, name + "{}", 0, "{}", 0, length),
                            )
                        )
                    )
                    total_bin_size += length
                    heappush(job_heap, (total_bin_size, bin_id))
                    chunk_num += 1
        ## remove the hefty list
        del job_heap
        ## write the resulting jobs
        job_list: List[str] = []
        ## write input for each bin in Bed6 format
        for bucket, intervals in bin2chunks.items():
            if not intervals:
                continue
            batch_prefix: str = f"batch{self.chunk_num}"
            bed_file: str = os.path.join(self.bed_dir, f"{batch_prefix}.bed")
            with open(bed_file, "w") as h:
                for interval in intervals:
                    ## record each chunk twice, once for each strand
                    for strand in STRANDS:
                        # upd_name: str = f'{name}{strand}'
                        # upd_interval: Tuple[Union[str, int]] = (
                        #     chrom, start, end, upd_name, 0, strand
                        # )
                        # line: str = '\t'.join(map(str, upd_interval))
                        line: str = interval.format(strand, strand)
                        h.write(line + "\n")
            ## add the resulting command to the job list
            cmd: str = (
                f"{EXEC_SCRIPT} {self.unmasked_twobit} {bed_file} "
                f"--round_to {self.round_to} --min_prob {self.min_prob} "
                f"-o {self.tmp_dir} -b {batch_prefix} "
                f"--twobittofa_binary {self.twobittofa_binary} "
                f"--wigtobigwig_binary {self.wigtobigwig_binary}"
            )
            job_list.append(cmd)
            self.chunk_num += 1
        with open(self.job_list, "w") as h:
            for job in job_list:
                h.write(job + "\n")

    def run_jobs(self) -> None:
        """
        Controls parallel job execution.
        This is a medley of parallel execution-related methods from toga_main.py
        simplified for the needs of the SpliceAI annotaiton mode.
        """
        ## get parallel process manager based on the requested strategy
        if self.parallel_strategy == "para":
            strategy = ParaStrategy()
        elif self.parallel_strategy == "custom":
            strategy = CustomStrategy()
        else:
            strategy = NextflowStrategy()
        job_manager: ParallelJobsManager = ParallelJobsManager(strategy)
        local_executor: bool = self.parallel_strategy == "local"

        ## generate Nextflow configuration fil
        if self.nextflow_config_file is None:
            nf_contents: str = Constants.NEXTFLOW_STUB.format(
                self.max_number_of_retries
            )
            nf_file: str = os.path.join(self.nextflow_dir, "execute_joblist.nf")
            self.nextflow_exec_script = nf_file
            with open(nf_file, "w") as h:
                h.write(nf_contents + "\n")

        project_name: str = self.project_name
        project_path: str = os.path.join(self.nextflow_dir, project_name)
        manager_data: Dict[str, str] = {
            "project_name": project_name,
            "project_path": project_path,
            "logs_dir": project_path,
            "nextflow_dir": self.nextflow_dir,
            "NF_EXECUTE": self.nextflow_exec_script,
            "local_executor": local_executor,
            "keep_nf_logs": self.keep_tmp,
            # 'nexflow_config_file': nextflow_config,
            "nextflow_config_dir": self.nextflow_dir,
            "temp_wd": self.tmp_dir,
            "queue_name": self.cluster_queue_name,
            "logger": self.logger,
        }
        if self.nextflow_config_file is not None:
            manager_data["nexflow_config_file"] = self.nextflow_config_file
        self._to_log("Launching the parallel SpliceAI jobs")
        try:
            job_manager.execute_jobs(
                self.job_list,
                manager_data,
                project_name,
                wait=True,
                memory_limit=self.memory_limit,  ## TODO: Potentially should be set to custom values
                queue_name=self.cluster_queue_name,
                clean=self.keep_tmp,
                cpu=1,
                process_num=self.chunk_num + 1,
                executor=self.parallel_strategy,
            )
        except KeyboardInterrupt:
            # TogaUtil.terminate_parallel_processes([job_manager])
            self._to_log("Aborting the parallel step", "warning")
            job_manager.terminate_process()

    def aggregate_jobs(self) -> None:
        """
        Aggregates individual job results and converts them into into bigWig format,
        one per each splice site and each strand
        """
        self._to_log("Aggregating SpliceAI prediction results")
        for template in FILE_NAME_TEMPLATES:
            final_file: str = os.path.join(self.output, f"spliceAi{template}.bw")
            tmp_aggr_wig: str = os.path.join(self.tmp_dir, f"spliceAi{template}.wig")
            files_to_aggr: List[str] = [
                x for x in os.listdir(self.tmp_dir) if template in x and "spliceAi" not in x
            ]
            # files_to_aggr.sort(
            #     key=lambda x: int(x.replace(template, "").replace("batch", "").replace(".wig", ""))
            # )
            # with open(tmp_aggr_wig, "w"):
            #     pass
            # for file in files_to_aggr:
            #     filepath: str = os.path.join(self.tmp_dir, file)
            #     if "batch0" in file:
            #         add_cmd: str = f"cat {filepath} > {tmp_aggr_wig}"
            #     else:
            #         add_cmd: str = f"cat {filepath} >> {tmp_aggr_wig}"
                # add_cmd = f"grep -v fixedStep {filepath} >> {tmp_aggr_wig}"
                # _ = self._exec(add_cmd, f"Adding file {filepath} to the main Wiggle stub failed")
                        ## sort the files by chromosome
            chrom2output: Dict[str, List[str]] = defaultdict(list)
            for file in files_to_aggr:
                chrom_name: str = "@".join(file.split("@")[:-1])
                chrom2output[chrom_name].append(file)
            for chrom, chrom_files_to_aggr in chrom2output.items():
                chrom_files_to_aggr.sort(key=lambda x: batch_num(x, template))
                for file in chrom_files_to_aggr:
                    filepath: str = os.path.join(self.tmp_dir, file)
                    if not os.path.exists(tmp_aggr_wig) or os.stat(tmp_aggr_wig).st_size == 0:
                        add_cmd: str = f"cat {filepath} > {tmp_aggr_wig}"
                    else:
                        add_cmd: str = f"cat {filepath} >> {tmp_aggr_wig}"
                    _ = self._exec(add_cmd, f"Adding file {filepath} to the main Wiggle stub failed")
            convert_cmd: str = f"{self.wigtobigwig_binary} {tmp_aggr_wig} {self.chrom_sizes} {final_file}"
            _ = self._exec(
                convert_cmd, 
                "Wiggle to BigWig convertion failed for file %s" % tmp_aggr_wig,
            )

    def _execute_step(self, step: str):
        """
        Defines whether the current step is to be executed based on 'resume' and 'halt' options
        """
        if step not in RESUME_ORDER:
            self._die("Improper step name provided")
        regular_start: bool = self.resume_from == "all"
        regular_finish: bool = self.halt_at == "all"
        step: int = RESUME_ORDER[step]
        resume_step: int = RESUME_ORDER[self.resume_from]
        halt_step: int = RESUME_ORDER[self.halt_at]
        exec_started: bool = resume_step <= step or regular_start
        exec_not_finished: bool = step < halt_step or regular_finish
        return exec_started and exec_not_finished
