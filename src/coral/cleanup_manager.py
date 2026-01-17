import os
import shutil
from .utils import log

class PipelineCleaner:
    def __init__(self, genomes, alignments, pileup, base_dir=None, verbose=True):
        self.genomes = genomes
        self.alignments = alignments
        self.pileup = pileup
        self.base_dir = base_dir
        self.verbose = verbose

    def _log(self, message):
        log(message, verbose=self.verbose)

    def _safe_rm(self, path):
        if os.path.isfile(path):
            os.remove(path)
            self._log(f"Removed file: {path}")
        elif os.path.isdir(path):
            shutil.rmtree(path)
            self._log(f"Removed directory: {path}")

    def clean_bams_folder(self):
        path = self.alignments[0].bam_dir
        if os.path.exists(path):
            self._log("Removing BAM folder...")
            self._safe_rm(path)

    def clean_raw_bams(self):
        self._log("Removing raw BAM files...")
        for alignment in self.alignments:
            if os.path.exists(alignment.raw_bam):
                self._safe_rm(alignment.raw_bam)

    def clean_pileup(self):
        if os.path.exists(self.pileup.pileup_path):
            self._log("Removing pileup file...")
            self._safe_rm(self.pileup.pileup_path)

    def clean_intervals(self):
        if not self.base_dir:
            self._log("Warning: base_dir not provided, skipping interval cleanup.")
            return
        path = os.path.join(self.base_dir, "Intervals")
        if os.path.exists(path):
            self._log("Removing interval coverage files...")
            self._safe_rm(path)

    def clean_genome_folders(self):
        self._log("Removing genome FASTA folders...")
        for genome in self.genomes + []:
            if os.path.exists(genome.output_dir):
                self._safe_rm(genome.output_dir)

    def clean_fasta(self):
        self._log("Removing FASTA files...")
        for genome in self.genomes:
            if os.path.exists(genome.fasta_path):
                self._safe_rm(genome.fasta_path)

    def clean_fastq(self):
        self._log("Removing fragmentation FASTQs...")
        for genome in self.genomes:
            if os.path.exists(genome.fastq_path):
                self._safe_rm(genome.fastq_path)

    def run(self, bams=False, pileup=False, intervals=False, genomes=False, fastas=False, fastqs=False, bam_folder=False):
        if bam_folder:
            self.clean_bams_folder()
        else:
            if bams:
                self.clean_raw_bams()
        if pileup:
            self.clean_pileup()
        if intervals:
            self.clean_intervals()
        if genomes:
            self.clean_genome_folders()
        else:
            if fastas:
                self.clean_fasta()
            if fastqs:
                self.clean_fastq()
        self._log("Cleanup complete.")

