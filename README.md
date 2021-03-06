YAMA
====

Yet Another Methylome Aligner

Copyright (C) 2013 Thomas J. Hardcastle <tjh48@cam.ac.uk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  You should have received a copy of 
the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

What is YAMA?

YAMA is a *-nix based perl script that performs unbiased alignment and calling of bisulphite treated sequencing (BS-Seq) data. It is released under the GPL-3 licence.

Yama is the lord of death as described in the Vedas, a religously inclined llama, a "largely forgettable" 1979 album of flugelhorn music, and a security module for the Linux kernel. If you are looking for these, or any other Yama, you are probably in the wrong place.

Why is YAMA?

YAMA is designed to produce output optimised for downstream visualisation, methylation locus detection, and differential methylation analyses using the R package segmentSeq. For a general introduction to processing of BS-Seq data, see Hardcastle (2013), Plant Methods, 9:16 (doi:10.1186/1746-4811-9-16).

Principal features of YAMA:

YAMA performs fast unbiased alignment (using C-T conversion) from FASTQ files
Alignment with either Bowtie or Bowtie2 aligner is supported.
(Almost) all options of Bowtie(2) may be specified in the alignment, allowing for gapped alignment, flexibility of seed length, choice of mismatches, et cetera.
Paired-end (Bowtie2 only) and single-end data are supported.
The output produced may be context (CpG, CHG and CHH) dependent (default) or independent.
Alignment of multireads is permitted, with the base-calling being weighted by the number of locations to which the multiread aligns.



Example input (single-end reads)

perl yama.pl -x referenceGenome.fasta -Q sequencedRead.fastq -o myDir -n myOutputs --bowtie2="-p 16 -k 50"

Example input (paired-end reads)

perl yama.pl -x referenceGenome.fasta -1 pair1.fastq -2 pair2.fastq -o myDir -n myOutputs --bowtie2="-p 16 -X 1500 -I 300 -k 50 --no-mixed --no-discordant"

Essential parameters

-P   a fastq file containing the bisulphite converted sequenced reads

-D   the reference genome on which you wish to identify methylation

-n   prefix to output files

Optional parameters/flags

-Q A comma-separated list of fastq files.

-1 A comma-separated list of fastq files, paired with -2.

-2 A comma-separated list of fastq files, paired with -1.

--useBowtie - use bowtie rather than bowtie2 aligner.

--bowtie2 - options for bowtie2 aligner - consult bowtie2 manual for details.

--bowtie - options for bowtie aligner - consult bowtie manual for details. If bowtie options are provided and bowtie2 options are not, the bowtie aligner will be used.

-o - output directory. Defaults to current working directory.

-i - the second element of the pair will not be inverted before C2T conversion. This will be needed in the rare case where the paired end sequencing is of the form 'forward-forward' rather than the usual 'forward-reverse'

-r - reuse temporary files from a previous run. YAMA will always try and do this for the genome converted file, but will by default ignore temporary files created for previous alignments. Use with caution - if a run has previously aborted, the files are likely to be incomplete.

--no-split - don't split the final output files by methylation context. Defaults to off.

-t - name of the temporary directory to be used by YAMA. If not given, defaults to 'myDir/yama_tmp_myOutputs'

-p - the number of parallel processors used by the 'sort' system call. If not specified, will steal the number from the bowtie or bowtie2 options if it is given here.

-S - the memory limit for the 'sort' system call. Defaults to 50% of available memory; see the man pages of 'sort' for options.

--keep-all-temp - keep all temporary files produced by YAMA. These are likely to be large, and so use of this is generally discouraged unless debugging.

--recover-from-temp - recover a crashed process from the last set of temporary files.

YAMA will create a conversion of the reference genome in the same location as the reference genome file called:
referenceGenome.fasta_methConvert (and associated bowtie2-index files)
These are the C-to-T (and G-to-A) reference files to which the fastq files will be aligned. If you do not have permission to write to this directory, YAMA will fail - copy the reference genome to some location where you have permission to write.

A set of output files will be created by YAMA in the temporary directory. These are

myOutputs_bestalignments	  Processed alignment files - mostly for debugging purposes

myOutputs_basecount	  Tab-delimited text file defining each cytosine base hit by a read and the identified methylation/non-methylation count due to that read, and the number of times that read maps to the genome.

myOutputs_basecount_duplication		Data on reads that map multiple times to the genome. Currently unused by any downstream process, but may be useful later.

myOutputs_basecount_sort			Identical to the _basecount file, but sorted.

myOutputs_basecount_sort_duplication	Identical to the _basecount_duplication file, but sorted.

The output files are (by default)

myOutputs_CG_methCalls
myOutputs_CHG_methCalls
myOutputs_CHH_methCalls


or if (--no-split) is used

myOutputs_methCalls
