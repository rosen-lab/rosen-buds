#! /bin/sh
# Auto-generated sequence alignment script.

# auto-populated script-specific fields
script_name='<script.name>'
script_tag='<script.tag>'
genome_assembly='<alignment.genome.assembly>'
genome_annotation='<alignment.genome.annotation>'

# auto-populated sample-constant fields
project_name='<project.name>'
sample_datatype='<sample.datatype>'
sample_subtype='<sample.subtype>'
sample_name='<sample.name>'
sample_rep='<sample.replicate>' # sample number amongst all biological/sample replicates
lib_barcode='<sample.library.barcode>' # sample barcode(s) (combine with hyphen)
lib_fragsize='<sample.library.fragment-size>' # median fragment size
lib_rep='<sample.library.replicate>' # number amongst all library replicates
seq_run='<run.flowcell>' # run/flowcell id
seq_date='<run.date>' # date of sequencing run
seq_mask='<run.mask>' # sequencing run mask
seq_lane='<sample.sequencing.lanes>' # lanes sample was loaded onto
seq_mach='<run.sequencer>' # sequencing machine
seq_type='<sequencer.platform>' # ILLUMINA, PACBIO, etc.
seq_rep='<sample.sequencing.replicate>' # number amongst all sequencing replicates

# default values for auto-populated fields
script_name="${script_name:-${SCRIPT_NAME:-$(command basename -- "${0}" '.sh')}}"
script_tag="${script_tag:-$(command date '%s')}"
genome_annotation="${genome_annotation:-alyubets}"

# source rosen-buds shell utilities
. "${SOFTWARE:-/srv/git}/rosen-buds/sh/functions/log-simple.sh"

# warn that arguments to this script are purposefully ignored
[ "${#}" != 0 ] &&
	warning $(command printf 'Ignoring argument(s): %s' "${*}")
set -- # clear the positional parameters, just in case

# assert that `hisat2` is available
hisat2_version='2.0.5'
PATH="${SOFTWARE:-/opt}/hisat2-${hisat2_version}/bin${PATH:+:${PATH}}"
hisat2="$(command -v hisat2)"
[ -z "${hisat2}" ] &&
	fatal 2 $(command printf 'Command `hisat2` not in PATH: %s' "${PATH}")
{
	"${hisat2}" --version 2>&1 |
		command grep --quiet "hisat2-align-s version ${hisat2_version}"
} 1>/dev/null 2>&1 ||
	fatal 2 $(command printf 'Incompatible version of `hisat2` (want %s): %s' "${hisat2_version}" "${hisat2}")

# assert that the target project is available
project_dir="${PROJECTS:-${HOME}/Projects}/${project_name}"
[ -d "${project_dir}" ] ||
	fatal 2 $(command printf 'Project directory does not exist: %s' "${project_dir}")
[ -r "${project_dir}" ] ||
	fatal 13 $(command printf 'Cannot read from project directory (permission denied): %s' "${project_dir}")
[ -w "${project_dir}" ] ||
	fatal 13 $(command printf 'Cannot write to project directory (permission denied): %s' "${project_dir}")

# assert that the project data directory is available
data_dir="Data/${sample_datatype}/${sample_subtype}/${sample_name}/${sample_rep}/${lib_rep}/${seq_rep}"
project_data_dir="${project_dir}/${data_dir}"
[ -d "${project_data_dir}" ] ||
	fatal 2 $(command printf 'Project data directory does not exist: %s' "${project_data_dir}")
[ -r "${project_data_dir}" ] ||
	fatal 13 $(command printf 'Cannot read from project data directory (permission denied): %s' "${project_data_dir}")
[ -w "${project_data_dir}" ] ||
	fatal 13 $(command printf 'Cannot write to project data directory (permission denied): %s' "${project_data_dir}")

# assert that the raw reads are available
reads="${project_data_dir}/reads.fastq.gz"
mate1s="${project_data_dir}/mate1s.fastq.gz"
mate2s="${project_data_dir}/mate2s.fastq.gz"
if [ -f "${reads}" ]; then
	[ -r "${reads}" ] ||
		fatal 13 $(command printf 'Cannot read single-end reads file (permission denied): %s' "${reads}")
	set -- '-U' "${reads}"
elif [ -f "${mate1s}" ] && [ -f "${mate2s}" ]; then
	[ -r "${mate1s}" ] ||
		fatal 13 $(command printf 'Cannot read paired-end mate 1 reads file (permission denied): %s' "${mate1s}")
	[ -r "${mate2s}" ] ||
		fatal 13 $(command printf 'Cannot read paired-end mate 2 reads file (permission denied): %s' "${mate2s}")
	[ -z "${lib_fragsize}" ] &&
		fatal 22 $(command printf 'Mean fragment size required for pair-end reads.')
	[ "${lib_fragsize}" -gt 0 ] 1>/dev/null 2>&1 ||
		fatal 22 $(command printf 'Mean fragment size must be a positive integer, got: %s' "${lib_fragsize}")
	set -- '-1' "${mate1s}" '-2' "${mate2s}" \
		'--minins' "$(( lib_fragsize - (lib_fragsize / 2) ))" \
		'--maxins' "$(( lib_fragsize + (lib_fragsize / 2) ))"
else
	fatal 2 $(command printf 'Missing sequencing reads file(s) for sample: %s' "${project_data_dir}")
fi

# assert that the HISAT2 genome-specific files are available
genome_dir="${DATA:-${HOME}/Data}/Genomes/${genome_assembly}/hisat2-${hisat2_version}"
[ -d "${genome_dir}" ] ||
	fatal 2 $(command printf 'Missing HISAT2 genome directory: %s' "${genome_dir}")
input_index="${genome_dir}/${genome_annotation}"
{ ext='ht2l' && [ -f "${input_index}.1.${ext}" ]; } ||
	{ ext='ht2' && [ -f "${input_index}.1.${ext}" ]; } ||
		fatal 2 $(command printf 'Missing HISAT 2 reference genome index file: %s' "${input_index}.*.${ext}")
for index_file in "${input_index}."{1,2,3,4,5,6,7,8}".${ext}"; do
	[ -f "${index_file}" ] ||
		fatal 2 $(command printf 'Missing HISAT 2 reference genome index file: %s' "${index_file}")
	[ -r "${index_file}" ] ||
		fatal 13 $(command printf 'Cannot read HISAT 2 reference genome index file (premission denied): %s' "${index_file}")
done
input_splices="${genome_dir}/${genome_annotation}.splices"
[ -f "${input_splices}" ] ||
	fatal 2 $(command printf 'Missing HISAT2 known splices file: %s' "${input_splices}")
[ -r "${input_splices}" ] ||
	fatal 2 $(command printf 'Cannot read HISAT2 known splices file (permission denied): %s' "${input_splices}")

# assert that the target output is available
output_sam="${project_data_dir}/alignment.sam.${script_tag}"
[ -f "${output_sam}" ] &&
	fatal 17 $(command printf 'Output SAM file already exists: %s' "${output_sam}")

# assert that the temporary (incidental) output is available
temp_dir="${project_data_dir}/Temporary"
[ -d "${temp_dir}" ] ||
	fatal 2 $(command printf 'Temporary directory does not exist: %s' "${temp_dir}")
[ -r "${temp_dir}" ] ||
	fatal 13 $(command printf 'Cannot read from temporary directory (permission denied): %s' "${temp_dir}")
[ -w "${temp_dir}" ] ||
	fatal 13 $(command printf 'Cannot write to project data directory (permission denied): %s' "${temp_dir}")
temp_output_metrics="${temp_dir}/${script_name}.metrics.${script_tag}"
[ -f "${temp_output_metrics}" ] && # this *should* never succeed
	fatal 17 $(command printf 'Metrics file already exists: %s' "${temp_output_metrics}")
temp_output_splices="${temp_dir}/${script_name}.splices.${script_tag}"
[ -f "${temp_output_splices}" ] && # this *should* never succeed
	fatal 17 $(command printf 'Novel splices file already exists: %s' "${temp_output_metrics}")

# align reads
export PATH
"${hisat2}" \
	"${@}" \
	-x "${input_index}" \
	-S "${output_sam}" \
		--rg-id "${sample_datatype}.${sample_subtype}.${sample_name}.${sample_rep}.${lib_rep}.${seq_rep}" \
		--rg "SM:${sample_subtype}.${sample_name}.${sample_rep}" \
		--rg "LB:${sample_subtype}.${sample_name}.${sample_rep}.${lib_rep}" \
		--rg "BC:${lib_barcode}" \
		--rg "PI:${lib_fragsize}" \
		--rg "DT:${seq_date}" \
		--rg "PL:${seq_type}" \
		--rg "PM:${seq_mach}" \
		--rg "PU:${seq_run}.${seq_lane}.${lib_barcode}" \
	--known-splicesite-infile "${input_splices}" \
	--novel-splicesite-outfile "${temp_output_splices}" \
	--met "$(( 60 ))" \
		--met-stderr \
		--met-file "${temp_output_metrics}" \
	--threads "${NSLOTS:-1}" \
	--time --seed "${random_seed:-${RANDOM}}"

# link the outputs of this instance as the new default outputs
for loop_path in "${output_sam}"; do
	loop_dir="$(command dirname -- "${loop_path}")"
	loop_file="$(command basename -- "${loop_path}" "${script_tag}")"
	command ln --symbolic --force "${loop_file}" "${loop_dir}/${loop_file}"
done
