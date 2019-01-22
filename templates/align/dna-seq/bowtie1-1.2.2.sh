#! /bin/sh
# Auto-generated sequence alignment script.

# auto-populated script-specific fields
script_name='<script.name>'
script_tag='<script.tag>'
genome_assembly='<alignment.genome.assembly>'

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

# assert that `bowtie` is available
bowtie1_version='1.2.2'
PATH="${SOFTWARE:-/opt}/bowtie1-${bowtie1_version}/bin${PATH:+:${PATH}}"
bowtie1="$(command -v bowtie)"
[ -z "${bowtie1}" ] &&
	fatal 2 $(command printf 'Command `bowtie` not in PATH: %s' "${PATH}")
{
	"${bowtie1}" --version 2>&1 |
		command grep --quiet "bowtie-align version ${bowtie1_version}"
} 1>/dev/null 2>&1 ||
	fatal 2 $(command printf 'Incompatible version of `bowtie` (want %s): %s' "${bowtie1_version}" "${bowtie1}")

# assert that the target project is available
project_dir="${PROJECTS:-${HOME}/Projects}/${project_name}"
[ -d "${project_dir}" ] ||
	fatal 2 $(command printf 'Project directory does not exist: %s' "${project_dir}")
[ -r "${project_dir}" ] ||
	fatal 13 $(command printf 'Cannot read from project directory (permission denied): %s' "${project_dir}")
[ -w "${project_dir}" ] ||
	fatal 13 $(command printf 'Cannot write to project directory (permission denied): %s' "${project_dir}")

# assert that the project data directory is available
data_dir="Data/${sample_datatype}/${sample_subtype}/${sample_name}/${sample_rep}/${lib_rep}/${sample_rep}"
project_data_dir="${project_dir}/${data_dir}"
[ -d "${project_data_dir}" ] ||
	fatal 2 $(command printf 'Project data directory does not exist: %s' "${project_data_dir}")
[ -r "${project_data_dir}" ] ||
	fatal 13 $(command printf 'Cannot read from project data directory (permission denied): %s' "${project_data_dir}")
[ -w "${project_data_dir}" ]||
	fatal 13 $(command printf 'Cannot write to project data directory (permission denied): %s' "${project_data_dir}")

# assert that the raw reads are available
reads="${project_data_dir}/reads.fastq.gz"
mate1s="${project_data_dir}/mate1s.fastq.gz"
mate2s="${project_data_dir}/mate2s.fastq.gz"
if [ -f "${reads}" ]; then
	[ -r "${reads}" ] ||
		fatal 13 $(command printf 'Cannot read single-end reads file (permission denied): %s' "${reads}")
	set -- "${reads}"
elif [ -f "${mate1s}" ] && [ -f "${mate2s}" ]; then
	[ -r "${mate1s}" ] ||
		fatal 13 $(command printf 'Cannot read paired-end mate 1 reads file (permission denied): %s' "${mate1s}")
	[ -r "${mate2s}" ] ||
		fatal 13 $(command printf 'Cannot read paired-end mate 2 reads file (permission denied): %s' "${mate2s}")
	[ -z "${lib_fragsize}" ] &&
		fatal 22 $(command printf 'Mean fragment size required for pair-end reads.')
	[ "${lib_fragsize}" -gt 0 ] 1>/dev/null 2>&1 ||
		fatal 22 $(command printf 'Mean fragment size must be a positive number, got: %s' "${lib_fragsize}")
	set -- '-1' "${mate1s}" '-2' "${mate2s}" \
		'--minins' "$(( lib_fragsize - (lib_fragsize / 2) ))" \
		'--maxins' "$(( lib_fragsize + (lib_fragsize / 2) ))"
else
	fatal 2 $(command printf 'Missing sequencing reads file(s) for sample: %s' "${project_data_dir}")
fi

# assert that the Bowtie 1 genome-specific files are available
genome_dir="${DATA:-${HOME}/Data}/Genomes/${genome_assembly}/bowtie1-${bowtie1_version}"
[ -d "${genome_dir}" ] ||
	fatal 2 $(command printf 'Missing Bowtie 1 genome directory: %s' "${genome_dir}")
input_index="${genome_dir}/genome"
{ ext='ebwtl' && [ -f "${input_index}.1.${ext}" ]; } ||
	{ ext='ebwt' && [ -f "${input_index}.1.${ext}" ]; } ||
		fatal 2 $(command printf 'Missing Bowtie 1 reference genome index files: %s' "${input_index}.*.${ext}")
for index_file in "${input_index}."{1,2,3,4}".${ext}" "${input_index}.rev."{1,2}".${ext}"; do
	[ -f "${index_file}" ] ||
		fatal 2 $(command printf 'Missing Bowtie 1 reference genome index file: %s' "${index_file}")
	[ -r "${index_file}" ] ||
		fatal 13 $(command printf 'Cannot read Bowtie 1 reference genome index file (premission denied): %s' "${index_file}")
done
unset ext

# assert that the target output is available
output_sam="${project_data_dir}/alignment.sam.${script_tag}"
[ -f "${output_sam}" ] &&
	fatal 17 $(command printf 'Output SAM file already exists: %s' "${output_sam}")

# (finally) perform alignment
export PATH
"${bowtie1}" \
	"${input_index}" "${@}" \
	--sam "${output_sam}" \
		--sam-RG "ID:${sample_datatype}.${sample_subtype}.${sample_name}.${sample_rep}.${lib_rep}.${seq_rep}" \
		--sam-RG "SM:${sample_subtype}.${sample_name}.${sample_rep}" \
		--sam-RG "LB:${sample_subtype}.${sample_name}.${sample_rep}.${lib_rep}" \
		--sam-RG "BC:${lib_barcode}" \
		--sam-RG "PI:${lib_fragsize}" \
		--sam-RG "DT:${seq_date}" \
		--sam-RG "PL:${seq_type}" \
		--sam-RG "PM:${seq_mach}" \
		--sam-RG "PU:${seq_run}.${seq_lane}.${lib_barcode}" \
	--threads "${NSLOTS:-1}" --best \
	--time --seed "${random_seed:-${RANDOM}}"

# link the outputs of this instance as the new default outputs
for loop_path in "${output_sam}"; do
	loop_dir="$(command dirname -- "${loop_path}")"
	loop_file="$(command basename -- "${loop_path}" "${script_tag}")"
	command ln --symbolic --force -- "${loop_file}" "${loop_dir}/${loop_file}"
done
