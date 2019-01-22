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

# assert that `STAR` is available
star_version='2.6.1d'
PATH="${SOFTWARE:-/opt}/star-${star_version}/bin${PATH:+:${PATH}}"
star="$(command -v STAR)"
[ -z "${star}" ] &&
	fatal 2 $(command printf 'Command `STAR` not in PATH: %s' "${PATH}")
{
	"${star}" --version 2>&1 |
		command grep --quiet "STAR_${star_version}"
} 1>/dev/null 2>&1 ||
	fatal 2 $(command printf 'Incompatible version of `STAR` (want %s): %s' "${star_version}" "${star}")

# assert that zutils `zcat` is available (efficient multi-format decompressor)
zutils_version='1.4'
decompress="$(command -v zcat)"
[ -z "${decompress}" ] &&
	fatal 2 $(command printf 'Command `zcat` not in PATH: %s' "${PATH}")
decompress_version="$("${decompress}" --version 2>&1 | command grep --only-matching 'zcat (zutils) [0-9][0-9]*.[0-9][0-9]*')"
[ -z "${decompress_version}" ] &&
	fatal 2 $(command printf 'Unrecognized version of `zcat` (want zutils: %s): %s' "${zutils_version}" "${PATH}")
(command printf '%s' "${decompress_version}" | command grep --quiet "zcat (zutils) ${zutils_version}") ||
	warning $(command printf 'Unrecognized `zcat` version (prefer: %s): %s' "${zutils_version}" "${decompress}")

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
	set -- "${reads}"
elif [ -f "${mate1s}" ] && [ -f "${mate2s}" ]; then
	[ -r "${mate1s}" ] ||
		fatal 13 $(command printf 'Cannot read paired-end mate 1 reads file (permission denied): %s' "${mate1s}")
	[ -r "${mate2s}" ] ||
		fatal 13 $(command printf 'Cannot read paired-end mate 2 reads file (permission denied): %s' "${mate2s}")
	set -- "${mate1s}" "${mate2s}"
else
	fatal 2 $(command printf 'Missing sequencing reads file(s) for sample: %s' "${project_data_dir}")
fi

# assert that the STAR genome-specific files are available
genome_dir="${DATA:-${HOME}/Data}/Genomes/${genome_assembly}/star-${star_version}"
[ -d "${genome_dir}" ] ||
	fatal 2 $(command printf 'Missing STAR genome directory: %s' "${genome_dir}")
[ -r "${genome_dir}" ] ||
	fatal 13 $(command printf 'Cannot read from genome directory (permission denied): %s' "${genome_dir}")
input_sjdb="${genome_dir}/${genome_annotation:-alyubets}.sjdb"
[ -f "${input_sjdb}" ] ||
	fatal 2 $(command printf 'Missing SJDB annotation file: %s' "${input_sjdb}")
[ -r "${input_sjdb}" ] ||
	fatal 13 $(command printf 'Cannot read SJDB annotation file (permission denied): %s' "${input_sjdb}")
[ "${library_readsize}" -gt 0 ] 1>/dev/null 2>&1 ||
	fatal 22 $(command printf 'Read size must be a positive number, got: %s' "${library_readsize}")

# assert that the target output is available
output_dir="${project_data_dir}/alignment.${script_tag}"
[ -d "${output_dir}" ] &&
	fatal 17 $(command printf 'Output directory already exists: %s' "${output_dir}")
command mkdir -- "${output_dir}"

# assert that the temporary (incidental) output is available
temp_dir="${project_data_dir}/Temporary"
[ -d "${temp_dir}" ] ||
	fatal 2 $(command printf 'Temporary directory does not exist: %s' "${temp_dir}")
[ -r "${temp_dir}" ] ||
	fatal 13 $(command printf 'Cannot read from temporary directory (permission denied): %s' "${temp_dir}")
[ -w "${temp_dir}" ] ||
	fatal 13 $(command printf 'Cannot write to project data directory (permission denied): %s' "${temp_dir}")
temp_output_stats="${temp_dir}/${script_name}.stats.${script_tag}"
[ -f "${temp_output_stats}" ] && # this *should* never succeed
	fatal 17 $(command printf 'Statistics file already exists: %s' "${temp_output_stats}")
temp_output_sjs="${temp_dir}/${script_name}.sjs.${script_tag}"
[ -f "${temp_output_sjs}" ] && # this *should* never succeed
	fatal 17 $(command printf 'Splice junctions file already exists: %s' "${temp_output_sjs}")

# align reads
export PATH
"${star}" \
	--readFilesIn "${@}" \
		--readFilesCommand "${decompress}" \
	--genomeDir "${genome_dir}" \
		--genomeLoad 'NoSharedMemory' \
	--sjdbGTFfile "${input_sjdb}" \
		--sjdbGTFchrPrefix 'chr' \
		--sjdbGTFfeatureExon 'exon' \
		--sjdbGTFtagExonParentTranscript 'transcript_id' \
		--sjdbGTFtagExonParentGene 'gene_id' \
		--sjdbOverhang "$(( $(command printf '%s' "${seq_mask}" | command sed 's/^.*Y\([0-9][0-9]*\).*$/\1/') - 1 ))" \
	--outFileNamePrefix "${output_dir}/" \
		--outSAMtype 'SAM' \
			--outSAMmode 'Full' \
			--outSAMorder 'Paired' \
			--outSAMattributes 'AS' 'NM' 'nM' 'MD' 'XS' 'NH' \
			--outSAMattrRGline "ID:${sample_datatype}.${sample_subtype}.${sample_name}.${sample_rep}.${lib_rep}.${seq_rep}" \
				"SM:${sample_subtype}.${sample_name}.${sample_rep}" \
				"LB:${sample_subtype}.${sample_name}.${sample_rep}.${lib_rep}" \
				"BC:${lib_barcode}" \
				"PI:${lib_fragsize}" \
				"DT:${seq_date}" \
				"PL:${seq_type}" \
				"PM:${seq_mach}" \
				"PU:${seq_run}.${seq_lane}.${lib_barcode}" \
	--outTmpDir "${output_dir}/Temporary" \
	--runThreadN ${NSLOTS:-1}

# link useful STAR sub-directory files back into the data directory
output_sam="${project_data_dir}/alignment.sam.${script_tag}"
command cp -- "${output_dir}/Aligned.out.sam" "${output_sam}"
command cp -- "${output_dir}/Log.final.out" "${temp_output_stats}"
command cp -- "${output_dir}/SJ.out.tab" "${temp_output_sjs}"

# link the outputs of this instance as the new default outputs
for loop_path in "${output_dir}" "${output_sam}"; do
	loop_dir="$(command dirname -- "${loop_path}")"
	loop_file="$(command basename -- "${loop_path}" "${script_tag}")"
	command ln --symbolic --force -- "${loop_file}" "${loop_dir}/${loop_file}"
done
