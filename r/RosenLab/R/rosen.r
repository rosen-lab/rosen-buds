# Some useful base definitions for various Rosen Lab functions

ROOT <- normalizePath(Sys.getenv('ROSENLAB_BASE'))

RAW_DATA_REPOSITORY <- file.path(RosenLab:::ROOT,'Data')
GENOME_DATA_REPOSITORY <- file.path(RosenLab:::RAW_DATA_REPOSITORY,'Genomes')
MOTIF_DATA_REPOSITORY <- file.path(RosenLab:::RAW_DATA_REPOSITORY,'Motifs')
SEQUENCING_DATA_REPOSITORY <- file.path(RosenLab:::RAW_DATA_REPOSITORY,'Sequencing')

PROJECT_REPOSITORY <- file.path(RosenLab:::ROOT,'Project')

SOFTWARE_REPOSITORY <- file.path(RosenLab:::ROOT,'Software')
