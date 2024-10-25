# Arun Chavan
# Started: 2024-10-24

# background ==================================================================

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(fs)

gse <- "GSE66362"

basedir <- here("drer-aging", gse)
fqdir <- here("~/palmer_scratch/geo-data", gse)

srp_meta <- read_tsv(here(basedir, "data/sample-metadata.tsv")) %>% 
  select(run, title, library_layout)

# functions ===================================================================

generate_kallisto_script <- function(
    sample_id, 
    fastq.files, 
    library_layout,
    index = "~/csf1/genome/drer/GRCz11/results/Danio_rerio.GRCz11.cdna.all.kallisto.idx"
) {
  
  kallisto_command <- if (library_layout == "SINGLE") {
    glue(
      "
kallisto quant \\\\\n
-i {index} \\\\\n
-t 6 -b 100 -o {basedir}/data/02_kallisto/{sample_id} --single -l 200 -s 20 \\\\\n
{fastq.files}
      "
    )
  } else if (library_layout == "PAIRED") {
    glue(
      "
kallisto quant \\\\\n
-i {index} \\\\\n
-t 6 -b 100 -o {basedir}/data/02_kallisto/{sample_id} \\\\\n
{fastq.files}
      "
    )
  }
  
  script <- glue(
    "
#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_{sample_id}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load Python
module load kallisto

mkdir -p {basedir}/data/02_kallisto/{sample_id}

cd {fqdir}

{kallisto_command}

  "
  )
  
  return(script)
}

# generate fastq file lists ===================================================
srp_meta$fastq_files <- NA_character_
for (i in seq_len(nrow(srp_meta))) {
  srp_meta$fastq_files[[i]] <- if (srp_meta$library_layout[[i]] == "SINGLE") {
    here(fqdir, srp_meta$run[[i]], paste0(srp_meta$run[[i]], ".fastq"))
  } else if (srp_meta$library_layout[[i]] == "PAIRED") {
    paste0(
      here(fqdir, srp_meta$run[[i]], paste0(srp_meta$run[[i]], "_1.fastq")),
      " ",
      here(fqdir, srp_meta$run[[i]], paste0(srp_meta$run[[i]], "_2.fastq"))
    )
  } else {
    NA_character_
  }
}

# generate kallisto scripts ===================================================

srp_meta$kallisto_script <- NA_character_

for (i in seq_len(nrow(srp_meta))) {
  srp_meta$kallisto_script[[i]] <- generate_kallisto_script(
    sample_id = srp_meta$run[[i]], 
    fastq.files = srp_meta$fastq_files[[i]], 
    library_layout = srp_meta$library_layout[[i]]
  )
}

# writeLines(srp_meta$kallisto_script[[4]])

# write script files ==========================================================
fs::dir_create(here(basedir, "code/02b_kallisto_scripts"))
for (i in 1:nrow(srp_meta)) {
  writeLines(
    as.character(srp_meta$kallisto_script[i]), 
    here(basedir, glue("code/02b_kallisto_scripts/{srp_meta$run[i]}.sh"))
  )
}

# write run_all script ========================================================
runall <- glue("sbatch {srp_meta$run}.sh")

runall <- str_flatten(runall, collapse = "\n")

runall <- str_c("#!/usr/bin/env bash",
                glue("cd /home/arc78/csf1/drer-aging/{gse}/code/02b_kallisto_scripts"),
                runall, 
                sep = "\n\n")

writeLines(runall, here(basedir, "code/02c_kallisto_run-all-scripts.sh"))


# end =========================================================================
