# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(fs)

basedir <- here("bulk/BafA1.MEF")

md <- read_csv(here(basedir, "data/samples.csv"))

# functions ===================================================================

generate_kallisto_script <- function(sample_id,
                                     fqdir = "~/palmer_scratch/data/bulk-rna/csf1/BafA1.MEF/01_fastq",
                                     fastq.files,
                                     index = "~/csf1/genome/mmus/data/ensembl/v111/kallisto-index/Mus_musculus.GRCm39.cdna.all.kallisto.idx",
                                     basedir = "~/csf1/bulk/BafA1.MEF",
                                     library_layout = "SINGLE") {
  kallisto_command <- if (library_layout == "SINGLE") {
    glue(
      "
kallisto quant \\\\\n
-i {index} \\\\\n
-t 6 -b 100 -o {basedir}/data/02_kallisto/{sample_id} --single -l 200 -s 20 \\\\\n
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

module load kallisto

mkdir -p {basedir}/data/02_kallisto/{sample_id}

cd {fqdir}

{kallisto_command}

  "
  )
  
  return(script)
}

# generate scripts ============================================================

md <- md |>
  mutate(fastq_file = paste0(sample_id_old, ".fastq.gz")) |>
  mutate(kallisto = map2_chr(
    sample_id,
    fastq_file,
    ~ generate_kallisto_script(sample_id = .x, fastq.files = .y)
  ))

# write script files ==========================================================
fs::dir_create(here(basedir, "code/02b_kallisto_scripts"))
for (i in 1:nrow(md)) {
  writeLines(
    md$kallisto[[i]],
    here(basedir, glue("code/02b_kallisto_scripts/{md$sample_id[[i]]}.sh"))
  )
}

# write runall script =========================================================
runall <- glue("sbatch {md$sample_id}.sh")
runall <- str_flatten(runall, "\n")

runall <- str_c("#!/usr/bin/env bash",
                "cd ~/csf1/bulk/BafA1.MEF/code/02b_kallisto_scripts",
                runall, 
                sep = "\n\n")

writeLines(runall, here(basedir, "code/02c_kallisto_run-all-scripts.sh"))

# end =========================================================================

