pipelineFolder="/users/fischerd/git/Snakebite-miRNA"
projectFolder="/scratch/project_2001310/TestProject"

# This conda module is just to make snakemake available
module load snakemake

# Here are the tmp folder for Apptainer/Singularity/Docker defined, depends on your system
export APPTAINER_TMPDIR="/scratch/project_2001310/tmp"
export APPTAINER_CACHEDIR="/scratch/project_2001310/tmp"

snakemake -s $pipelineFolder/Snakebite-miRNA.smk \
          --configfile $projectFolder/Snakebite-miRNA_config.yaml \
          --cluster-config $projectFolder/Snakebite-miRNA_server-config.yaml \
          --forceall --rulegraph | dot -T png > $projectFolder/workflow.png

snakemake -s $pipelineFolder/Snakebite-miRNA.smk \
          -j 200 \
          --latency-wait 60 \
          --use-singularity \
          --singularity-args "-B /scratch,/projappl,/users,/dev/shm:/tmp" \
          --configfile $projectFolder/Snakebite-miRNA_config.yaml \
          --cluster-config $pipelineFolder/Snakebite-miRNA_server-config.yaml \
          --cluster "sbatch -t {resources.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={threads} --mem-per-cpu={resources.mem} -p {cluster.partition} -D {cluster.working-directory}" \
          --scheduler greedy \
          $@ 
