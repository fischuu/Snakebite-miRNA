pipelineFolder="/users/fischerd/git/Pipeline-miRNA"
projectFolder="/scratch/project_2001310/TestProject"

# This conda module is just to make snakemake available
module load bioconda/3
source activate Snakemake

export SINGULARITY_TMPDIR="/scratch/project_2001310/tmp"
export SINGULARITY_CACHEDIR="/scratch/project_2001310/tmp"

snakemake -s $pipelineFolder/miRNA-pipeline.smk \
          --configfile $projectFolder/miRNA-pipeline_config.yaml \
          --cluster-config $pipelineFolder/miRNA-pipeline_puhti-config.yaml \
          --forceall --rulegraph | dot -T png > $projectFolder/workflow.png

snakemake -s $pipelineFolder/miRNA-pipeline.smk \
          -j 200 \
          --latency-wait 60 \
          --use-singularity \
          --singularity-args "-B /scratch,/projappl,/users,/dev/shm:/tmp" \
          --configfile $projectFolder/miRNA-pipeline_config.yaml \
          --cluster-config $pipelineFolder/miRNA-pipeline_puhti-config.yaml \
          --cluster "sbatch -t {resources.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={threads} --mem-per-cpu={resources.mem} -p {cluster.partition} -D {cluster.working-directory}" \
          --scheduler greedy \
          $@ 