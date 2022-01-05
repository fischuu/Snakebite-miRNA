#run snakemake on the cluster
#
#$1 is target file 

# This conda module is just to make snakemake available
module load bioconda/3

snakemake -s miRNA-pipeline.smk \
          --configfile /scratch/project_2001310/TwinCalf/miRNA-pipeline_config.yaml \
          --forceall --rulegraph | dot -T png > ./workflow.png

snakemake -s miRNA-pipeline.smk \
          -j 200 \
          --latency-wait 60 \
          --use-singularity \
          --singularity-args "-B /scratch,/projappl,/scratch/myProject/tmp:/tmp" \
          --configfile /scratch/project_2001310/TwinCalf/miRNA-pipeline_config.yaml \
          --cluster-config miRNA-pipeline_puhti-config.yaml \
          --cluster "sbatch -t {resources.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={threads} --mem-per-cpu={resources.mem} -p {cluster.partition} -D {cluster.working-directory}" $1 
