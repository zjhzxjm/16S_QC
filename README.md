# 16S_QC

Copy snakefile and cluster.json to your raw data directory. Modify the LIBRARY, BARCODE_NAME to your project.
The raw data filename must be "{LIBRARY}.R1.fastq.gz and {LIBRARY}.R2.fastq.gz"

LIBRARY = "H16A28P250-1"

run the shell:

export PATH=/home/xujm/soft/Anaconda/anaconda3/bin:$PATH

source activate snakemake

snakemake -j 80 --cluster-config cluster.json -c "qsub -cwd -l vf={cluster.vf} -e {cluster.err} -o {cluster.std}" --jn {params.job_name}{jobid} -T

