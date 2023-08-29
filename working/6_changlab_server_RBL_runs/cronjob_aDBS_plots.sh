export PATH=$PATH:/data_store2/software/fsl/bin:/usr/local/cuda/bin:/data_store2/software/fsl/bin:/usr/local/freesurfer/bin:/usr/local/freesurfer/fsfast/bin:/usr/local/freesurfer/tktools:/usr/local/freesurfer/mni/bin:/usr/local/sge/bin:/usr/local/sge/bin/lx-amd64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/data_store2/server_tools/bin:/data_store2/software/node-v12.13.1-linux-x64/bin:/snap/bin:/usr/local/sge/bin/lx-amd64

export SGE_ROOT=/usr/local/sge

/data_store2/server_tools/bin/submit_job -q pia-batch.q -c 8 -m 30 -o /home/rleriche/cronjob_aDBS_log.txt -x /data_store2/MATLAB/R2022a/bin/matlab /home/rleriche/Analysis-rcs-data/working/6_changlab_server_RBL_runs/plot_aDBS_server_wrapper.m
