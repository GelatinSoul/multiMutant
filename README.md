# multiMutant

Sequential version is in the /regular/ folder.

After throwing together that script, I tried forking processes and running in parallel, but unfortunately it appears that rMutant creates temporary files, so multiple running instances of rMutant seem to fight over the file access. To resolve this, I decided to run each mutation inside a docker container, from which the result could be copied back to the host. Not finished implementing this yet!

Depending on how effective this proves (particularly for larger workloads) it might be doable to try running on AWS Batch or a comparable service too. Alternatively, if the file access issues prove to be inconsequential, traditional forking could be done for simplicity's sake on small workloads and combined with containers when dealing with larger job sizes.

### Todo
- make appropriate changes to script run inside docker container
  - change arguments, work size
  - copy resulting files to host
 - fix libstdc++5 being unresolved in rMutant on the container - change base image or add package in dockerfile?
- change container-spawning script to take an image name as an argument
