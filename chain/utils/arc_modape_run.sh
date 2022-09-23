#!/bin/bash
#
# Runs the provided script inside the container. Usage: ./arc_modape_run.sh policy_dataset_test.sh
#
docker run -v $PWD:/arc_modape_chain -v /var/storage:/var/storage \
  --rm -it arc-modape bash -c '/arc_modape_chain/$0 |& tee -a /var/storage/modape_$(date -d "today" +"%Y%m%d%H%M").log' $1
