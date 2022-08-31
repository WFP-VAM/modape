#!/bin/bash
docker run -v $PWD:/arc_modape_chain -v /var/storage:/var/storage \
  --rm -it arc-modape bash
