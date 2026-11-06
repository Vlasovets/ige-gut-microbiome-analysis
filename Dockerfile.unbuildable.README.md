# `Dockerfile.unbuildable` — retained for the record, do not use

The Dockerfile shipped with the submitted version **cannot build**. Its first
instruction is

    COPY install_apt.sh /opt/

and `install_apt.sh`, `install_rpacks.sh` and `install_pip.sh` are all absent from
the repository and from its git history. The README instructed readers to
`docker pull ovlasovets/rpy2:latest` instead, so the contents of the image that
actually produced the published results cannot be verified from this repository.

It is renamed rather than deleted so the record of what was submitted is intact.

## Use instead

    conda env create -f revision/env/slr-env.yml
    conda activate slr-env

`revision/env/versions.json` records the exact versions captured at run time,
including git SHAs for the two packages installed from non-release branches
(SpiecEasi `lowrank`, NetCoMi `develop`). TaxSEA requires R >= 4.5 and therefore
lives in a separate environment; see `revision/slurm/00c_create_taxsea_env.sh`.
