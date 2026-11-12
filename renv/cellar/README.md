# Pinned pcmt source

`pcmt_b92c1ccb4f8982fb63d35d7b0986cb3fc5af13bd.tar.gz` is an R source package built from the exact git revision
`b92c1ccb4f8982fb63d35d7b0986cb3fc5af13bd` of
<https://github.com/ZiangNiu6/pcmt>. It includes the original GPL (>= 3)
license declaration, source, documentation, and tests.

The lockfile retains the original GitHub revision. During restoration, renv
prefers this bundled source archive over downloading that revision. Reviewers
therefore do not need GitHub credentials for `pcmt`. All other package pins
are unchanged. See [renv's package cellar documentation](https://rstudio.github.io/renv/articles/cellar.html).

Run `./reproduction-code/setup_environment.sh` from the repository root, or
`./setup_environment.sh` from a standalone copy of `reproduction-code/`.
The setup script verifies the SHA-256 below, selects this cellar explicitly,
and bootstraps the lockfile's exact renv version into `renv/bootstrap/` if it
is not already available. Other dependencies still require their public
CRAN, Bioconductor, or GitHub sources during a clean restore.

SHA-256: `82b6dc0686a376848a44b9baacca9ac3229509902fd5ff277417ba93bbce8d0e`

Build command: `R CMD build --no-build-vignettes --no-manual pcmt` on an
archive of that commit. No source changes were made.
