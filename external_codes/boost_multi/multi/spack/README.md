# Spack recipe for Boost.Multi

This directory is a self-contained [Spack](https://spack.io) package repository for
**boost-multi** (a header-only library). It lets anyone install the library with Spack,
and serves as the source for upstreaming the recipe to Spack's central package collection.

```
spack/
├── repo.yaml                          # declares this as a Spack repo (namespace: correaa)
└── packages/
    └── boost-multi/
        └── package.py                 # the recipe (header-only Package: just installs headers)
```

## 1. Install Spack (not in apt — it's a git checkout)

```bash
git clone -c feature.manyFiles=true --depth=2 https://github.com/spack/spack.git ~/spack
. ~/spack/share/spack/setup-env.sh      # add this line to ~/.bashrc to make it permanent
spack compiler find                     # register the system compiler (avoids building gcc from source)
```

## 2. Use this recipe locally (for you, or anyone who clones this repo)

```bash
spack repo add /path/to/boost-multi/spack
spack info boost-multi                   # sanity-check the recipe parses
spack spec boost-multi@develop           # show the (tiny) dependency tree
spack install boost-multi@develop        # build branch 'develop' (no checksum needed)
spack load boost-multi                   # put the headers on the compiler include path
```

A tagged release needs a source hash; generate it automatically and paste it over the
placeholder `sha256` in `package.py`:

```bash
spack checksum boost-multi 0.90.0
```

## 3. Make it available to *everyone* (upstream to Spack)

Spack has **no separate registry** — its package collection *is* the
[`spack/spack`](https://github.com/spack/spack) GitHub monorepo. Publishing = a pull request,
analogous to the Conan `conan-center-index` flow.

1. **Fork** https://github.com/spack/spack and clone your fork.
2. **Copy the recipe** into the built-in repo (drop the local `repo.yaml`; `builtin` already
   has its own namespace):
   ```bash
   mkdir -p var/spack/repos/builtin/packages/boost-multi
   cp /path/to/boost-multi/spack/packages/boost-multi/package.py \
      var/spack/repos/builtin/packages/boost-multi/package.py
   ```
3. **Make sure tagged versions have real hashes** (`spack checksum boost-multi`), not the
   placeholder.
4. **Pass Spack's linters/checks** locally:
   ```bash
   spack style                          # formatting / flake8
   spack audit packages boost-multi     # metadata sanity (homepage, license, ...)
   spack install boost-multi            # it must actually build
   ```
5. **Open a pull request** against `spack/spack`. A maintainer reviews and merges. After that,
   anyone can simply run `spack install boost-multi` with no extra setup.

Until the PR is merged, others can use the library immediately by pointing Spack at this repo:
`spack repo add <path>/boost-multi/spack`.

## Notes

- The recipe is a plain `Package` (header-only): it installs only the headers, like the Conan
  `header-library` recipe, so it pulls **no** build dependencies (no cmake/ncurses/etc.).
  Dependents that `depends_on("boost-multi")` get the include path automatically.
- If you instead want installed CMake config files for `find_package(boost-multi)`, switch the
  base class to `CMakePackage` — at the cost of a `cmake` build dependency.
- Keep Spack (and its package recipes) current with `git -C ~/spack pull`.
