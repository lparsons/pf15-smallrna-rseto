# Changelog

All notable changes to this workflow will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## 1.0.0 (2024-05-21)


### Features

* add genelist summary ([b833b56](https://github.com/lparsons/pf15-smallrna-rseto/commit/b833b56adf368c0d6cf71508d454fda79a07e486))
* add MultiQC report ([9c11602](https://github.com/lparsons/pf15-smallrna-rseto/commit/9c11602af237904f9a67e878cafbc3540e8f3034))
* align reads to pf15 genome ([5591194](https://github.com/lparsons/pf15-smallrna-rseto/commit/5591194e26f508b9a0fb037bfa8ae0a69b476f36))
* annotate homologous regions with gene ids ([6033f5b](https://github.com/lparsons/pf15-smallrna-rseto/commit/6033f5b559c2ee6e39e7095bb8d8192e66bff265))
* convert blast homology to bed format ([aa3cf7f](https://github.com/lparsons/pf15-smallrna-rseto/commit/aa3cf7f5a5a472f6fba9a9ed48ccb927137235e4))
* filter homologous regions by gene expression ([011740a](https://github.com/lparsons/pf15-smallrna-rseto/commit/011740a522cc36ceafdd66436165180dc9208e6a))
* filter putative smallrna on sample coverage ([ddbf5a7](https://github.com/lparsons/pf15-smallrna-rseto/commit/ddbf5a79f297b1a44a23c0b5827e8188278da4fe))
* find homologous regions with BLAST ([9b7a0c6](https://github.com/lparsons/pf15-smallrna-rseto/commit/9b7a0c6e38b4b59c6d083c96edfa863ad6b86015))
* generate gene id to transcript id mapping ([c9e7bbf](https://github.com/lparsons/pf15-smallrna-rseto/commit/c9e7bbfada2d860e9087bdea2182139918c243e6))
* generate pf15 genome coverage ([83d671a](https://github.com/lparsons/pf15-smallrna-rseto/commit/83d671a95f38dc4a4b8731af57a1c9a4f56c92f6))
* generate report ([6eabbfb](https://github.com/lparsons/pf15-smallrna-rseto/commit/6eabbfbc8ee5431f86a3ff6be537f7b8085855c3))
* generate unqiue gene lists ([1acf2e5](https://github.com/lparsons/pf15-smallrna-rseto/commit/1acf2e5bbe7c01e8ee9115815a5f3d8092fb3acd))
* homology overlap downregulated filter ([85f490b](https://github.com/lparsons/pf15-smallrna-rseto/commit/85f490b0b83b1b9576704a670960de8e8dc246e3))
* homology overlap with smallrna ([cc4b819](https://github.com/lparsons/pf15-smallrna-rseto/commit/cc4b819f5e41b3663850bbf49319930bfb709dd2))
* organize reference config by organism ([2ea2ed8](https://github.com/lparsons/pf15-smallrna-rseto/commit/2ea2ed8b6a84cd4f1b301b7419eb84261e89f37a))
* trim adapters from fastq ([4bd8c5a](https://github.com/lparsons/pf15-smallrna-rseto/commit/4bd8c5ab618348c090817f87d1cefa5601d9085d))


### Bug Fixes

* ensure homologous region covers the entire putative small RNA ([84c8be5](https://github.com/lparsons/pf15-smallrna-rseto/commit/84c8be5aee040127a5a831c0d8ab177aa86a1472))
* ensure homologous region is entirely within putative small RNA ([84c8be5](https://github.com/lparsons/pf15-smallrna-rseto/commit/84c8be5aee040127a5a831c0d8ab177aa86a1472))
* update release-please config ([a3dfbe3](https://github.com/lparsons/pf15-smallrna-rseto/commit/a3dfbe3f2d1a0871a7d68f0fae4c26521a8ecdf4))
* update test files to ensure final results ([5e97050](https://github.com/lparsons/pf15-smallrna-rseto/commit/5e970502354fd8a5a9145117099f988270cf126b))


## [2024-05-08]

### Added

- update report to include small rnas and column definitions

## [2024-05-07]

### Added

- initial report delivered


[unreleased]: https://github.com/lparsons/pf15-smallrna-rseto/compare/d70c35adea3753c75293a081e5695a3e3dc50d1b...HEAD
[2024-05-08]: https://github.com/lparsons/pf15-smallrna-rseto/compare/1d53dd1976e3ad680064d513650daa19078ebcfe...d70c35adea3753c75293a081e5695a3e3dc50d1b
[2024-05-07]: https://github.com/lparsons/pf15-smallrna-rseto/commit/1d53dd1976e3ad680064d513650daa19078ebcfe
