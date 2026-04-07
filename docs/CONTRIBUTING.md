# `Clinical-Genomics/oncorefiner`: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving Clinical-Genomics/oncorefiner.

We try to manage the required tasks for Clinical-Genomics/oncorefiner using GitHub issues, you probably came to this page when creating one.
Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome!
Contributions to the code are even more welcome ;)

## Table of contents

- [General](#general)
  - [Contribution workflow](#contribution-workflow)
    - [Pull Requests](#pull-requests)
      - [PR title conventions](#pr-title-conventions)
      - [Review](#review)
  - [Installation and dependencies for development](#installation-and-dependencies-for-development)
  - [Running tests](#running-tests)
    - [Lint tests](#lint-tests)
    - [Pipeline tests](#pipeline-tests)
  - [Patch](#patch)
  - [Nextflow version bumping](#nextflow-version-bumping)
  - [Update nf-core template](#update-nf-core-template)
  - [GitHub Codespaces](#github-codespaces)
  - [Adding citations](#adding-citations)
    - [1. `CITATIONS.md`](#1-citationsmd)
    - [2. `subworkflows/local/utils_nfcore_oncorefiner_pipeline/main.nf`](#2-subworkflowslocalutils_nfcore_oncorefiner_pipelinemainnf)
    - [3. `README.md`](#3-readmemd)
  - [Images and figures](#images-and-figures)
- [Coding conventions](#coding-conventions)
  - [Architecture & structure](#architecture--structure)
  - [Adding a new step](#adding-a-new-step)
  - [Channels](#channels)
    - [Naming schemes](#naming-schemes)
  - [Parameters](#parameters)
  - [Publishing](#publishing)
  - [Configuration](#configuration)
  - [Writing tests](#writing-tests)
  - [Style](#style)

## General

### Contribution workflow

If you'd like to write some code for Clinical-Genomics/oncorefiner, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [Clinical-Genomics/oncorefiner issues](https://github.com/Clinical-Genomics/oncorefiner/issues) to avoid duplicating work. If there isn't one already, please create one so that others know you're working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [Clinical-Genomics/oncorefiner repository](https://github.com/Clinical-Genomics/oncorefiner) to your GitHub account
3. Make the necessary changes / additions within your forked repository following [Pipeline conventions](#pipeline-contribution-conventions)
4. Use `nf-core pipelines schema build` and add any new parameters to the pipeline JSON schema (requires [nf-core tools](https://github.com/nf-core/tools) >= 1.10).
5. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

#### Pull Requests

When opening a pull request to suggest changes to the code, please make sure to follow the [Pipeline contribution conventions](#pipeline-contribution-conventions) for the code and to fill in the necessary information in the pull request template as well as address all points in the `PR checklist`.

##### PR title conventions

We have implemented a standardised PR title format to make it easier to understand the type of change being proposed at a glance.
Addionally, there is an automated check for every PR that will only allow mergere if the title adheres to one of the following formats:

- feat: A new feature
- fix: A bug fix
- docs: Documentation only changes
- style: Changes that do not affect the meaning of the code (white-space, formatting, missing semi-colons, etc)
- refactor: A code change that neither fixes a bug nor adds a feature
- perf: A code change that improves performance
- test: Adding missing tests or correcting existing tests
- build: Changes that affect the build system or external dependencies (example scopes: gulp, broccoli, npm)
- ci: Changes to our CI configuration files and scripts (example scopes: Travis, Circle, BrowserStack, SauceLabs)
- chore: Other changes that don't modify src or test files
- revert: Reverts a previous commit

##### Review

When reviewing a PR, make sure to check that:

- The code follows the [Pipeline contribution conventions](#pipeline-contribution-conventions).
- The information in the PR (and related issue) is clear and sufficient to understand the change and the motivation for it - title, description and entry in `CHANGELOG.md`, if applicable.
- All the items in the `PR checklist` have been addressed, the changes are well documented and the tests are passing.

Be positive and constructive in your review, and whenever possible offer suggestions for improvement rather than just pointing out issues.

### Installation and dependencies for development

In order to run the pipeline, develop and test your changes locally, we recommend that you set up:

- A conda environment with `nextflow` and `nf-core` tools. - For this, follow the instructions from the [nf-core documentation](https://nf-co.re/docs/nf-core-tools/installation). Additional information about [Installation of nf-core dependencies](https://nf-co.re/docs/usage/getting_started/installation/) is also available, if needed.
- Install Docker (https://www.docker.com/products/docker-desktop/) and make sure the daemon is running when you want to run the tests locally.

Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Running tests

You have the option to test your changes locally by running the pipeline. For receiving warnings about process selectors and other `debug` information, it is recommended to use the debug profile. Execute all the tests with the following command:

```bash
nf-test test --profile debug,test,docker --verbose
```

You can also run `nf-test test <path>` for a single given test.

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

There are typically two types of tests that run:

#### Lint tests

`nf-core` has a [set of guidelines](https://nf-co.re/developers/guidelines) which all pipelines must adhere to.
To enforce these and ensure that all pipelines stay in sync, we have developed a helper tool which runs checks on the pipeline code. This is in the [nf-core/tools repository](https://github.com/nf-core/tools) and once installed can be run locally with the `nf-core pipelines lint <pipeline-directory>` command.

If any failures or warnings are encountered, please follow the listed URL for more documentation.

#### Pipeline tests

Each `nf-core` pipeline should be set up with a minimal set of test-data.
`GitHub Actions` then runs the pipeline on this data to ensure that it exits successfully.
If there are any test failures then the automated check has status set to fail.
These tests are run both with the latest available version of `Nextflow` and also the minimum required version that is stated in the pipeline code.

### Patch

:warning: Only in the unlikely and regretful event of a release happening with a bug.

- On your own fork, make a new branch `patch` based on `upstream/main` or `upstream/master`.
- Fix the bug, and bump version (X.Y.Z+1).
- Open a pull-request from `patch` to `main`/`master` with the changes.

### Nextflow version bumping

If you are using a new feature from core Nextflow, you may bump the minimum required version of nextflow in the pipeline with: `nf-core pipelines bump-version --nextflow . [min-nf-version]`

### Update nf-core template

Since this is not an nf-core pipeline, the nf-core template is not automatically updated in the `TEMPLATE` branch. Follow these step to update the template:

1. Update the `TEMPLATE` branch by running `nf-core pipelines sync`. Fix any merge conflicts and open a PR to then merge the changes.
1. Open a PR to merge the `TEMPLATE` branch into `dev` to update the template files in the main codebase.

### GitHub Codespaces

This repo includes a devcontainer configuration which will create a GitHub Codespaces for Nextflow development! This is an online developer environment that runs in your browser, complete with VSCode and a terminal.

To get started:

- Open the repo in [Codespaces](https://github.com/Clinical-Genomics/oncorefiner/codespaces)
- Tools installed
  - nf-core
  - Nextflow

Devcontainer specs:

- [DevContainer config](.devcontainer/devcontainer.json)

### Adding citations

When adding a new tool to the pipeline, update the following three locations:

#### 1. `CITATIONS.md`

Add an entry for the tool in alphabetical order under `## Pipeline tools`. If the tool has a publication, include a `>` citation block:

```markdown
- [ToolName](https://link-to-paper-or-repo)

  > Author A, Author B. Title. Journal. Year;vol(issue):pages. doi:...
```

#### 2. `subworkflows/local/utils_nfcore_oncorefiner_pipeline/main.nf`

Add citation text and bibliography entries inside `toolCitationText()` and `toolBibliographyText()`. Both functions are structured identically — group the tool's entry under the relevant category variable (e.g. `align_text`, `qc_bam_text`, `preprocessing_text`, `snv_annotation_text`). Mirror any conditional logic that gates the tool's execution (e.g. skip params, analysis type, or input content) so the citation only appears when the tool actually runs:

```groovy
// toolCitationText()
qc_bam_text = [
    ...,
    (condition) ? "ToolName (Author et al., Year)," : ""
]

// toolBibliographyText()
qc_bam_text = [
    ...,
    (condition) ? "<li>Author A, Author B. Title. Journal. Year. doi:...</li>" : ""
]
```

For tools that run only when the input samplesheet contains a particular file type, use a helper function rather than a param check — see `hasSpringInput()` as an example.

#### 3. `README.md`

Add the tool to the relevant numbered section in the **Pipeline summary**. If the tool belongs to a new category not yet represented, add a new numbered section in the appropriate position.

### Images and figures

For overview images and other documents we follow the nf-core [style guidelines and examples](https://nf-co.re/developers/design_guidelines).

## Coding conventions

To make the `Clinical-Genomics/oncorefiner` code and processing logic more understandable for new contributors and to ensure quality, we semi-standardise the way the code and other contributions are written.

### Architecture & structure

- **Use subworkflows** — Don't add logic to `workflows/oncorefiner.nf` that is specific to a subworkflow. Create new subworkflows as needed under `workflows/` and import them into `workflows/oncorefiner.nf`.
- **Reuse over duplication** — `ENSEMBLVEP_VEP` and `BCFTOOLS_VIEW` are intentionally included multiple times under different aliases. Follow this pattern before creating a near-identical subworkflow.
- **nf-core modules take precedence** — prefer a module from modules/nf-core/ over writing a local one. Only add to modules/local/ when no nf-core module exists or the use case is too pipeline-specific.

### Adding a new step

If you wish to contribute a new step, please use the following coding standards:

1. Define the corresponding input channel into your new process from the expected previous process channel.
2. Write the process block (see below).
3. Define the output channel if needed (see below).
4. Add any new parameters to `nextflow.config` with a default (see below).
5. Add any new parameters to `nextflow_schema.json` with help text (via the `nf-core pipelines schema build` tool).
6. Add sanity checks and validation for all relevant parameters.
7. Perform local tests to validate that the new code works as expected.
8. If applicable, add a new test in the `tests` directory.
9. Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name clean up and module plots are in the appropriate order. If applicable, add a [MultiQC](https://https://multiqc.info/) module.
10. Add a description of the output files and if relevant any appropriate images from the MultiQC report to `docs/output.md`.

### Channels

- **Conditional channels**: always initialize to `channel.empty()` before any `if` block that may or may not assign them. Never leave a channel potentially undefined.

#### Naming conversions

Please use the following naming schemes, to make it easy to understand what is going where.

- initial process channel: `ch_output_from_<process>`
- intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

### Parameters

- `params` must only be accessed in the main unnamed workflow (`workflow` in `main.nf`). Subworkflows and named workflows receive all values as explicit `val_*` arguments. Never reference `params` directly inside a subworkflow.
- Parameters should be initialised/defined with default values within the `params` scope in `nextflow.config`. Don't hardcode values that a user might reasonably want to change. Once added, run nf-core pipelines schema build to register them in nextflow_schema.json.

<!-- TODO: Add information about parameters for skipping tools when that logic is decided. -->

### Publishing

- Build a single `ch_publish` channel inside each subworkflow by mixing all publishable outputs into `[destination, value]` tuples.
- The emit name must always be `publish = ch_publish` — never the bare shorthand.
- Group channels that share a destination with `mix` first, then apply **one** `.map` per destination group — never one map per channel.
- If your subworkflow calls inner subworkflows, always mix their `.out.publish` into the outer `ch_publish`. Never discard it.
- Remove the corresponding `publishDir` entry from `conf/modules/` when adding a process to `ch_publish`.

### Configuration

- Process-level options go in `conf/modules/<subworkflow_name>.config`, not inline in the subworkflow `.nf` file.
- Only `ext.args`, `ext.args2`, and `ext.prefix` belong in module configs. Don't add business logic there.
- Conditional behavior (e.g. save as CRAM vs BAM) is handled in the subworkflow via `channel.empty()` gating — not via config-level flags.
- Process resource requirements (CPUs / memory / time) go in `conf/base.config` using `withLabel:` selectors so they can be shared across processes. Use `${task.cpus}` and `${task.memory}` in `script:` blocks to apply them dynamically.

### Writing tests

- Every subworkflow should have a test at `subworkflows/local/<name>/tests/main.nf.test`.
- Use `-stub` in the `when:` block only when real test data is difficult to generate. Prefer running with real data where it is reasonably available.
- Snapshot files (`*.nf.test.snap`) are committed alongside tests — update them when outputs change.
- Pipeline-level tests live in `tests/` and cover the default and `test` profiles.

### Style

- Both `take:` and `emit:` block entries require an inline type comment. Use `name // type: [mandatory|optional] description` for `take:` and `name = value // channel: [type description]` for `emit:`. Always include the comment — never leave an entry uncommented.

  ```groovy
  take:
      ch_vcf                // channel: [mandatory] [ val(meta), path(vcf) ]
      ch_reduced_penetrance // channel: [optional]  [ path(penetrance) ]
      val_aligner           // string:  [mandatory] aligner name (bwa/bwamem2/bwameme)
      process_with_sort     // Boolean

  emit:
      vcf     = ch_vcf      // channel: [ val(meta), path(vcf) ]
      publish = ch_publish  // channel: [ val(destination), val(value) ]
  ```

- Avoid using the `.set {ch_*}` operator to create new channels. Use `ch_* = <...>` whenever possible.
- Intermediate publish channels in `workflows/oncorefiner.nf` follow the `ch_<subworkflow_name>_publish` naming convention and are assigned immediately after the subworkflow call, not inline in the emit block.
- Initialize all `ch_*_publish` variables at the top of the `main:` block alongside `ch_multiqc_files`.
