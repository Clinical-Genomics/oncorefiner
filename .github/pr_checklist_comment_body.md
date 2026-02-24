## PR checklist

    - [ ] Fill in description of the PR and link to any relevant issues.
    - [ ] If you've fixed a bug or added code that should be tested, add tests!
    - [ ] If you've added a new tool, follow the pipeline conventions in the [contribution docs](https://github.com/nf-core/oncorefiner/tree/master/.github/CONTRIBUTING.md).
    - [ ] Make sure your code lints (`nf-core pipelines lint`).
    - [ ] Ensure the test suite passes (`nextflow run . -profile test,docker --outdir <OUTDIR>`).
    - [ ] Check for unexpected warnings in debug mode (`nextflow run . -profile debug,test,docker --outdir <OUTDIR>`).
    - [ ] Usage Documentation in `docs/usage.md` is updated.
    - [ ] Output Documentation in `docs/output.md` is updated.
    - [ ] `README.md` is updated (including new tool citations and authors/contributors).
