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

#### 3. `README.md`

Add the tool to the relevant numbered section in the **Pipeline summary**. If the tool belongs to a new category not yet represented, add a new numbered section in the appropriate position.
