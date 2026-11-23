# User guide

For researchers using the website. No programming needed.

## Start the website (on your own computer)

In VS Code: **Terminal → New Terminal**, then:

```powershell
.venv\Scripts\flask.exe --app microplatform.web run
```

Wait for `Running on http://127.0.0.1:5000`, open that address in your
browser, and keep the terminal open. To stop the website, press **Ctrl+C** in the
terminal.

## An analysis in four steps

1. **Upload.** On the home page you can pick your question ("Which bacteria
   differ between my groups?"). The upload box then names the files that
   question needs. Drop your files: which one is the counts table, the labels
   file or the metabolite table is detected automatically.
2. **Data check.** You see what was found: samples, taxa, taxonomy depth,
   sparsity, read depth and label type. Red problems block the analysis and
   say how to fix them. Amber warnings (for example, a sample with far more
   reads than the rest, or samples low on reads or on taxa detected) do not block it.
3. **Plan.** A task is recommended, with the reasons. You can choose
   another. If a task needs a file you have not uploaded, it says which one,
   and you can add a labels file or a metabolite table right there. You can
   also choose the **taxonomy level** (see below); species is the default.
4. **Results.** After you click Run, the analysis runs on the server and the
   page shows "Running…". You can close the page: it updates by itself, and
   the link shown there brings you back to the results at any time. When it is
   done you get a summary, figures, the processed table and a record of every
   setting used. (If the lab has switched email on, you can also leave your
   address to receive the link; it is deleted once the email is sent.)

## Taxonomy level

MIPMLP groups the bacteria at one taxonomy level before any tool runs, and every tool of the
analysis uses that table. Species is the default (the MIPMLP paper's setting). The Plan step offers
every level from phylum down to the deepest one your table names for most taxa; if your table names
species for only some taxa, "species" groups the rest at their deepest named rank.

- A coarser level merges bacteria into fewer features: findings name genera (or families…) instead
  of species.
- **Phylum and class** leave few features, so iMic and GIMIC (which turn each sample into an image of
  the taxonomy tree) are less informative there. At phylum, GIMIC cannot measure distances at all and
  says so.
- **miMic below species** runs its test and shows its results, but not its cladogram and its
  real-vs-shuffled figure: the authors' plotting code needs species level.
- The data overview (composition, diversity, sparsity) does not change with the level: it describes
  your uploaded table. To see composition at a finer level, use "Show composition at" in the overview.
- A table that is already preprocessed (used as it is) is not grouped again, so it has no level.

## Which files each task needs

| Task | Tool | Files |
|---|---|---|
| Just preprocess | MIPMLP | counts table |
| Find bacteria that differ between groups | miMic | counts + labels (two groups, or a continuous label) |
| Predict the outcome | iMic | counts + labels (two groups) |
| See how samples group | GIMIC | counts (labels optional, for colouring) |
| Link bacteria to metabolites | LOCATE | counts + metabolite table (at least 20 samples in both) |

## File formats

**Counts table** (CSV or TSV): one row per sample, one column per bacterium, with
read counts or relative abundances. The first column holds the sample IDs (named
`ID`). Each column's lineage is given in one of two ways.

Either a **last row** named `taxonomy`:

```
ID,       ASV1,          ASV2,          …
sample1,  120,           0,             …
sample2,  37,            5,             …
taxonomy, k__Bacteria;…, k__Bacteria;…, …
```

or the **column names** themselves (a very common export; empty ranks such as `p__`
are fine, and the same lineage may appear in two columns):

```
ID,         k__Bacteria;p__Firmicutes;…, k__Bacteria;p__;c__;…, …
ERR1549010, 0.0312,                      0,                     …
ERR1549011, 0,                           0.0047,                …
```

Both give exactly the same analysis; the data check says which one was found.

**YaMAS output (16S)**: upload YaMAS's two files together, as YaMAS wrote them.
`otu.csv` has one row per OTU (its ID first) and one column per sample, with read
counts; `taxonomy.csv` has the same OTU IDs and a `Taxon` column (other columns, such
as `Confidence`, are ignored).

```
OTU,     sample1, sample2, …
a1f3c9…, 120,     0,       …
```

```
OTU,     Taxon,                          Confidence
a1f3c9…, d__Bacteria; p__Firmicutes; …,  0.98
```

They are recognised by their shape (any file names, any order) and joined into one
counts table. Each lineage is written in the standard form with the lab's
`normalize_taxonomy_name` (ranks as `k__…;p__…;…`, `D_0__` and `d__` prefixes
converted, classifier confidences such as `(0.71)` removed, empty taxonomy as
"Unassigned"). The **read counts are kept** (not turned into percentages), so the data
check and the quality figures see the sequencing depth, and MIPMLP does the
normalisation. OTUs that end up with the same lineage stay separate columns: MIPMLP
groups them (its sub-PCA uses the differences among them), and richness in the data
overview counts OTUs. The pair counts as one file: with a labels file and a metabolite
table, that is four files. Every sample is kept as uploaded; if some sample IDs end in
`_1` or `_2`, the data check notes how many (in case they mark copies of the same sample,
which you would remove before uploading).

**MetaPhlAn 4 output (shotgun)**: upload the table as MetaPhlAn wrote it (a merged table with one
column per sample, the same table turned so samples are rows, or one sample's profile). It is
recognised by MetaPhlAn's `#mpa_…` line or its `clade_name` header. MetaPhlAn lists every rank
(kingdom, phylum, … species, SGB); only the **species** are kept, since the other ranks are sums of
the species (or, for SGBs, splits of them) and would count the same reads again. Lineages are written
in the standard form, and an `UNCLASSIFIED` row is left out (the data check gives its share).

- **Read counts (preferred, not open yet):** MetaPhlAn's estimated reads (`metaphlan … -t
  rel_ab_w_read_stats`). They will be kept as counts, so the depth checks work and MIPMLP normalises
  them as it does 16S counts. This is switched off until it has been confirmed on a real counts file;
  until then such a file gets a message saying so.
- **Relative abundance (%)**, MetaPhlAn's default: used as it is, and the data check says what
  that limits. There is no sequencing depth, and MIPMLP's log (log10 of the value + 0.1) brings
  species below about 0.1% close to the value for "absent".

**Labels file** (optional): sample IDs in the first column, and the outcome in
a column named `Tag`. Binary (sick/healthy), categorical (a few groups) and
continuous (a number such as age) are all recognised.

```
ID,      Tag
sample1, 1
sample2, 0
```

**Metabolite table** (optional, for linking bacteria to metabolites): sample
IDs in the first column, the same as in the counts table, then one column per
metabolite with the measured concentrations or intensities. Do not
log-transform them: 0 means "not detected", and the platform normalises the
values the way the LOCATE paper does. A table is recognised as metabolites
when it has at least 3 numeric columns and no `Tag` column.

```
ID,      Acetate, Butyrate, …
sample1, 5210.4,  830.2,    …
sample2, 4480.0,  0,        …
```

Templates of all three are on the website (in every error message, and under
"What should my files look like?"). Excel files must be saved as CSV first
(File → Save As → CSV UTF-8).

## What you get

- **The processed table** (CSV): MIPMLP's output, one row per sample.
- **Data quality** (read counts only): sequencing depth per sample, and **reads vs taxa
  detected**, one dot per sample. Samples below 25% of the median on either axis (or under 1,000
  reads) fall in a shaded band and are coloured and named; low on both usually means a failed
  sample. Few taxa with normal reads can be real, so check before removing a sample.
- **Figures**: each as PNG, SVG (with editable text) and an interactive HTML
  file that works offline, each with a CSV of the exact data plotted.
- **"Download everything"**: one ZIP with a README, the reproducibility
  record (every setting and version), the results, the figures and your input
  files.
- **A data overview**, the same for every task, at the end of the results. It
  holds the standard figures that describe a microbiome dataset. They are
  **mostly descriptive**. When your samples form groups (a labels file with
  two or more groups), the two diversity figures add the lab's group tests.
  These are **exploratory**: computed on your uploaded data, not corrected
  across the whole report, and not the tools' findings.
  - **Composition of each sample** (in the lab's figure style): on the left a
    stacked bar per sample, by phylum, sorted by Firmicutes; on the right one
    bar with the mean of all samples. The 12 phyla with the highest mean get a
    colour; the rest are "Other", in grey. Each bar is the sample's whole
    total (100%): abundance your table does not name to a phylum is its own
    light-grey segment, **Not named at phylum**, as at the finer levels (only
    if your table has any). Your file's name is the heading of the image (PNG
    at 300 dpi, and SVG).
  - **Average composition:** a pie, one per group. The colours match the bars.
  - **Composition at a finer level:** above the composition figures, "Show composition at" switches
    the bars and the pie between phylum (the default) and each finer level your table names (class,
    order, family, genus, species). Every level is drawn with the run, so switching re-runs nothing,
    and it does not change the analysis level (see "Taxonomy level"). Below phylum, the 12 most
    abundant taxa get a colour and the rest are in two greys, so each bar is still 100% of the sample:
    **Other (named, rare)** (grey) sums the other named taxa, and **Not named at <level>** (light grey)
    is the abundance your table does not name down to that level (including "unclassified" or
    "uncultured"). At genus and species the greys can be large: many rare taxa, and taxa named only to
    a coarser level. The caption says how much the 12 cover and how large each grey is. Names are kept
    as your table writes them (GTDB suffixes included); where two taxa share a name, the family is added.
  - **Look at one taxon:** under the composition figures, pick a level and start typing a name: the list
    holds every taxon your table names at that level, most abundant first ("uncultured", "unclassified" and
    similar are not offered: they are not one taxon). "Show" draws one bar per sample, the taxon's share of
    that sample's whole total, highest to lowest (hover for the reads, for a table of counts); with groups,
    also a box per group with every sample as a dot. The line under it says in how many samples the taxon
    is found, and its mean and median share. It is descriptive only, with no test (to find bacteria that
    differ between groups, run miMic). Taxa are counted and named as in the composition figures. The data
    (CSV) has each sample's share and its reads (or, for relative abundances, the table's own value).
  - **Alpha diversity:** Shannon index and richness, as a box per group with
    every sample shown. With groups: the test is chosen by normality
    (Shapiro-Wilk in each group: t-test or ANOVA if every group looks normal,
    otherwise Mann-Whitney U or Kruskal-Wallis; with three or more groups an
    overall test plus every pair). Brackets show Benjamini-Hochberg (FDR)
    corrected p-values (*** < 0.001, ** < 0.01, * < 0.05, ns). "Tests (CSV)"
    has every p-value, corrected p-value and Cohen's d. Richness is counted
    without rarefaction, so it also reflects sequencing depth.
  - **Beta diversity:** a PCoA map of Bray-Curtis distances, coloured by
    group. With groups: PERMANOVA (999 permutations) with its pseudo-F and p,
    an effect size η² (and PERMANOVA's R²), dashed ellipses where about 95% of
    each group's points would fall, and PERMDISP, which says whether the
    groups also differ in spread (if they do, part of a PERMANOVA result can
    come from spread rather than position). "Tests (CSV)" has every number.
  - **Sparsity:** taxa detected per sample (rows) and in how many samples each
    taxon appears (columns).

  All of these are computed from your own table (counts or relative
  abundances), not from the log-transformed table the tools use. A phylum
  always has the same colour, in every figure and every dataset (the lab's
  palette by default; in the package, `palette="colorblind"` gives a
  colour-blind-safer set). A table that is already log-transformed (negative
  values) gets no composition figures. In these figures only, one phylum gets one bar: GTDB
  variants (`Bacillota_A_368345`, `Bacillota_C`) are shown as their base phylum (`Bacillota`), and
  an old and a new name of one phylum (Firmicutes and Bacillota) as the one your table uses more;
  the caption says when this happened. The analysis keeps the original names.

## How good are LOCATE's predictions?

LOCATE's results report **one training run**, which is what you can expect in
practice: the median correlation between measured and predicted metabolite
levels on held-out samples. The LOCATE paper's per-dataset figures show, for
each metabolite, the best of 16 runs. On two of the paper's public datasets
that made the median about 0.13 higher (FRANZOSA: 0.37 for one run vs 0.50
best-of-16; YACHIDA: 0.14 vs 0.28), so expect lower numbers here than in the
paper's figures. Details: `validation/locate/`.

## Try it without data

Use **Example** in the top menu: a 40-sample slice of the miMic authors'
example data. Its data check flags one real oddity: sample a498 has about
10 million reads, 90% of them from one genus.

## Citing

The results ZIP lists the papers for the methods used. They are also on the
**Methods** and **About** pages.
