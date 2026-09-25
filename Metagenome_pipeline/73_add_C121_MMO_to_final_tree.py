from pathlib import Path
import re


path = Path(
    "scripts/92_plot_final_tree_MCA0421_quality_refined.R"
)


if not path.exists():

    raise SystemExit(
        f"ERROR: Stage 92 not found: {path}"
    )


text = path.read_text()


## ================================================================== ##
## Safety
## ================================================================== ##

if (
    "FINAL_C121_MMO_TRACK_FILE" in text
    or
    "C121_XMIN" in text
    or
    "MMO_XMIN" in text
):

    raise SystemExit(
        "ERROR: Stage 92 already appears to contain the new "
        "Cluster_00121/MMO tracks. No changes made."
    )


backup = Path(
    "scripts/92_plot_final_tree_MCA0421_quality_refined.R.pre_final_C121_MMO"
)


if not backup.exists():

    backup.write_text(
        text
    )


## ================================================================== ##
## 1. Add final track-file path
## ================================================================== ##

marker = '''TRACK_OUT <- file.path(
    OUTDIR,
    "final_650_tip_functional_tracks_refined.tsv"
)
'''


addition = marker + '''

FINAL_C121_MMO_TRACK_FILE <- file.path(
    OUTDIR,
    "extra_tracks",
    "final_C121_MMO_tree_tracks.tsv"
)
'''


if marker not in text:

    raise SystemExit(
        "ERROR: could not find TRACK_OUT path block."
    )


text = text.replace(
    marker,
    addition,
    1
)


## ================================================================== ##
## 2. Add NEW colors
##
## Existing colors are not altered.
##
## None of these exact colors occurs in the current palettes.
##
## Cluster_00121:
##     teal
##
## MMO:
##     particulate  = lime green
##     sMMO         = lavender
##     both         = near-black
##     not detected = light grey
##
## Completeness remains the only blue gradient.
## ================================================================== ##

mca_pattern = re.compile(
    r'(MCA_COLOR\s*<-\s*c\(\s*'
    r'"Present"\s*=\s*"#1B7837"[^\)]*'
    r'\)\s*)',
    re.S
)


match = mca_pattern.search(
    text
)


if not match:

    raise SystemExit(
        "ERROR: could not find MCA_COLOR definition."
    )


new_colors = match.group(1) + '''

## Cluster_00121

CLUSTER00121_COLOR <- c(

    "Present" =
        "#00A6A6"       # teal
)


## Methane monooxygenase repertoire

MMO_COLORS <- c(

    "Particulate MMO" =
        "#7CAE00",      # lime green

    "sMMO" =
        "#C77CFF",      # lavender

    "Both" =
        "#111111",      # near-black

    "Not detected" =
        "#D9D9D9"       # light grey
)

'''


text = (
    text[:match.start()]
    +
    new_colors
    +
    text[match.end():]
)


## ================================================================== ##
## 3. Join the final two-track table to ann
## ================================================================== ##

write_marker = '''write_tsv(
    ann,
    TRACK_OUT
)
'''


if write_marker not in text:

    raise SystemExit(
        "ERROR: could not find write_tsv(ann, TRACK_OUT)."
    )


join_block = '''## ================================================================== ##
## Final Cluster_00121 + methane-monooxygenase tracks
## ================================================================== ##

if (
    !file.exists(
        FINAL_C121_MMO_TRACK_FILE
    )
) {

    stop(
        "Missing final Cluster_00121/MMO track file: ",
        FINAL_C121_MMO_TRACK_FILE,
        "\\nRun Stage 93H first."
    )
}


final_extra_tracks <- read_tsv(
    FINAL_C121_MMO_TRACK_FILE,
    show_col_types = FALSE
)


if (
    nrow(
        final_extra_tracks
    ) != 650
) {

    stop(
        "Expected 650 rows in final extra-track table; found ",
        nrow(
            final_extra_tracks
        )
    )
}


if (
    n_distinct(
        final_extra_tracks$genome
    ) != 650
) {

    stop(
        "Final extra-track table does not contain 650 unique genomes."
    )
}


ann <- ann %>%
    left_join(

        final_extra_tracks %>%
            select(
                genome,
                Cluster00121_plot,
                MMO_status
            ),

        by = "genome"
    ) %>%
    mutate(

        Cluster00121_plot =
            if_else(
                Cluster00121_plot ==
                    "Present",
                "Present",
                NA_character_
            ),

        MMO_status =
            factor(
                MMO_status,
                levels =
                    c(
                        "Particulate MMO",
                        "sMMO",
                        "Both",
                        "Not detected"
                    )
            )
    )


if (
    any(
        is.na(
            ann$MMO_status
        )
    )
) {

    stop(
        "One or more final tree tips lacks an MMO status."
    )
}


''' + write_marker


text = text.replace(
    write_marker,
    join_block,
    1
)


## ================================================================== ##
## 4. Add columns to tip_geom
## ================================================================== ##

tip_pattern = re.compile(
    r'(MCA0421_plot\s*,)'
)


matches = list(
    tip_pattern.finditer(
        text
    )
)


if not matches:

    raise SystemExit(
        "ERROR: could not locate MCA0421_plot in tip_geom."
    )


## Use the LAST MCA0421_plot occurrence before ring geometry.
##
## The relevant select() occurs after construction of p_base.

ring_geom_pos = text.find(
    "## Ring geometry"
)


eligible = [
    m
    for m in matches
    if m.start() <
        ring_geom_pos
]


if not eligible:

    raise SystemExit(
        "ERROR: could not identify tip_geom MCA0421_plot field."
    )


m = eligible[-1]


text = (
    text[:m.end()]
    +
    '''

                Cluster00121_plot,

                MMO_status,'''.replace(
                    "\n                ",
                    "\n                "
                )
    +
    text[m.end():]
)


## ================================================================== ##
## 5. Ring geometry
##
## Existing:
##
##   Cyc2
##   Module 20
##   MCA0421
##   completeness
##   contamination
##
## New:
##
##   Cyc2
##   Module 20
##   MCA0421
##   Cluster_00121
##   MMO repertoire
##   completeness
##   contamination
##
## Cyc2 remains the first track, therefore genus shading is untouched.
## ================================================================== ##

geom_pattern = re.compile(
    r'COMP_XMIN\s*<-\s*'
    r'MCA_XMAX\s*\+\s*RING_GAP\s*'
    r'COMP_XMAX\s*<-\s*'
    r'COMP_XMIN\s*\+\s*RING_WIDTH',
    re.S
)


replacement_geom = '''C121_XMIN <-
    MCA_XMAX +
    RING_GAP


C121_XMAX <-
    C121_XMIN +
    RING_WIDTH


MMO_XMIN <-
    C121_XMAX +
    RING_GAP


MMO_XMAX <-
    MMO_XMIN +
    RING_WIDTH


COMP_XMIN <-
    MMO_XMAX +
    RING_GAP


COMP_XMAX <-
    COMP_XMIN +
    RING_WIDTH'''


text, n = geom_pattern.subn(
    replacement_geom,
    text,
    count=1
)


if n != 1:

    raise SystemExit(
        "ERROR: could not replace completeness ring geometry."
    )


## ================================================================== ##
## 6. Add ring outline objects
## ================================================================== ##

comp_outline_marker = '''comp_outline <- make_outline(
    outline_cells,
    COMP_XMIN,
    COMP_XMAX
)
'''


if comp_outline_marker not in text:

    ## Compact version used by one previous Stage-92 revision.
    comp_outline_marker = '''comp_outline <- make_outline(outline_cells, COMP_XMIN, COMP_XMAX)
'''


if comp_outline_marker not in text:

    raise SystemExit(
        "ERROR: could not find comp_outline definition."
    )


outline_addition = '''c121_outline <- make_outline(
    outline_cells,
    C121_XMIN,
    C121_XMAX
)


mmo_outline <- make_outline(
    outline_cells,
    MMO_XMIN,
    MMO_XMAX
)


''' + comp_outline_marker


text = text.replace(
    comp_outline_marker,
    outline_addition,
    1
)


## ================================================================== ##
## 7. Add ring data before completeness
## ================================================================== ##

comp_ring_index = text.find(
    "comp_ring <- tip_geom"
)


if comp_ring_index == -1:

    raise SystemExit(
        "ERROR: could not find comp_ring definition."
    )


ring_data = '''## Cluster_00121 - positive cells only

c121_ring <- tip_geom %>%
    filter(
        Cluster00121_plot ==
            "Present"
    ) %>%
    select(
        genome,
        Cluster00121_plot
    ) %>%
    left_join(
        outline_cells,
        by = "genome"
    ) %>%
    mutate(
        xmin =
            C121_XMIN,
        xmax =
            C121_XMAX
    )


## Methane monooxygenase repertoire - all 650 genomes

mmo_ring <- tip_geom %>%
    select(
        genome,
        MMO_status
    ) %>%
    left_join(
        outline_cells,
        by = "genome"
    ) %>%
    mutate(
        xmin =
            MMO_XMIN,
        xmax =
            MMO_XMAX
    )


'''


text = (
    text[:comp_ring_index]
    +
    ring_data
    +
    text[comp_ring_index:]
)


## ================================================================== ##
## 8. Add ring-outline plot layers immediately before comp_outline
## ================================================================== ##

plot_start = text.find(
    "## Plot"
)


if plot_start == -1:

    raise SystemExit(
        "ERROR: plot section not found."
    )


comp_layer_match = re.search(
    r'\s*geom_rect\(\s*'
    r'data\s*=\s*comp_outline\s*,',
    text[
        plot_start:
    ],
    re.S
)


if not comp_layer_match:

    raise SystemExit(
        "ERROR: could not find comp_outline plot layer."
    )


insert_at = (
    plot_start
    +
    comp_layer_match.start()
)


outline_layers = '''

    geom_rect(
        data =
            c121_outline,
        aes(
            xmin =
                xmin,
            xmax =
                xmax,
            ymin =
                ymin,
            ymax =
                ymax
        ),
        inherit.aes =
            FALSE,
        fill =
            NA,
        color =
            OUTLINE_COLOR,
        linewidth =
            OUTLINE_LINEWIDTH
    ) +


    geom_rect(
        data =
            mmo_outline,
        aes(
            xmin =
                xmin,
            xmax =
                xmax,
            ymin =
                ymin,
            ymax =
                ymax
        ),
        inherit.aes =
            FALSE,
        fill =
            NA,
        color =
            OUTLINE_COLOR,
        linewidth =
            OUTLINE_LINEWIDTH
    ) +

'''


text = (
    text[:insert_at]
    +
    outline_layers
    +
    text[insert_at:]
)


## ================================================================== ##
## 9. Insert the two colored tracks immediately before completeness
## ================================================================== ##

complete_heading = re.search(
    r'    ## -+ ##\n'
    r'    ## 4\. Completeness\n'
    r'    ## -+ ##',
    text
)


if not complete_heading:

    raise SystemExit(
        "ERROR: could not find completeness track heading."
    )


track_layers = '''    ## -------------------------------------------------------------- ##
    ## 4. Cluster_00121
    ##
    ## Presence only. Absence remains white / empty.
    ## -------------------------------------------------------------- ##

    geom_rect(
        data =
            c121_ring,
        aes(
            xmin =
                xmin,
            xmax =
                xmax,
            ymin =
                ymin,
            ymax =
                ymax,
            fill =
                Cluster00121_plot
        ),
        inherit.aes =
            FALSE,
        color =
            NA
    ) +


    scale_fill_manual(

        values =
            CLUSTER00121_COLOR,

        breaks =
            "Present",

        labels =
            "Present",

        drop =
            FALSE,

        name =
            "Cluster_00121"
    ) +


    ggnewscale::new_scale_fill() +


    ## -------------------------------------------------------------- ##
    ## 5. Methane monooxygenase repertoire
    ##
    ## All 650 genomes.
    ##
    ## Particulate MMO = complete canonical pMMO OR complete pXMO.
    ## -------------------------------------------------------------- ##

    geom_rect(
        data =
            mmo_ring,
        aes(
            xmin =
                xmin,
            xmax =
                xmax,
            ymin =
                ymin,
            ymax =
                ymax,
            fill =
                MMO_status
        ),
        inherit.aes =
            FALSE,
        color =
            NA
    ) +


    scale_fill_manual(

        values =
            MMO_COLORS,

        breaks =
            c(
                "Particulate MMO",
                "sMMO",
                "Both",
                "Not detected"
            ),

        labels =
            c(
                "Particulate MMO (pmo/pxm)",
                "sMMO",
                "Both",
                "Not detected"
            ),

        drop =
            FALSE,

        name =
            "Methane\\nmonooxygenase"
    ) +


    ggnewscale::new_scale_fill() +


    ## -------------------------------------------------------------- ##
    ## 6. Completeness
    ## -------------------------------------------------------------- ##'''


text = (
    text[:complete_heading.start()]
    +
    track_layers
    +
    text[complete_heading.end():]
)


## Renumber contamination comment only.

text = re.sub(
    r'## 5\. Contamination',
    '## 7. Contamination',
    text,
    count=1
)


## ================================================================== ##
## 10. Add report output
## ================================================================== ##

report_marker = '''cat(
    "\\nMCA0421 positive genomes:\\n"
)
'''


if report_marker in text:

    report_pos = text.find(
        report_marker
    )

    ## Find insertion point after the following print() block by using
    ## the next known heading / cat statement.

    next_marker = '''cat(
    "\\nExample MAG species labels:\\n"
)
'''


    next_pos = text.find(
        next_marker,
        report_pos
    )


    if next_pos != -1:

        report_addition = '''cat(
    "\\nCluster_00121 positive genomes:\\n"
)

print(
    sum(
        ann$Cluster00121_plot ==
            "Present",
        na.rm = TRUE
    )
)


cat(
    "\\nMethane monooxygenase repertoire:\\n"
)

print(
    table(
        ann$MMO_status,
        useNA = "ifany"
    )
)


'''


        text = (
            text[:next_pos]
            +
            report_addition
            +
            text[next_pos:]
        )


## ================================================================== ##
## Save
## ================================================================== ##

path.write_text(
    text
)


print()
print("=" * 72)
print("FINAL TREE PATCH COMPLETE")
print("=" * 72)
print()
print("Added:")
print("  Cluster_00121")
print("  methane monooxygenase repertoire")
print()
print("New ring order:")
print("  Cyc2")
print("  Module 20 / Cluster_00035")
print("  MCA0421")
print("  Cluster_00121")
print("  methane monooxygenase")
print("  completeness")
print("  contamination")
print()
print("UNCHANGED:")
print("  tree")
print("  genus shading")
print("  all existing colors")
print("  all existing labels")
print("  MAG bolding")
print("  figure size")
print("  first-ring position")
print()
print("Backup:")
print(f"  {backup}")
print()
