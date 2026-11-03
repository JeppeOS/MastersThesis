#!/usr/bin/env python3

from pathlib import Path
import math
import pandas as pd

from PIL import (
    Image,
    ImageDraw,
    ImageFont,
)


## ================================================================== ##
## PURPOSE
##
## Assemble the Stage-38 validated structure renders into four final
## thesis-ready figures.
##
## No scientific values are recalculated here.
##
## Each panel contains:
##
##   1. genome/taxon + source
##   2. protein-family annotation
##   3. pTM + US-align TM-score/RMSD relative to panel A
##
## Rendering:
##
##   - protein: AlphaFold pLDDT colors
##   - Heme C: dark red sticks
##   - Fe: orange/gold spheres
##
## The final footer explicitly states that TM-score/RMSD were obtained
## from protein-only US-align comparisons, while the displayed models
## are the original holo AlphaFold3 structures.
## ================================================================== ##


BASE = (
    Path.home()
    / "methanotrophs"
    / "methanotroph_project"
    / "jeppe"
)


MANIFEST = (
    BASE
    / "structure_standardization"
    / "03_plot_manifest"
    / "final_structure_plot_manifest.tsv"
)


PANEL_ROOT = (
    BASE
    / "structure_standardization"
    / "04_rendered"
    / "panels"
)


OUTDIR = (
    BASE
    / "structure_standardization"
    / "05_final_figures"
)


OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


## ================================================================== ##
## Input
## ================================================================== ##


df = pd.read_csv(
    MANIFEST,
    sep="\t"
)


if len(df) != 15:

    raise RuntimeError(
        f"Expected 15 structures; found {len(df)}."
    )


## ================================================================== ##
## Font handling
## ================================================================== ##


## ================================================================== ##
## Robust font discovery
##
## Try:
##   1. common system font paths
##   2. fonts installed inside the active Conda environment
##   3. recursive search for DejaVu / Liberation Sans
##
## This avoids assuming a specific Linux distribution.
## ================================================================== ##


import os


FONT_SEARCH_ROOTS = [

    Path(
        "/usr/share/fonts"
    ),

    Path(
        "/usr/local/share/fonts"
    ),

    Path.home()
    / ".fonts",

]


conda_prefix = os.environ.get(
    "CONDA_PREFIX"
)


if conda_prefix:

    FONT_SEARCH_ROOTS.extend(
        [

            Path(
                conda_prefix
            )
            / "fonts",

            Path(
                conda_prefix
            )
            / "share"
            / "fonts",

        ]
    )


REGULAR_NAMES = [

    "DejaVuSans.ttf",
    "LiberationSans-Regular.ttf",
    "LiberationSans.ttf",
    "Arial.ttf",
]


BOLD_NAMES = [

    "DejaVuSans-Bold.ttf",
    "LiberationSans-Bold.ttf",
    "Arial-Bold.ttf",
    "Arialbd.ttf",
]


def find_font(
    filenames
):

    ## -------------------------------------------------------------- ##
    ## First try direct recursive search under likely roots
    ## -------------------------------------------------------------- ##

    for root in FONT_SEARCH_ROOTS:

        if not root.exists():

            continue


        for filename in filenames:

            hits = list(
                root.rglob(
                    filename
                )
            )


            if hits:

                return hits[
                    0
                ]


    ## -------------------------------------------------------------- ##
    ## Generic fallback: look for any Sans regular/bold TTF
    ## -------------------------------------------------------------- ##

    want_bold = any(
        "Bold" in x
        or "bold" in x
        for x in filenames
    )


    for root in FONT_SEARCH_ROOTS:

        if not root.exists():

            continue


        for hit in root.rglob(
            "*.ttf"
        ):

            name = hit.name.lower()


            if "sans" not in name:

                continue


            is_bold = (
                "bold" in name
            )


            if is_bold == want_bold:

                return hit


    raise RuntimeError(
        "Could not locate a usable TrueType sans-serif font.\n"
        "Searched:\n"
        + "\n".join(
            str(x)
            for x in FONT_SEARCH_ROOTS
        )
    )


FONT_REGULAR = find_font(
    REGULAR_NAMES
)


FONT_BOLD = find_font(
    BOLD_NAMES
)


print(
    "Regular font:",
    FONT_REGULAR
)

print(
    "Bold font:",
    FONT_BOLD
)



def font(
    size,
    bold=False
):

    return ImageFont.truetype(
        str(
            FONT_BOLD
            if bold
            else FONT_REGULAR
        ),
        size=size
    )


## ================================================================== ##
## Layout
##
## Designed for full-width thesis placement.
## ================================================================== ##


FIG_W = 2300

LEFT = 80
RIGHT = 80

TOP = 65

COL_GAP = 65

TITLE_H = 70
SUBTITLE_H = 55

TITLE_TO_GRID = 18

PANEL_W = (
    FIG_W
    - LEFT
    - RIGHT
    - COL_GAP
) // 2


IMAGE_W = PANEL_W - 50

## Stage-38 images are 1600 × 1250.
IMAGE_H = int(
    IMAGE_W
    * 1250
    / 1600
)


TEXT_TOP_GAP = 6

LINE1_H = 52
LINE2_H = 48
LINE3_H = 46

PANEL_BOTTOM_PAD = 16

ROW_H = (
    IMAGE_H
    + TEXT_TOP_GAP
    + LINE1_H
    + LINE2_H
    + LINE3_H
    + PANEL_BOTTOM_PAD
)


ROW_GAP = 16

LEGEND_TOP_GAP = 8
LEGEND_H = 235

BOTTOM = 55


## ================================================================== ##
## Typography
## ================================================================== ##


FONT_TITLE = font(
    44,
    bold=True
)


FONT_SUBTITLE = font(
    27
)


FONT_LINE1 = font(
    28,
    bold=True
)


FONT_LINE2 = font(
    25
)


FONT_LINE3 = font(
    23
)


FONT_LEGEND_TITLE = font(
    29,
    bold=True
)


FONT_LEGEND = font(
    26
)


FONT_NOTE = font(
    19
)


## ================================================================== ##
## Colours
## ================================================================== ##


BG = (
    255,
    255,
    255
)


TEXT = (
    0,
    0,
    0
)


SUBTEXT = (
    60,
    60,
    60
)


LIGHT_GREY = (
    225,
    225,
    225
)


PLDDT = [

    (
        "Very high (>90)",
        (
            0,
            83,
            214
        )
    ),

    (
        "Confident (70–90)",
        (
            101,
            203,
            243
        )
    ),

    (
        "Low (50–70)",
        (
            255,
            219,
            19
        )
    ),

    (
        "Very low (<50)",
        (
            255,
            125,
            69
        )
    ),
]


HEME_RED = (
    184,
    26,
    26
)


FE_GOLD = (
    255,
    171,
    0
)


## ================================================================== ##
## Text helpers
## ================================================================== ##


def text_width(
    draw,
    text,
    selected_font
):

    box = draw.textbbox(
        (
            0,
            0
        ),
        text,
        font=selected_font
    )


    return (
        box[
            2
        ]
        -
        box[
            0
        ]
    )


def fitted_font(
    draw,
    text,
    max_width,
    initial_size,
    minimum_size,
    bold=False
):

    size = initial_size


    while size >= minimum_size:

        candidate = font(
            size,
            bold=bold
        )


        if text_width(
            draw,
            text,
            candidate
        ) <= max_width:

            return candidate


        size -= 1


    return font(
        minimum_size,
        bold=bold
    )


def draw_centered(
    draw,
    text,
    y,
    selected_font,
    fill=TEXT
):

    width = text_width(
        draw,
        text,
        selected_font
    )


    draw.text(

        (
            (
                FIG_W
                -
                width
            )
            / 2,
            y
        ),

        text,

        font=selected_font,

        fill=fill
    )



## ================================================================== ##
## Shared whitespace cropping
##
## For each figure set, determine the bounding box containing visible
## non-white structure pixels across ALL panels.
##
## The union of those boxes is then used for every panel in that set.
##
## This preserves common structure scale/orientation while removing
## white space between the protein and its label.
## ================================================================== ##


SHARED_CROP = {}


def visible_bbox(
    image,
    threshold=248
):

    rgb = image.convert(
        "RGB"
    )


    pixels = rgb.load()


    xs = []
    ys = []


    for y in range(
        rgb.height
    ):

        for x in range(
            rgb.width
        ):

            r, g, b = pixels[
                x,
                y
            ]


            if (
                r < threshold
                or g < threshold
                or b < threshold
            ):

                xs.append(
                    x
                )

                ys.append(
                    y
                )


    if not xs:

        return (
            0,
            0,
            rgb.width,
            rgb.height
        )


    return (
        min(
            xs
        ),
        min(
            ys
        ),
        max(
            xs
        )
        + 1,
        max(
            ys
        )
        + 1
    )


for figure_set, group in df.groupby(
    "figure_set",
    sort=False
):

    boxes = []


    for _, row in group.iterrows():

        path = (
            PANEL_ROOT
            / figure_set
            / (
                f"{row['panel']}.png"
            )
        )


        image = Image.open(
            path
        ).convert(
            "RGB"
        )


        boxes.append(
            visible_bbox(
                image
            )
        )


    left = min(
        b[
            0
        ]
        for b in boxes
    )


    top = min(
        b[
            1
        ]
        for b in boxes
    )


    right = max(
        b[
            2
        ]
        for b in boxes
    )


    bottom = max(
        b[
            3
        ]
        for b in boxes
    )


    ## -------------------------------------------------------------- ##
    ## Add some breathing room so the structure does not touch edges.
    ## -------------------------------------------------------------- ##

    PAD_X = 55

    PAD_TOP = 45

    ## Smaller bottom padding is deliberate:
    ## this is what brings the annotation closer to the structure.
    PAD_BOTTOM = 20


    sample_path = (
        PANEL_ROOT
        / figure_set
        / (
            f"{group.iloc[0]['panel']}.png"
        )
    )


    sample = Image.open(
        sample_path
    )


    left = max(
        0,
        left - PAD_X
    )


    right = min(
        sample.width,
        right + PAD_X
    )


    top = max(
        0,
        top - PAD_TOP
    )


    bottom = min(
        sample.height,
        bottom + PAD_BOTTOM
    )


    SHARED_CROP[
        figure_set
    ] = (
        left,
        top,
        right,
        bottom
    )


    print(
        "Shared crop:",
        figure_set,
        SHARED_CROP[
            figure_set
        ]
    )



## ================================================================== ##
## Image helper
## ================================================================== ##


def load_panel(
    figure_set,
    panel
):

    path = (
        PANEL_ROOT
        / figure_set
        / f"{panel}.png"
    )


    if not path.is_file():

        raise RuntimeError(
            f"Missing rendered panel:\n{path}"
        )


    image = Image.open(
        path
    ).convert(
        "RGB"
    )


    ## -------------------------------------------------------------- ##
    ## Apply ONE shared crop box to every panel in the figure set.
    ## -------------------------------------------------------------- ##

    image = image.crop(
        SHARED_CROP[
            figure_set
        ]
    )


    ## -------------------------------------------------------------- ##
    ## Resize based on cropped aspect ratio.
    ##
    ## Width remains standardized across all panels.
    ## -------------------------------------------------------------- ##

    new_h = int(
        IMAGE_W
        * image.height
        / image.width
    )


    image = image.resize(

        (
            IMAGE_W,
            new_h
        ),

        Image.Resampling.LANCZOS
    )


    return image


## ================================================================== ##
## Panel positions
##
## 2 panels:
##
##   A B
##
## 4 panels:
##
##   A B
##   C D
##
## 5 panels:
##
##   A B
##   C D
##    E
##
## E is centered in the final row.
## ================================================================== ##


def panel_positions(
    n_panels,
    grid_top,
    row_height
):

    positions = []


    if n_panels not in (
        2,
        4,
        5
    ):

        raise RuntimeError(
            f"Unsupported panel count: {n_panels}"
        )


    for i in range(
        n_panels
    ):

        row = i // 2

        col = i % 2


        if (
            n_panels == 5
            and i == 4
        ):

            x = (
                FIG_W
                -
                PANEL_W
            ) // 2

        else:

            x = (
                LEFT
                +
                col
                * (
                    PANEL_W
                    +
                    COL_GAP
                )
            )


        y = (
            grid_top
            +
            row
            * (
                row_height
                +
                ROW_GAP
            )
        )


        positions.append(
            (
                x,
                y
            )
        )


    return positions


## ================================================================== ##
## Draw one panel
## ================================================================== ##


def draw_panel(
    canvas,
    draw,
    row,
    x,
    y
):

    panel_img = load_panel(
        row[
            "figure_set"
        ],
        row[
            "panel"
        ]
    )


    image_x = (
        x
        +
        (
            PANEL_W
            -
            IMAGE_W
        )
        // 2
    )


    canvas.paste(
        panel_img,
        (
            image_x,
            y
        )
    )


    text_y = (
        y
        +
        panel_img.height
        +
        TEXT_TOP_GAP
    )


    max_text_width = (
        PANEL_W
        -
        24
    )


    ## ------------------------------------------------------------- ##
    ## Line 1
    ## ------------------------------------------------------------- ##


    line1 = str(
        row[
            "display_line_1"
        ]
    )


    line1_font = fitted_font(

        draw,
        line1,
        max_text_width,
        initial_size=35,
        minimum_size=25,
        bold=True
    )


    draw.text(

        (
            x
            +
            12,
            text_y
        ),

        line1,

        font=line1_font,

        fill=TEXT
    )


    ## ------------------------------------------------------------- ##
    ## Line 2
    ## ------------------------------------------------------------- ##


    line2 = str(
        row[
            "display_line_2"
        ]
    )


    line2_font = fitted_font(

        draw,
        line2,
        max_text_width,
        initial_size=31,
        minimum_size=23
    )


    draw.text(

        (
            x
            +
            12,
            text_y
            +
            LINE1_H
        ),

        line2,

        font=line2_font,

        fill=TEXT
    )


    ## ------------------------------------------------------------- ##
    ## Line 3
    ## ------------------------------------------------------------- ##


    line3 = str(
        row[
            "display_line_3"
        ]
    )


    line3_font = fitted_font(

        draw,
        line3,
        max_text_width,
        initial_size=28,
        minimum_size=21
    )


    draw.text(

        (
            x
            +
            12,
            text_y
            +
            LINE1_H
            +
            LINE2_H
        ),

        line3,

        font=line3_font,

        fill=SUBTEXT
    )


## ================================================================== ##
## Legend
## ================================================================== ##


def draw_legend(
    draw,
    legend_top
):

    draw.line(

        (
            LEFT,
            legend_top,
            FIG_W - RIGHT,
            legend_top
        ),

        fill=LIGHT_GREY,

        width=2
    )


    y_title = (
        legend_top
        +
        22
    )


    draw.text(

        (
            LEFT,
            y_title
        ),

        "AlphaFold pLDDT confidence",

        font=FONT_LEGEND_TITLE,

        fill=TEXT
    )


    ## ------------------------------------------------------------- ##
    ## pLDDT row
    ## ------------------------------------------------------------- ##


    y_items = (
        y_title
        +
        42
    )


    x = LEFT


    SWATCH = 28
    TEXT_GAP = 10
    ITEM_GAP = 43


    for label, colour in PLDDT:

        draw.rectangle(

            (
                x,
                y_items,
                x + SWATCH,
                y_items + SWATCH
            ),

            fill=colour,

            outline=(
                100,
                100,
                100
            ),

            width=1
        )


        draw.text(

            (
                x
                +
                SWATCH
                +
                TEXT_GAP,
                y_items
                -
                1
            ),

            label,

            font=FONT_LEGEND,

            fill=TEXT
        )


        label_w = text_width(
            draw,
            label,
            FONT_LEGEND
        )


        x += (
            SWATCH
            +
            TEXT_GAP
            +
            label_w
            +
            ITEM_GAP
        )


    ## ------------------------------------------------------------- ##
    ## Heme legend
    ## ------------------------------------------------------------- ##


    y_heme = (
        y_items
        +
        52
    )


    x = LEFT


    draw.line(

        (
            x,
            y_heme + 14,
            x + 42,
            y_heme + 14
        ),

        fill=HEME_RED,

        width=8
    )


    draw.text(

        (
            x + 55,
            y_heme
        ),

        "Heme C ligand",

        font=FONT_LEGEND,

        fill=TEXT
    )


    x += 250


    r = 13


    draw.ellipse(

        (
            x,
            y_heme + 1,
            x + 2 * r,
            y_heme + 1 + 2 * r
        ),

        fill=FE_GOLD,

        outline=(
            120,
            90,
            0
        )
    )


    draw.text(

        (
            x + 38,
            y_heme
        ),

        "Fe",

        font=FONT_LEGEND,

        fill=TEXT
    )


    ## ------------------------------------------------------------- ##
    ## Method note
    ## ------------------------------------------------------------- ##


    note = (
        "Displayed structures are mature holo AlphaFold3 models. "
        "TM-score and RMSD are from protein-only US-align comparisons "
        "against panel A; Heme C ligands were excluded from alignment."
    )


    note_font = fitted_font(

        draw,
        note,
        FIG_W
        -
        LEFT
        -
        RIGHT,
        initial_size=23,
        minimum_size=18
    )


    draw.text(

        (
            LEFT,
            y_heme
            +
            50
        ),

        note,

        font=note_font,

        fill=SUBTEXT
    )


## ================================================================== ##
## Assemble each figure
## ================================================================== ##


for figure_set, group in df.groupby(
    "figure_set",
    sort=False
):

    group = group.sort_values(
        "panel"
    ).reset_index(
        drop=True
    )


    n_panels = len(
        group
    )


    n_rows = math.ceil(
        n_panels
        / 2
    )


    grid_top = (
        TOP
        +
        TITLE_H
        +
        SUBTITLE_H
        +
        TITLE_TO_GRID
    )


    ## ------------------------------------------------------------- ##
    ## Calculate actual displayed image height after shared crop.
    ##
    ## load_panel() applies the same crop to every panel within this
    ## figure set, so panel A gives the correct height for the set.
    ## ------------------------------------------------------------- ##

    sample_panel = load_panel(
        figure_set,
        group.iloc[
            0
        ][
            "panel"
        ]
    )


    actual_image_h = sample_panel.height


    actual_row_h = (
        actual_image_h
        + TEXT_TOP_GAP
        + LINE1_H
        + LINE2_H
        + LINE3_H
        + PANEL_BOTTOM_PAD
    )


    grid_h = (
        n_rows
        * actual_row_h
        +
        (
            n_rows
            -
            1
        )
        * ROW_GAP
    )


    legend_top = (
        grid_top
        +
        grid_h
        +
        LEGEND_TOP_GAP
    )


    FIG_H = (
        legend_top
        +
        LEGEND_H
        +
        BOTTOM
    )


    canvas = Image.new(

        "RGB",

        (
            FIG_W,
            FIG_H
        ),

        BG
    )


    draw = ImageDraw.Draw(
        canvas
    )


    ## ------------------------------------------------------------- ##
    ## Title
    ## ------------------------------------------------------------- ##


    title = str(
        group[
            "figure_display_title"
        ].iloc[
            0
        ]
    )


    subtitle = str(
        group[
            "figure_display_subtitle"
        ].iloc[
            0
        ]
    )


    title_font = fitted_font(

        draw,
        title,
        FIG_W - LEFT - RIGHT,
        initial_size=50,
        minimum_size=36,
        bold=True
    )


    subtitle_font = fitted_font(

        draw,
        subtitle,
        FIG_W - LEFT - RIGHT,
        initial_size=32,
        minimum_size=24
    )


    draw_centered(

        draw,
        title,
        TOP,
        title_font,
        fill=TEXT
    )


    draw_centered(

        draw,
        subtitle,
        TOP
        +
        TITLE_H,
        subtitle_font,
        fill=SUBTEXT
    )


    ## ------------------------------------------------------------- ##
    ## Panels
    ## ------------------------------------------------------------- ##


    positions = panel_positions(
        n_panels,
        grid_top,
        actual_row_h
    )


    for (
        (_, row),
        (
            x,
            y
        )
    ) in zip(
        group.iterrows(),
        positions
    ):

        draw_panel(

            canvas,
            draw,
            row,
            x,
            y
        )


    ## ------------------------------------------------------------- ##
    ## Legend/footer
    ## ------------------------------------------------------------- ##


    draw_legend(
        draw,
        legend_top
    )


    ## ------------------------------------------------------------- ##
    ## Save
    ## ------------------------------------------------------------- ##


    out_png = (
        OUTDIR
        / (
            f"{figure_set}_"
            "standardized_structure_figure.png"
        )
    )


    out_pdf = (
        OUTDIR
        / (
            f"{figure_set}_"
            "standardized_structure_figure.pdf"
        )
    )


    canvas.save(

        out_png,

        dpi=(
            300,
            300
        )
    )


    canvas.save(

        out_pdf,

        resolution=300
    )


    print()

    print(
        "=" * 100
    )

    print(
        figure_set
    )

    print(
        "=" * 100
    )

    print(
        f"Panels: {n_panels}"
    )

    print(
        f"PNG: {out_png}"
    )

    print(
        f"PDF: {out_pdf}"
    )


print()

print(
    "=" * 100
)

print(
    "STAGE 39 COMPLETE"
)

print(
    "=" * 100
)

print()

print(
    f"Final figures:\n{OUTDIR}"
)

print()

print(
    "No scientific values were recalculated during assembly."
)
