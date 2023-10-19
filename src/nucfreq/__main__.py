import argparse
import itertools
import math
import os
import subprocess
import shutil
import sys
from functools import partial
from pathlib import Path

## Dependencies
import pandas as pd
from p_tqdm import p_umap
import numpy as np
import pysam
import re
import pandas as pd
import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns


########################## Helper Functions ##########################
def roundUp(x, WINDOWSIZE):
    """Round up length positions to nearest whole window"""
    return int(math.ceil(x / WINDOWSIZE)) * WINDOWSIZE


def mygrouper(n, iterable):
    args = [iter(iterable)] * n
    return list([e for e in t if e != None] for t in itertools.zip_longest(*args))


def nucplot(
    infile=None,
    outfile=None,
    d=False,
    legend=False,
    zerostart=False,
    a=False,
    repeatmasker=None,
    regions=[],
    bed=None,
    obed=None,
    minobed=2,
    ylim=None,
    font_size=16,
    freey=False,
    height=4,
    width=16,
    dpi=600,
    threads=8,
    header=False,
    psvsites=None,
    soft=False,
    minclip=1000,
):
    M = 0  # M  BAM_CMATCH      0
    I = 1  # I  BAM_CINS        1
    D = 2  # D  BAM_CDEL        2
    N = 3  # N  BAM_CREF_SKIP   3
    S = 4  # S  BAM_CSOFT_CLIP  4
    H = 5  # H  BAM_CHARD_CLIP  5
    P = 6  # P  BAM_CPAD        6
    E = 7  # =  BAM_CEQUAL      7
    X = 8  # X  BAM_CDIFF       8
    B = 9  # B  BAM_CBACK       9
    NM = 10  # NM       NM tag  10
    conRef = [M, D, N, E, X]  # these ones "consume" the reference
    conQuery = [M, I, S, E, X]  # these ones "consume" the query
    conAln = [M, I, D, N, S, E, X]  # these ones "consume" the alignments

    # sys.stderr.write("Packages loaded\n")

    def getSoft(read, group=0):
        rtn = []
        cigar = read.cigartuples
        start = cigar[0]
        end = cigar[-1]
        if start[0] in [S, H]:
            rtn.append(
                (
                    read.reference_name,
                    "start",
                    start[1],
                    read.reference_start,
                    read,
                    group,
                )
            )
        if end[0] in [S, H]:
            rtn.append(
                (read.reference_name, "end", end[1], read.reference_end, read, group)
            )
        return rtn

    soft = []
    bam = pysam.AlignmentFile(infile, threads=threads)
    refs = {}
    regions = []
    if regions is not None or bed is not None:
        # sys.stderr.write("Reading in the region or bed argument(s).\n")
        if regions is not None:
            for region in regions:
                match = re.match("(.+):(\d+)-(\d+)", region)
                assert match, region + " not valid!"
                chrm, start, end = match.groups()
                refs[chrm] = [int(start), int(end)]
                regions.append((chrm, int(start), int(end)))

        if bed is not None:
            for line in open(bed):
                if line[0] == "#":
                    continue
                line = line.strip().split()
                chrm, start, end = line[0:3]
                refs[chrm] = [int(start), int(end)]
                regions.append((chrm, int(start), int(end)))

    else:
        for read in bam.fetch(until_eof=True):
            ref = read.reference_name
            # read.query_qualities = [60] * len(read.query_sequence)
            soft += getSoft(read)
            if ref not in refs:
                if a:
                    refs[ref] = [0, 2147483648]
                else:
                    refs[ref] = [2147483648, 0]

            start = read.reference_start
            end = read.reference_end
            if refs[ref][0] > start:
                refs[ref][0] = start
            if refs[ref][1] < end:
                refs[ref][1] = end
        for contig in refs:
            regions.append((contig, refs[contig][0], refs[contig][1]))

    def getCovByBase(contig, start, end):
        # coverage = bam.count_coverage(contig, quality_threshold=None, start=start, stop=end, read_callback="nofilter")
        coverage = bam.count_coverage(
            contig,
            start=start,
            stop=end,
            read_callback="nofilter",
            quality_threshold=None,
        )
        assert len(coverage) == 4
        cov = {}
        cov["A"] = coverage[0]
        cov["C"] = coverage[1]
        cov["T"] = coverage[2]
        cov["G"] = coverage[3]
        return cov

    #
    # creates a table of nucleotide frequescies
    #
    # nf = []
    nf = {"contig": [], "position": [], "A": [], "C": [], "G": [], "T": [], "group": []}
    GROUPS = 0
    for contig, start, end in regions:
        cov = getCovByBase(contig, start, end)
        contiglen = len(cov["A"])
        if contiglen > 0:
            nf["contig"] += [contig] * contiglen
            nf["group"] += [GROUPS] * contiglen
            nf["position"] += list(range(start, start + contiglen))
            nf["A"] += cov["A"]
            nf["C"] += cov["C"]
            nf["G"] += cov["G"]
            nf["T"] += cov["T"]
            if soft:
                for read in bam.fetch(contig, start, end):
                    soft += getSoft(read, group=GROUPS)
            GROUPS += 1

    df = pd.DataFrame(nf)
    sort = np.flip(np.sort(df[["A", "C", "G", "T"]].values), 1)
    df["first"] = sort[:, 0]
    df["second"] = sort[:, 1]
    df["third"] = sort[:, 2]
    df["fourth"] = sort[:, 3]
    df.sort_values(by=["contig", "position", "second"], inplace=True)

    soft = pd.DataFrame(
        soft, columns=["contig", "side", "value", "position", "read", "group"]
    )
    soft = soft[soft.value >= minclip]
    soft.sort_values(by=["contig", "position"], inplace=True)

    RM = None
    colors = sns.color_palette()
    cmap = {}
    counter = 0
    if repeatmasker is not None:
        names = [
            "score",
            "perdiv",
            "perdel",
            "perins",
            "qname",
            "start",
            "end",
            "left",
            "strand",
            "repeat",
            "family",
            "rstart",
            "rend",
            "rleft",
            "ID",
        ]
        lines = []
        for idx, line in enumerate(repeatmasker):
            if idx > 2:
                lines.append(line.strip().split()[0:15])

        RM = pd.DataFrame(lines, columns=names)
        RM.start = RM.start.astype(int)
        RM.end = RM.end.astype(int)
        RM["label"] = RM.family.str.replace("/.*", "")
        for idx, lab in enumerate(sorted(RM.label.unique())):
            cmap[lab] = colors[counter % len(colors)]
            counter += 1
        RM["color"] = RM.label.map(cmap)

        repeatmasker.close()

    # SET up the plot based on the number of regions
    HEIGHT = GROUPS * height
    # set text size
    matplotlib.rcParams.update({"font.size": font_size})
    # make axes
    fig, axs = plt.subplots(nrows=GROUPS, ncols=1, figsize=(width, HEIGHT))
    if GROUPS == 1:
        axs = [axs]
    # make space for the bottom label of the plot
    # fig.subplots_adjust(bottom=0.2)
    # set figure YLIM
    YLIM = int(max(df["first"]) * 1.05)

    # iterate over regions
    counter = 0
    for group_id, group in df.groupby(by="group"):
        if freey:
            YLIM = int(max(group["first"]) * 1.05)

        contig = list(group.contig)[0]

        truepos = group.position.values
        first = group["first"].values
        second = group["second"].values

        # df = pd.DataFrame(nf, columns=["contig", "position", "A", "C", "G", "T"])
        if obed:
            tmp = group.loc[
                group.second >= minobed,
                ["contig", "position", "position", "first", "second"],
            ]
            if counter == 0:
                tmp.to_csv(
                    obed,
                    header=["#contig", "start", "end", "first", "second"],
                    sep="\t",
                    index=False,
                )
            else:
                tmp.to_csv(obed, mode="a", header=None, sep="\t", index=False)

        # get the correct axis
        ax = axs[group_id]

        if RM is not None:
            rmax = ax
            rm = RM[
                (RM.qname == contig)
                & (RM.start >= min(truepos))
                & (RM.end <= max(truepos))
            ]
            assert len(rm.index) != 0, "No matching RM contig"
            # rmax.set_xlim(rm.start.min(), rm.end.max())
            # rmax.set_ylim(-1, 9)
            rmlength = len(rm.index) * 1.0
            rmcount = 0
            rectangles = []
            height_offset = max(YLIM, ylim) / 20
            for idx, row in rm.iterrows():
                width = row.end - row.start
                rect = patches.Rectangle(
                    (row.start, max(YLIM, ylim) - height_offset),
                    width,
                    height_offset,
                    linewidth=1,
                    edgecolor="none",
                    facecolor=row.color,
                    alpha=0.75,
                )
                rmax.add_patch(rect)
                # rectangles.append(rect)

            # rmax.add_collection(
            #    PatchCollection(
            #        rectangles,
            # match_original=True
            # facecolor=rm.color,
            # edgecolor='none'
            #    ))
            plt.show()

        (prime,) = ax.plot(
            truepos,
            first,
            "o",
            color="black",
            markeredgewidth=0.0,
            markersize=2,
            label="most frequent base pair",
        )
        (sec,) = ax.plot(
            truepos,
            second,
            "o",
            color="red",
            markeredgewidth=0.0,
            markersize=2,
            label="second most frequent base pair",
        )

        # inter = int( (max(truepos)-min(truepos))/50)
        # sns.lineplot(  (truepos/inter).astype(int)*inter, first, ax = ax, err_style="bars")

        maxval = max(truepos)
        minval = min(truepos)
        subval = 0

        title = "{}:{}-{}\n".format(contig, minval, maxval)
        if GROUPS > 1:
            ax.set_title(title, fontweight="bold")
        # sys.stderr.write(title)

        if zerostart:
            subval = minval - 1
            ax.set_xticks(
                [x for x in ax.get_xticks() if (x - subval > 0) and (x < maxval)]
            )
            maxval = maxval - minval

        if maxval < 1000000:
            xlabels = [format((label - subval), ",.0f") for label in ax.get_xticks()]
            lab = "bp"
        elif maxval < 10000000:
            xlabels = [
                format((label - subval) / 1000, ",.1f") for label in ax.get_xticks()
            ]
            lab = "kbp"
        else:
            xlabels = [
                format((label - subval) / 1000, ",.1f") for label in ax.get_xticks()
            ]
            lab = "kbp"
            # xlabels = [format( (label-subval)/1000000, ',.2f') for label in ax.get_xticks()]
            # lab = "Mbp"

        if ylim is not None:
            ax.set_ylim(0, ylim)
        else:
            ax.set_ylim(0, YLIM)

        ax.set_xlabel("Assembly position ({})".format(lab), fontweight="bold")
        ax.set_ylabel("Sequence read depth", fontweight="bold")

        # Including this causes some internal bug in matplotlib when the font-size changes
        # ylabels = [format(label, ",.0f") for label in ax.get_yticks()]
        # ax.set_yticklabels(ylabels)
        ax.set_xticklabels(xlabels)

        # Hide the right and top spines
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")

        if counter == 0 and legend:
            lgnd = plt.legend(loc="upper right")
            for handle in lgnd.legendHandles:
                handle._sizes = [300.0]

        if not soft.empty:
            tmpsoft = soft[soft.group == group_id]
            if len(tmpsoft) > 0:
                axsoft = ax.twinx()
                axsoft.invert_yaxis()
                bins = width * 5
                color = "darkgreen"
                sns.distplot(
                    tmpsoft.position,
                    bins=bins,
                    kde=False,
                    ax=axsoft,
                    hist_kws={
                        "weights": tmpsoft.value / 1000,
                        "alpha": 0.25,
                        "color": color,
                    },
                )
                bot, top = axsoft.get_ylim()
                axsoft.set_ylim(1.1 * bot, 0)
                axsoft.set_xlim(minval, maxval)
                # Hide the right and top spines
                axsoft.spines["top"].set_visible(False)
                # Only show ticks on the left and bottom spines
                axsoft.yaxis.set_ticks_position("right")
                # axsoft.xaxis.set_ticks_position('bottom')
                axsoft.tick_params(axis="y", colors=color)
                axsoft.set_ylabel("Clipped Bases (kbp)", color=color)

        if psvsites is not None:  # and len(psvsites)>0):
            cuts = {}
            for idx, line in enumerate(open(psvsites).readlines()):
                try:
                    vals = line.strip().split()
                    cuts[idx] = list(map(int, vals))
                    # make plot
                    x = np.array(cuts[idx]) - 1
                    idxs = np.isin(truepos, x)
                    y = second[idxs]
                    ax.plot(x, y, alpha=0.5)  # , label="group:{}".format(idx) )
                except Exception as e:
                    # print("Skipping because error: {}".format(e), file=sys.stderr)
                    continue

        # outpath = os.path.abspath(outfile)
        # if(counter == 0):
        #  outf =   outpath
        # else:
        #  name, ext = os.path.splitext(outpath)
        #  outf = "{}_{}{}".format(name, counter + 1, ext)

        counter += 1

    plt.tight_layout()
    plt.savefig(outfile, dpi=dpi)
    return


def run(
    row,
    png_dir,
    collage_dir,
    pdf_dir,
    datapoints_dir,
    rm_files,
    BAM,
    WINDOWSIZE,
    OUTPUT,
    DIMENSIONS,
    YAXIS,
    conv2pdf,
):
    """Run NucFreq commands per-window"""
    chromosome, length = row[0], row[1]
    chrompngdir = png_dir / chromosome
    chrompngdir.mkdir(parents=True, exist_ok=True)
    range_groups = mygrouper(
        3, [i for i in range(WINDOWSIZE, length + WINDOWSIZE, WINDOWSIZE)]
    )
    for n, group in enumerate(range_groups):
        # 1. Make tmp bed file with coordinates
        tmpbed = OUTPUT / f"{chromosome}.{n}.tmp.bed"
        try:
            os.remove(tmpbed)
        except:
            pass
        with open(tmpbed, "w") as oh:
            windows = "\n".join([f"{chromosome}\t{i-WINDOWSIZE+1}\t{i}" for i in group])
            oh.write(windows)
            pass
        # 2. Run command
        obed = datapoints_dir / f"{chromosome}_{n}.datapoints.bed"
        colorf = datapoints_dir / f"{chromosome}_{n}.colors"
        opng = chrompngdir / f"{chromosome}_{n}.png"
        try:
            if opng.exists():
                continue
            elif rm_files:
                rm_file = [i for i in rm_files if chromosome in i.stem][0]
                nucplot(
                    infile=BAM,
                    outfile=opng,
                    ylim=YAXIS,
                    threads=5,
                    repeatmasker=rm_file,
                    bed=tmpbed,
                    a=True,
                    obed=obed,
                )
            else:
                nucplot(
                    infile=BAM,
                    outfile=opng,
                    ylim=YAXIS,
                    threads=5,
                    bed=tmpbed,
                    a=True,
                    obed=obed,
                )
        except subprocess.CalledProcessError as e:
            print(e)
            return
        # Delete tmp file
        os.remove(tmpbed)
        continue
    # --------------------------------------------------------------
    if conv2pdf:
        # -- Imports --
        from PIL import Image

        # -- Logic --
        files = sorted(
            [str(f) for f in chrompngdir.iterdir()],
            key=lambda x: (
                x.split("/")[-1].strip(".png").split("_")[0],
                int(x.split("/")[-1].strip(".png").split("_")[1]),
            ),
        )

        images = [Image.open(f) for f in files]

        pdf_path = pdf_dir / f"{chromosome}.pdf"

        images[0].save(
            pdf_path, "PDF", resolution=100.0, save_all=True, append_images=images[1:]
        )
    return


########################### Main Function ###########################
def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-f",
        "--fai",
        type=Path,
        action="store",
        required=True,
        help="Reference assembly .fai file",
    )
    parser.add_argument(
        "-b",
        "--bam",
        type=Path,
        action="store",
        required=True,
        help="Input BAM file",
    )
    parser.add_argument(
        "-r",
        "--rm",
        type=Path,
        action="store",
        default=None,
        help="Directory containing RepeatMasker .out files to add to plots",
    )
    parser.add_argument(
        "-w",
        "--windowsize",
        type=int,
        action="store",
        required=True,
        default=100000,
        help="Window size (default: 100000)",
    )
    parser.add_argument(
        "-d",
        "--dimensions",
        action="store",
        default="1x4",
        help="Figure dimensions. Give as nxn or 3xn or 3x4 - WIDTHxHEIGHT",
        metavar="",
        type=str,
    )
    parser.add_argument(
        "-y",
        "--yaxis",
        action="store",
        default=100,
        help="Y-axis limit for NucFreq plots",
        metavar="",
        type=int,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        action="store",
        default=Path().cwd(),
        help="Output directory (default: cwd)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        action="store",
        default=8,
        help="Number of threads (default: 8)",
    )
    parser.add_argument(
        "--nucplot",
        type=Path,
        action="store",
        default="src/nucfreq/NucPlot.py",
        help="Nucplot path",
    )
    parser.add_argument(
        "--conv2pdf",
        action="store_true",
        default=False,
        help="Convert individual png files into a single PDF file",
    )
    args = parser.parse_args()
    # -- Create output directories --
    png_dir = args.output / "PNG"
    datapoints_dir = args.output / "datapoints"
    collage_dir = args.output / "PNGcollage"
    pdf_dir = args.output / "PDFcollage"

    # collage_dir.mkdir(parents=True, exist_ok=True)
    datapoints_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(parents=True, exist_ok=True)
    pdf_dir.mkdir(parents=True, exist_ok=True)
    # -- Load Chromosome Lengths --
    chrom_lengths = pd.read_csv(
        args.fai,
        sep="\t",
        header=None,
    )
    chrom_lengths[1] = chrom_lengths[1].apply(lambda x: roundUp(x, args.windowsize))
    # -- Run chromosomes --
    p_umap(
        partial(
            run,
            png_dir=png_dir,
            collage_dir=collage_dir,
            pdf_dir=pdf_dir,
            datapoints_dir=datapoints_dir,
            rm_files=Path(args.rm).iterdir() if args.rm else None,
            BAM=args.bam,
            WINDOWSIZE=args.windowsize,
            OUTPUT=args.output,
            DIMENSIONS=args.dimensions,
            YAXIS=args.yaxis,
            conv2pdf=args.conv2pdf,
        ),
        [list(r) for r in chrom_lengths.itertuples(index=False)],
        **{"num_cpus": args.threads},
    )
    # -- Remove intermediate png files --
    if args.conv2pdf:
        shutil.rmtree(png_dir)
    return


if __name__ == "__main__":
    sys.exit(main())
