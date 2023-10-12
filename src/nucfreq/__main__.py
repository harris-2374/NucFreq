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


########################## Helper Functions ##########################
def roundUp(x, WINDOWSIZE):
    """Round up length positions to nearest whole window"""
    return int(math.ceil(x / WINDOWSIZE)) * WINDOWSIZE


def mygrouper(n, iterable):
    args = [iter(iterable)] * n
    return list([e for e in t if e != None] for t in itertools.zip_longest(*args))


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
    nucfreqpath,
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
                subprocess.run(
                    [
                        f"python {nucfreqpath} -y {YAXIS} -a -t 5 -r {rm_file} --bed {tmpbed} --obed {obed} {BAM} {opng}"
                    ],
                    stderr=subprocess.DEVNULL,
                    shell=True,
                    check=True,
                )
            else:
                subprocess.run(
                    [
                        f"python {nucfreqpath} -y {YAXIS} -a -t 5 --bed {tmpbed} --obed {obed} {BAM} {opng}"
                    ],
                    stderr=subprocess.DEVNULL,
                    shell=True,
                    check=True,
                )
        except subprocess.CalledProcessError as e:
            print(e)
            exit()
        # Delete tmp file
        os.remove(tmpbed)
        continue
    # --------------------------------------------------------------
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
        default=200,
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
            nucfreqpath=args.nucplot,
        ),
        [list(r) for r in chrom_lengths.itertuples(index=False)],
        **{"num_cpus": 8},
        # **{"num_cpus": args.threads},
    )
    # -- Remove intermediate png files --
    shutil.rmtree(png_dir)
    return


if __name__ == "__main__":
    sys.exit(main())
