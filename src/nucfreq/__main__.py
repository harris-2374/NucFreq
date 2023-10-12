"""
Author: Andrew Harris
Python 3.8

python nucFreqWholeGenome.py -c Oge-1_final.sizes.bed -b /WMLab3/ajharris/BigCatHiFiReads/ccs_bam_files/GxD2/Minimap2/Oge1HiFi_to_Oge1Ref.sort.bam -w 1000000 -o /WMLab3/ajharris/BigCatHiFiReads/ccs_bam_files/GxD2/NucFreq/Oge1_output -r /WMLab3/ajharris/BigCatHiFiReads/ccs_bam_files/GxD2/NucFreq/Oge1_RepeatMasker/rm_chroms/

python nucFreqWholeGenome.py -c Fca-126_final.sizes.bed -b /WMLab3/ajharris/BigCatHiFiReads/ccs_bam_files/GxD2/Minimap2/Fca126HiFi_to_Fca126Ref.sort.bam  -w 1000000 -o /WMLab3/ajharris/BigCatHiFiReads/ccs_bam_files/GxD2/NucFreq/Fca126_output -r /WMLab3/ajharris/BigCatHiFiReads/ccs_bam_files/GxD2/NucFreq/Fca126_RepeatMasker/rm_chroms/
"""
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
):
    """Run NucFreq commands per-window"""
    chromosome, length = row.chromosome, row.length
    range_groups = mygrouper(
        5, [i for i in range(WINDOWSIZE, length + WINDOWSIZE, WINDOWSIZE)]
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
        opng = png_dir / f"{chromosome}_{n}.png"
        try:
            if opng.exists():
                continue
            elif rm_files:
                rm_file = [i for i in rm_files if chromosome in i.stem][0]
                subprocess.run(
                    [
                        f"python NucPlot.py -y {YAXIS} -a -t 5 -r {rm_file} --bed {tmpbed} --obed {obed} {BAM} {opng}"
                    ],
                    stderr=subprocess.DEVNULL,
                    shell=True,
                    check=True,
                )
            else:
                subprocess.run(
                    [
                        f"python NucPlot.py -y {YAXIS} -a -t 5 --bed {tmpbed} --obed {obed} {BAM} {opng}"
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
    exit()
    # --------------------------------------------------------------
    # -- Imports --
    import math
    from PIL import Image, ImageEnhance

    # -- Helper Functions --
    def dimensions():
        width, height = DIMENSIONS.split("x")
        try:
            height = int(height)
        except ValueError:
            if width == "n":
                height = math.ceil(len(files) / 2)
            else:
                height = math.ceil(len(files) / int(width))
        try:
            width = int(width)
        except ValueError:
            if height == "n":
                width = math.ceil(len(files) / 2)
            else:
                width = math.ceil(len(files) / int(height))
        return width, height

    def build_file_df(files):
        df = pd.DataFrame(
            pd.NA,
            columns=["files", "chromosome", "start", "stop"],
            index=range(len(files)),
        )
        df["files"] = [f for f in files]
        df["chromosome"] = df["files"].apply(
            lambda x: x.split("/")[-1].strip(".png").split("_")[0]
        )
        df["start"] = df["files"].apply(
            lambda x: x.split("/")[-1].strip(".png").split("_")[1]
        )
        df["stop"] = df["files"].apply(
            lambda x: x.split("/")[-1].strip(".png").split("_")[2]
        )
        return df

    # -- Logic --
    files = sorted(
        [str(f) for f in png_dir.iterdir()],
        key=lambda x: (
            x.split("/")[-1].strip(".png").split("_")[0],
            int(x.split("/")[-1].strip(".png").split("_")[1]),
        ),
    )
    # files = sorted([str(f) for f in png_dir.iterdir() if chromosome in f.name])
    df = build_file_df(files)
    df = df[df["chromosome"] == chromosome]
    collaged_png_files = []
    png_file_groups = mygrouper(5, df["files"])
    filenum = 1
    # -- Generate collaged files --
    for f in png_file_groups:
        if len(files) < 100:
            outfile = collage_dir / f"{chromosome}_{filenum:03d}.png"
            collaged_png_files.append(outfile)
        elif 100 < len(files) < 1000:
            outfile = collage_dir / f"{chromosome}_{filenum:04d}.png"
            collaged_png_files.append(outfile)
        else:
            outfile = collage_dir / f"{chromosome}_{filenum:05d}.png"
            collaged_png_files.append(outfile)
        width, height = dimensions()
        collage = Image.new(
            "RGBA", (1500 * width, 500 * height), color=(255, 255, 255, 255)
        )
        c = 0
        try:
            for j in range(0, 500 * height, 500):
                for i in range(0, 1500 * width, 1500):
                    file = f[c]
                    photo = Image.open(file).convert("RGBA")
                    photo = photo.resize((1500, 500))
                    filter = ImageEnhance.Sharpness(photo)
                    filtered_photo = filter.enhance(1.0)
                    collage.paste(filtered_photo, (i, j))
                    c += 1
            collage.save(outfile)
            filenum += 1
            continue
        except:
            collage.save(outfile)
            filenum += 1
            continue
    # -- Merge all chromosome png files into a single PDF --
    img_list = []
    for f in collaged_png_files[1:]:
        img_list.append(Image.open(f).convert("RGB"))
        continue
    pdfoutfile = pdf_dir / f"{chromosome}.pdf"
    img1 = Image.open(collaged_png_files[0]).convert("RGB")
    img1.save(pdfoutfile, save_all=True, append_images=img_list)
    return


########################### Main Function ###########################
def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-c",
        "--chromLengths",
        type=Path,
        action="store",
        required=True,
        help="Bed file with chromosome lengths",
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
        help="Output directory",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        action="store",
        default=8,
        help="Number of threads",
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

    collage_dir.mkdir(parents=True, exist_ok=True)
    datapoints_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(parents=True, exist_ok=True)
    pdf_dir.mkdir(parents=True, exist_ok=True)
    # -- Load Chromosome Lengths --
    chrom_lengths = pd.read_csv(
        args.chromLengths,
        sep="\t",
        header=None,
        names=["chromosome", "start", "length"],
    )
    chrom_lengths["length"] = chrom_lengths["length"].apply(
        lambda x: roundUp(x, args.windowsize)
    )
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
        ),
        [r for r in chrom_lengths.itertuples(index=False)],
        **{"num_cpus": args.threads},
    )
    # -- Remove intermediate png files --
    shutil.rmtree(png_dir)
    shutil.rmtree(collage_dir)

    return


if __name__ == "__main__":
    sys.exit(main())
