import shutil
import sys
import statistics
import os
from glob import glob
import subprocess
import pandas as pd
import numpy as np
import tempfile
from matplotlib import pyplot as plt
from time import sleep as sl
"""
reqs:
samtools
viralconsensus
bwa-mem2
seqkit
pandas
"""
def rundepth(bamfile,region,chroms):
    results = {}
    for chrom in chroms:

        with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmpbed:
            with open(region) as f:
                for line in f:
                    if line.startswith(chrom + "\t"):
                        tmpbed.write(line)
            tmpbed_name = tmpbed.name

        sorted_bam = bamfile.replace(".bam", ".sorted.bam")
        if not os.path.exists(sorted_bam):
            res = subprocess.run(["samtools", "sort", "-o", sorted_bam, bamfile], capture_output=True, text=True)
            if res.returncode != 0:
                print(f"Error running samtools sort: {res.stderr}")
                sys.exit(1)
            res2 = subprocess.run(["samtools", "index", sorted_bam], capture_output=True, text=True)
            if res2.returncode != 0:
                print(f"Error running samtools index: {res2.stderr}")
                sys.exit(1)

        view = subprocess.Popen(
            ["samtools", "view", "-b", "-f", "0x2", sorted_bam, chrom],
            stdout=subprocess.PIPE
        )
        depth_proc = subprocess.Popen(
            ["samtools", "depth", "-b", tmpbed_name, "-"],
            stdin=view.stdout,
            stdout=subprocess.PIPE,
            text=True
        )

        count = 0
        total = 0
        values = []

        for line in depth_proc.stdout:
            parts = line.split()
            if len(parts) < 3:
                continue

            d = int(parts[2])

            total += d
            count += 1
            values.append(d)

        if count > 0:
            mean = total / count
            median = statistics.median(values)
        else:
            mean = 0
            median = 0

        results[chrom] = {"mean": mean, "median": median}

    return results

def runstats(bamfile, region=False, min_mapq=0, min_maplen=0):

    filters = ["proper_pair"]

    if min_mapq > 0:
        filters.append(f"mapping_quality >= {min_mapq}")

    if min_maplen > 0:
        filters.append(f"sequence_length >= {min_maplen}")

    filter_expr = " and ".join(filters)

    sambamba_cmd = [
        "sambamba",
        "view",
        "-f", "bam",
        "-F", filter_expr
    ]

    if region:
        sambamba_cmd += ["-L", region]

    sambamba_cmd.append(bamfile)

    view = subprocess.Popen(
        sambamba_cmd,
        stdout=subprocess.PIPE
    )

    result = subprocess.run(
        ["samtools", "stats"],
        stdin=view.stdout,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"Error running samtools stats: {result.stderr}")
        sys.exit(1)

    return result

def get_inset_size(bamfile,enrichedregions,unenrichedregions,chroms):
    """samtools view -b your.bam chr1:100000-200000 | samtools sort - | samtools stats | grep ^IS"""

    sizes,unfilteredsizes = get_fragment_sizes(bamfile)
    enrichedresults =  rundepth(bamfile,enrichedregions,chroms)
    unenrichedresults  = rundepth(bamfile, unenrichedregions,chroms)

    if not sizes:
        outstats = {}
        for chrom in chroms:
            outstats[chrom] = ("NA","NA","NA")
        return outstats,[],[]

    outstats = {}
    for chrom in chroms:
        cenrichedmedian_cov = enrichedresults[chrom]["median"]
        cunenrichedmedian_cov = unenrichedresults[chrom]["median"]
        if cunenrichedmedian_cov < 1:
            ratio = (cenrichedmedian_cov + 1) / (cunenrichedmedian_cov + 1)
        else:
            ratio = cenrichedmedian_cov / cunenrichedmedian_cov
        outstats[chrom] = (enrichedresults[chrom]["median"],unenrichedresults[chrom]["median"],ratio)


    return outstats,sizes,unfilteredsizes

def read_lengths(fq_file):
    result = subprocess.run(
        ["seqkit", "fx2tab", "-l", fq_file],
        capture_output=True, text=True
    )
    lengths = []
    for line in result.stdout.strip().split("\n"):
        if not line.strip():  # skip empty lines
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        lengths.append(int(len(parts[1])))

    return lengths
def get_fragment_sizes(bamfile):

    cmd = [
        "samtools", "view",
        "-f", "0x2",
        "-F", "0x900",
        bamfile
    ]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    sizes = []
    unfiltered = []

    for line in proc.stdout:

        parts = line.split("\t")

        tlen = abs(int(parts[8]))
        mapq = int(parts[4])
        seqlen = len(parts[9])

        unfiltered.append(tlen)

        if mapq > 20 and seqlen > 50 and tlen > 50:
            sizes.append(tlen)

    return sizes, unfiltered

def read_lengths_stream(fq_file):

    cmd = ["seqkit", "fx2tab", "-l", fq_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    count = 0
    total = 0
    lengths = []

    for line in proc.stdout:
        parts = line.rstrip().split("\t")
        if len(parts) < 2:
            continue

        l = len(parts[1])
        total += l
        count += 1
        lengths.append(l)  # only needed for median

    avg = total / count if count else 0
    median = statistics.median(lengths) if lengths else 0

    return avg, median, count,lengths

def stats(fq1_file,fq2_file):
    r1avg,r1median,r1count,r1lengths = read_lengths_stream(fq1_file)
    r2avg,r2median,r2count,r2lengths = read_lengths_stream(fq2_file)
    bothcounts = r1count + r2count
    bothlengths = r1lengths + r2lengths
    avg_len = statistics.mean(bothlengths) if bothlengths else 0
    median_len = statistics.median(bothlengths) if bothcounts else 0
    r1avg = statistics.mean(r1lengths) if r1lengths else 0
    r2avg = statistics.mean(r2lengths) if r2lengths else 0
    return bothlengths, avg_len, median_len, r1avg, r1median, r2avg, r2median,r1lengths,r2lengths

def get_args():
    import argparse
    parser = argparse.ArgumentParser(description="Gather QC statistics from tNGS data")
    parser.add_argument("-s", "--sample", help="sample id, (must be read name minus suffix)",
                        required=True)
    parser.add_argument("-r","--raw", help="forward and reverse raw reads. Space delimited.",nargs=2,required=False)
    parser.add_argument("-t","--trimmed", help="forward and reverse trimmed reads. Space delimited.",nargs=2,required=False)
    parser.add_argument("-k", "--kraken", help="kraken2 report file", required=False)
    parser.add_argument("--castanetbam", help="castanet output bam file", required=False)
    parser.add_argument("--castanetdepth", help="castanet output depths file", required=False)
    parser.add_argument("-b", "--batch", help="batch name", required=True)
    parser.add_argument("-m","--mttarget", help="fasta file of mt genome to map to",required=False)
    parser.add_argument("--bedfile", help="bed file of enriched, unenriched and avoided regions",required=False)
    parser.add_argument("-o","--outdir", help="output directory",required=True)
    parser.add_argument("--keeptmp", help="keep temporary files",action="store_true",default=False)
    parser.add_argument("--humantargets",help="fasta file of human targets included in capture panel",default=False)
    parser.add_argument("--humanmappref",help="castanet mapping reference table for human targets ",default=False)
    parser.add_argument("--threads", help="number of threads to use", default=1)
    args = parser.parse_args()
    return args

def run_bwa_mt(r1,r2,mttarget,outbam):
    cmd = f"bwa-mem2 mem -t 4 {mttarget} {r1} {r2} | samtools view -h -Sb -F4 -F2048 | samtools sort - 1> {outbam}"
    res = subprocess.run(cmd, shell=True,capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Error running bwa-mem2 or samtools: {res.stderr}")
        sys.exit(1)
    if os.path.exists(outbam):
        return outbam
    else:
        print(f"Error: {outbam} not created")
        sys.exit(1)


def summarize_kraken_kingdoms(report_file,taxa):

    counts ={x:["0","0"] for x in taxa}
    total_reads = 0

    with open(report_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            parts = [p.strip() for p in parts]
            if len(parts) < 6:
                continue

            perc, clade_reads, direct_reads, rank, taxid, name = parts
            if name in taxa:
                counts[name] = [clade_reads,perc]
    return counts,total_reads

def get_regions(bedfile,outprefix=""):
    enrichedbed = outprefix + "enriched.bed"
    chroms = []
    with open(bedfile) as f, open(enrichedbed, "w") as out:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) > 2:
                chrom, start, end, category = line.strip().split("\t")[:4]
                if category == "enriched":
                    out.write(f"{chrom}\t{start}\t{end}\n")
                if chrom not in chroms:
                    chroms.append(chrom)

    unenrichedbed = outprefix + "unenriched.bed"
    with open(bedfile) as f, open(unenrichedbed, "w") as out:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) > 2:
                chrom, start, end, category = line.strip().split("\t")[:4]
                if category == "unenriched":
                    out.write(f"{chrom}\t{start}\t{end}\n")

    return enrichedbed,unenrichedbed,chroms

def insert_size_plot(sizes,title,outpath,unfilteredsizes=None):

    if unfilteredsizes:
        arr = np.array(unfilteredsizes)
        filtlarr = np.array(sizes)
    else:
        arr = np.array(sizes)
    min_val = int(arr.min())
    max_val = int(arr.max())
    start = (min_val // 10) * 10
    end = ((max_val // 10) + 1) * 10
    bin_edges = np.arange(start, end + 1, 10)

    xmax = max_val if max_val < 1200 else 1200

    plt.figure()
    if unfilteredsizes:
        if unfilteredsizes:
            plt.hist(np.minimum(arr, xmax), bins=np.arange(start, xmax + 10, 10), color="blue", edgecolor=None,
                     alpha=0.5)
            plt.hist(np.minimum(filtlarr, xmax), bins=np.arange(start, xmax + 10, 10), color="blue", edgecolor=None)
        else:
            plt.hist(arr, bins=np.arange(start, xmax + 10, 10), edgecolor="blue")
    plt.xlabel("Insert Size")
    plt.title(title)
    plt.xlim(0, xmax)
    plt.ylim(bottom=0.9)
    plt.tight_layout()
    plt.yscale("log")
    plt.savefig(outpath)
    plt.close()

def insertSizeStats(sizes):
    if len(sizes) == 0:
        return "NA","NA",["NA","NA","NA"]
    mean_val = statistics.mean(sizes)

    stdev_val = statistics.pstdev(sizes)
    quartiles = statistics.quantiles(sizes, n=4)
    quartiles = [f"{x:.2f}" for x in quartiles]
    if len(quartiles) < 3:
        quartiles = ["NA","NA","NA"]
    return f"{mean_val:.2f}",f"{stdev_val:.2f}",quartiles

def get_castanet_stats(depthfile):

    depthdata = pd.read_csv(depthfile, sep=",")
    if "n_reads_all" in depthdata.columns:
        total_reads = depthdata["n_reads_all"].sum()
    elif "reads_for_mapping" in depthdata.columns:
        total_reads = depthdata["reads_for_mapping"].sum()
    else:
        print("Error: neither n_reads_all nor reads_for_mapping column found in depth csv")
        return np.nan, np.nan
    dedup_reads = depthdata["n_reads_dedup"].sum()
    return total_reads,dedup_reads

def run_castanet(sample,outdir,trimmed,humantargets,threads,mapref):
    trimmedfolder = os.path.dirname(trimmed[0])
    if len(glob(trimmedfolder + "/" + "*")) != 2:
        sys.exit("to perform humantargets analysis with castanet reads must be in one folder per sample")
    castanetcmd= f"""python3 -m app.castanet_lite -ExpName {sample} -ExpDir {trimmedfolder} -SaveDir {outdir} -RefStem {humantargets} -DoKrakenPrefilter False -KrakenDbDir "" -DoConsensus False -DoTrimming False -NThreads {threads} -PostFilt False -MappingRefTable {mapref}
    """
    res = subprocess.run(castanetcmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Error running castanet: \nSTDOUT\n\n{res.stdout}\n\n\nSTDERR\n{res.stderr}")
        sys.exit(1)
    outdepth = outdir + "/" + sample + "/" + sample + "_depth.csv"
    depthplots = outdir + "/" + sample + "/Depth_output"
    if os.path.exists(outdepth):
        if not os.path.exists(outdir+"/plots"):
            os.mkdir(outdir+"/plots")
        for depthfile in glob(depthplots + "/human*"):
            shutil.copy(depthfile, outdir + "/plots/")
        return outdepth
    else:
        print(f"expected output depth file, {outdepth} was not found, check castanet outputs below: \nSTDOUT\n\n{res.stdout}\n\n\nSTDERR\n{res.stderr}")
        sys.exit(1)


def main():

    args = get_args()
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    outfile = args.outdir + "/" + args.sample + "_qc.csv"
    outf = open(outfile, "w")
    outstring = f"{args.batch},{args.sample}"
    outheader = f"batch,sampleid"

    if args.trimmed:
        bothlengths, avg_len, median_len, r1avg, r1median, r2avg, r2median, r1len, r2len = stats(args.trimmed[0],args.trimmed[1])
        lenbothlens = len(bothlengths)
        r1count=len(r1len)
        r2count=len(r2len)
    else:
        args.trimmed = ["NA", "NA"]
        bothlengths, avg_len, median_len, r1avg, r1median, r2avg, r2median, r1len, r2len = "NA","NA","NA","NA","NA","NA","NA",[],[]
        lenbothlens = "NA"
        r1count="NA"
        r2count="NA"
    print("Trimmed reads processed")
    if args.raw:
        rawbothlengths, _, _, _, _, _, _, _, _ = stats(args.raw[0], args.raw[1])
        lenrawbothlengths = len(rawbothlengths)
        print("Raw reads processed")
    else:
        args.raw = ["NA","NA"]
        lenrawbothlengths = "NA"


    outheader += ",R1,R2,trimmedR1,trimmedR2,rawreads,trimmedreads,avglen,medlen,r1count,r2count,r1avglen,r1medlen,r2avglen,r2medlen"
    outstring += f",{args.raw[0]},{args.raw[1]},{args.trimmed[0]},{args.trimmed[1]},{lenrawbothlengths},{lenbothlens},{avg_len},{median_len},{r1count},{r2count},{r1avg},{r1median},{r2avg},{r2median}"


    if args.trimmed[0] != "NA" and args.mttarget and args.bedfile:
        if not os.path.exists(args.mttarget) or not os.path.exists(args.bedfile):
            sys.exit("Please provide either args.bedfile and args.mttarget to perform enrichment analysis")
        mtbam = args.outdir + "/" + args.sample + ".mt_bwa.bam"
        run_bwa_mt(args.trimmed[0],args.trimmed[1],args.mttarget,args.outdir+"/"+args.sample+".mt_bwa.bam")
        print("BWA mapping done")
        enrichedbed,unenrichedbed,chroms = get_regions(args.bedfile,args.outdir+"/"+args.sample+"_")
        print("Enriched reads processed")
        enrichstats,sizes,unfilteredsizes = get_inset_size(mtbam,enrichedbed,unenrichedbed,chroms)
        # enrichstats,sizes = get_inset_size(mtbam,enrichedbed,unenrichedbed,chroms)
        # get_fragment_sizes(mtbam)
        mean_ins,stdev_ins, mtquartiles  =  insertSizeStats(sizes)
        mean_ins_uf,stdev_ins_uf, quartiles_uf = insertSizeStats(unfilteredsizes)
        title = f"{args.sample} COX1 insert size distribution" if args.sample else "Insert size distribution"
        outname = f"{args.sample}_COX1_insert_size_distribution.png" if args.sample else "insert_size_distribution.png"
        outpath = os.path.join(args.outdir, outname) if args.outdir else outname
        if sizes:
            insert_size_plot(sizes,title,outpath,unfilteredsizes=unfilteredsizes)
            # insert_size_plot(sizes,title,outpath)
        outstring += f",{mean_ins},{stdev_ins},{mtquartiles[0]},{mtquartiles[1]},{mtquartiles[2]}"
        outheader += ",enrichmentloci_avginsert,enrichmentloci_stdinsert,enrichmentloci_insert25,enrichmentloci_insert50,enrichmentloci_insert75"
        print("Insert sizes and enrichment processed")

        maxchrom,enriched,unenriched,enrichratio = "NA",0,0,0
        for chrom in chroms:
            if chrom in enrichstats:
                stenriched, stunenriched, stratio = enrichstats[chrom]
                if stenriched == "NA":
                    continue
                if stenriched > enriched:
                    maxchrom,enriched,unenriched,enrichratio = chrom,stenriched,stunenriched,stratio,
        if maxchrom == "NA":
            outstring += f",NA,NA,NA,NA"
            outheader += f",enrichedMedian,unenrichedMedian,enrichmentRatio,enrichmentloci"
        else:
            outstring += f",{enriched:.2f},{unenriched:.2f},{enrichratio:.2f},{maxchrom}"
            outheader += f",enrichedMedian,unenrichedMedian,enrichmentRatio,enrichmentloci"
    elif args.mttarget and args.trimmed[0] == "NA":
        sys.exit("to run enrichment please provide trimmed reads with --trimmed flag")
    else:
        outstring += f",NA,NA,NA,NA"
        outheader += f",enrichedMedian,unenrichedMedian,enrichmentRatio,enrichmentloci"
    if args.humantargets and args.humanmappref:
        human_depth = run_castanet(args.sample,args.outdir,args.trimmed,args.humantargets,args.threads,args.humanmappref)
        depthdf = pd.read_csv(human_depth)
        human = depthdf[depthdf["probetype"].str.contains("human")]
        if "batch" in human.columns:
            ind=["batch", "sampleid"]
        else:
            ind=["sampleid"]
        if "n_reads_all" in human.columns:
            allreads = "n_reads_all"
        elif "reads_for_mapping" in human.columns:
            allreads = "reads_for_mapping"
        else:
            sys.exit("Error: neither n_reads_all nor reads_for_mapping column found in castanet depth csv")
        humanall = human.pivot(
            index=ind,
            columns="probetype",
            values=allreads).fillna(0)
        humanall["allhuman"] = humanall.sum(axis=1)
        humanall.columns = [x + "_all" for x in humanall.columns]
        humandedup = human.pivot(
            index=ind,
            columns="probetype",
            values="n_reads_dedup").fillna(0)
        humandedup["allhuman"] = humandedup.sum(axis=1)
        humandedup.columns = [x + "_dedup" for x in humandedup.columns]
        humannc2 = human.pivot(
            index=ind,
            columns="probetype",
            values="prop_npos_cov2").fillna(0)
        humannc2.columns = [x + "_nc2" for x in humannc2.columns]
        humandedup.reset_index()
        humannc2.reset_index()



        allstats = pd.merge(humanall, humandedup, left_on=ind, right_on=ind)
        allstats = pd.merge(allstats, humannc2, left_on=ind, right_on=ind)
        outheaderls = [x for x in  allstats.columns if x not in ind]
        outstr = allstats[outheaderls].iloc[0].tolist()
        print("Human targets processed with castanet")
        outstring += "," + ",".join([str(x) for x in outstr])
        outheader += "," + ",".join([str(x) for x in outheaderls])

        shutil.rmtree(args.outdir+"/"+args.sample)
    elif args.humanmappref or args.humantargets:
        sys.exit(f"human target castanet analysis requires --humanmappref and --humantargets inputs")
    if args.kraken:
        taxa =["Eukaryota","Bacteria","Archaea","Viruses","Fungi","Caudoviricetes","Homo sapiens"]
        outheader += "," + ",".join([f"kraken:{x}" for x in taxa])
        if os.path.exists(args.kraken):
            krakencounts,total_reads = summarize_kraken_kingdoms(args.kraken,taxa)
            krakenres = ",".join([krakencounts[x][0] for x in taxa])
            outstring += f",{krakenres}"
        else:
            krakenres = ",".join(["NA" for x in taxa])
            outstring += f",{krakenres}"
        print("kraken results processed")
    if args.castanetbam:

        if len(args.castanetbam) == 0:
            sys.exit(f"Error: no bam file found in {args.castanet}")

        sizes,unfilteredsizes = get_fragment_sizes(args.castanetbam)
        if sizes:
            mean_val,stdev_val,quartiles = insertSizeStats(sizes)
            mean_valuf,stdev_uf,quartilesuf = insertSizeStats(unfilteredsizes)
            title = f"{args.sample} All reads insert size distribution" if args.sample else "All reads insert size distribution"
            outname = f"{args.sample}_all_reads_insert_size_distribution.png" if args.sample else "All_reads_insert_size_distribution.png"
            outpath = os.path.join(args.outdir, outname) if args.outdir else outname

            insert_size_plot(sizes, title, outpath,unfilteredsizes=unfilteredsizes)
            outheader += ",all_mapped_avginsert,all_mapped_stdinsert,all_mapped_insert25,all_mapped_insert50,all_mapped_insert75,filtered_mapped_avginsert,filtered_mapped_stdinsert,filtered_mapped_insert25,filtered_mapped_insert50,filtered_mapped_insert75"
            outstring += f",{mean_valuf},{stdev_uf},{quartilesuf[0]},{quartilesuf[1]},{quartilesuf[2]},{mean_val},{stdev_val},{quartiles[0]},{quartiles[1]},{quartiles[2]}"

        else:
            outheader += ",all_mapped__avginsert,all_mapped__stdinsert,all_mapped__insert25,all_mapped__insert50,all_mapped__insert75,filtered_mapped__avginsert,filtered_mapped__stdinsert,filtered_mapped__insert25,filtered_mapped__insert50,filtered_mapped__insert75"
            outstring += f",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA"

        ## get total mapped and total dedup reads from castanet


        if not args.keeptmp and args.mttarget:
            os.remove(mtbam)
            os.remove(unenrichedbed)
            os.remove(enrichedbed)
            os.remove(mtbam.replace(".bam", ".sorted.bam"))
            os.remove(mtbam.replace(".bam", ".sorted.bam.bai"))
            # os.remove(sortedbam)
            # os.remove(sortedbam + ".bai")
    if args.castanetdepth:
        if os.path.exists(args.castanetdepth):
            total_reads,dedup_reads = get_castanet_stats(args.castanetdepth)
            outheader += ",castanet_total_mapped_reads,castanet_dedup_reads"
            outstring += f",{total_reads},{dedup_reads}"
            print("castanet results processed")
        else:
            print("WARNING: castanet depth file not found, skipping castanet summaries")

    outf.write(f"{outheader}\n")
    outf.write(f"{outstring}\n")
    outf.close()

if __name__ == "__main__":
    main()
