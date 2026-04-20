#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import argparse


# ----------------------------
# Functions
# ----------------------------

def len_from_CIGAR(string):
    out = re.split(r'(?<=[a-zA-Z])(?=.)', string)
    runs = []
    for i in out:
        if set('MND') & set(i):
            numbers = re.findall(r'\d+', i)
            integers = [int(num) for num in numbers][0]
            runs.append(integers)
    return sum(runs)

def feature_by_coor(int,feature1,end1,feature2,end2,feature3,end3):
    if int < end1:
        feature=feature1
    elif end1 + 1 < int < end2:
        feature=feature2
    elif end2 + 1 < int < end3:
        feature=feature3

    elif int == end1:
        feature="".join([feature1,'_3*'])
    elif int == end2:
        feature="".join([feature2,'_3*'])
    elif int == end3:
        feature="".join([feature3,'_3*'])

    elif int == 1:
        feature="".join([feature1,'_5*'])
    elif int == end1+ 1:
        feature="".join([feature2,'_5*'])
    elif int == end2 + 1:
        feature="".join([feature3,'_5*'])
    else:
        print("ERROR- coordinate outside range")
        exit(1)
    return(feature)

def plot_3_5(df, col_5PRIME, col_3PRIME, ref_len, HHcoor, HDVcoor, title, dpi=300): #Plot 5' and 3' ends

    # Define bins: one per base (0–240 → 240 bins)
    bins = np.arange(0, ref_len + 1, 1)

    # Compute histograms
    hist_5, bin_edges = np.histogram(df[col_5PRIME], bins=bins)
    hist_3, _ = np.histogram(df[col_3PRIME], bins=bins)

    # Bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Plot
    plt.figure(figsize=(10, 6))

    # 5' ends (positive)
    plt.bar(bin_centers, hist_5, width=1.0,
            color='green', alpha=0.7, label="5' PRIME ENDS")

    # 3' ends (negative)
    plt.bar(bin_centers, -hist_3, width=1.0,
            color='red', alpha=0.7, label="3' PRIME ENDS")

    # add lines denoting differnt parts of transcript
    plt.axvline(x=HHcoor+1, linestyle='--', linewidth=1.5,
                color='purple', label="5'HH ribozyme cut")
    plt.axvline(x=HDVcoor+1, linestyle='--', linewidth=1.5,
                color='blue', label="3' HDV ribozyme cut")
    
    #Axis formatting
    plt.axhline(0, color='black', linewidth=1)
    plt.xlim(0, ref_len)
    plt.xticks(np.arange(0, ref_len + 1, 10), rotation=45, ha='right')
    yticks = plt.yticks()[0] # Make y-axis labels positive
    plt.yticks(yticks, [abs(int(y)) for y in yticks])
    
    #label Axes
    plt.title(title)
    plt.xlabel("Position (bp)")
    plt.ylabel("Read Count")
    plt.legend()
    plt.tight_layout()

    #save as png
    outfile = f"{title}.png"
    plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Saved: {outfile}")

def plot_depth(df, col_depth, ref_len, HHcoor, HDVcoor, title, dpi=300):

    df = df.sort_values("POSITION")

    x = df["POSITION"]
    y = df[col_depth]

    plt.figure(figsize=(10, 4))

    plt.bar(x, y, width=1.0,
            color='blue', alpha=0.7, label="Read Depths")

    plt.axvline(x=HHcoor+0.5, linestyle='--', linewidth=1.5,
                color='purple', label="5'HH ribozyme cut")

    plt.axvline(x=HDVcoor+0.5, linestyle='--', linewidth=1.5,
                color='blue', label="3' HDV ribozyme cut")

    plt.axhline(0, color='black', linewidth=1)
    plt.xlim(0, ref_len)

    plt.xticks(np.arange(0, ref_len + 1, 10), rotation=45, ha='right')

    plt.ylim(0, y.max())

    yticks = plt.yticks()[0]
    plt.yticks(yticks, [abs(y) for y in yticks])

    plt.title(title)
    plt.xlabel("Position (bp)")
    plt.ylabel("Read Count")
    plt.legend()
    plt.tight_layout()

    outfile = f"{title}.png"
    plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
    plt.close()

    print(f"[INFO] Saved: {outfile}")



def plot_3_5_norm(df, col_5PRIME, col_3PRIME, ref_len, HHcoor, HDVcoor, title, dpi=300):

    df = df.sort_values("POSITION")

    x = df["POSITION"]
    y_5 = df[col_5PRIME]
    y_3 = -df[col_3PRIME]

    plt.figure(figsize=(10, 6))

    plt.bar(x, y_5, width=1.0,
            color='green', alpha=0.7, label="5' PRIME ENDS")

    plt.bar(x, y_3, width=1.0,
            color='red', alpha=0.7, label="3' PRIME ENDS")

    plt.axvline(x=HHcoor+0.5, linestyle='--', linewidth=1.5,
                color='purple', label="5'HH ribozyme cut")

    plt.axvline(x=HDVcoor+0.5, linestyle='--', linewidth=1.5,
                color='blue', label="3' HDV ribozyme cut")

    plt.axhline(0, color='black', linewidth=1)
    plt.xlim(0, ref_len)

    plt.xticks(np.arange(0, ref_len + 1, 10), rotation=45, ha='right')

    max_val = max(y_5.max(), abs(y_3.min()))
    plt.ylim(-max_val, max_val)

    yticks = plt.yticks()[0]
    plt.yticks(yticks, [abs(y) for y in yticks])

    plt.title(title)
    plt.xlabel("Position (bp)")
    plt.ylabel("Normalized Read Count")
    plt.legend()
    plt.tight_layout()

    outfile = f"{title}.png"
    plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
    plt.close()

    print(f"[INFO] Saved: {outfile}")


# ----------------------------
# Main
# ----------------------------

def main():

    parser = argparse.ArgumentParser(description="Plot 5' and 3' read end distributions from SAM files")

    parser.add_argument("-s", "--sam", required=True, help="Input SAM file")
    parser.add_argument("-d", "--depth", required=True, help="Depth file")

    parser.add_argument("--ref_len", type=int, default=240)
    parser.add_argument("--hh", type=float, default=59)
    parser.add_argument("--hdv", type=float, default=136)

    args = parser.parse_args()

    file = args.sam

    # count SAM header lines
    count = 0
    with open(file, "r") as f:
        for line in f:
            if line.startswith("@"):
                count += 1
            else:
                break

    sam_cols = ["QNAME", "SEQNAME", "5PRIME", "MAPQ", "CIGAR", "SEQ"]

    # Import sam file
    sam = pd.read_csv(
        file,
        sep="\t",
        skiprows=count, #gets ride of sam header
        names=sam_cols,
        usecols=[0, 2, 3, 4, 5, 9]
    ).sort_values("5PRIME")

    # Calculate 3' End from read start coordinate and CIGAR string. 
    # Only works for merged/single end reads running sense direction.
    sam["ALIGN_LEN"] = sam["CIGAR"].apply(len_from_CIGAR)
    sam["3PRIME"] = sam["5PRIME"] + sam["ALIGN_LEN"] - 1

    #Make and export table with fragment counts
    sam['5PRIME_feature']=sam['5PRIME'].apply(lambda x: feature_by_coor(x,'A_HH',args.hh,'B_tRNA',args.hdv,'C_HDV',args.ref_len))
    sam['3PRIME_feature']=sam['3PRIME'].apply(lambda x: feature_by_coor(x,'A_HH',args.hh,'B_tRNA',args.hdv,'C_HDV',args.ref_len))
    f_title="_".join([sam.at[1,'SEQNAME'],'fragmentcounts.csv'])
    fragment_counts=pd.DataFrame(sam.groupby(by=['5PRIME_feature','3PRIME_feature']).size()).reset_index(drop=False)
    fragment_counts.columns=['5PRIME_feature','3PRIME_feature','COUNT']
    fragment_counts.to_csv(f_title)
    print(f'"[INFO] Saved: {f_title}')

    # Make histogram of 5' and 3' ends
    title = f"{sam.at[1, 'SEQNAME']}_end_distribution"
    plot_3_5(
        sam,
        "5PRIME",
        "3PRIME",
        args.ref_len,
        args.hh,
        args.hdv,
        title
    )

    # Reformat data in new tables to allow easy normalization
    hist5 = sam.groupby("5PRIME").size().reset_index(name="5PRIMECOUNT").rename(columns = {"5PRIME":"POSITION"}) #5' end counts per position
    hist3 =sam.groupby(by='3PRIME').size().reset_index(name="3PRIMECOUNT").rename(columns = {"3PRIME":"POSITION"}) # 3' end counts per position
    hist = pd.merge(hist5, hist3,on="POSITION", how="outer") # merge tables
    hist.fillna(0, inplace=True) # deal with NA by making 0

    #Import depth information
    depth = pd.read_csv(args.depth, sep="\t", names=["SEQNAME", "POSITION", "DEPTH"])
    hist_depth = pd.merge(hist, depth, on="POSITION", how="outer").fillna(0)
    hist_depth.to_csv("_".join([sam.at[1,'SEQNAME'],'EndsDepths.csv']))
    
    # avoid division by zero
    hist_depth["DEPTH"] = hist_depth["DEPTH"].replace(0, np.nan)
    hist_depth["NORM5PRIME"] =hist_depth["5PRIMECOUNT"] / hist_depth["DEPTH"]
    hist_depth["NORM3PRIME"] = hist_depth["3PRIMECOUNT"] / hist_depth["DEPTH"]
    hist_depth = hist_depth.fillna(0)
    if hist_depth[['NORM5PRIME','NORM3PRIME']].isin([np.inf, -np.inf]).any().any():
        print('WARNING: end counts exceed depth')
        exit(1)


    # plot depths
    depth_title= f"{sam.at[1, 'SEQNAME']}_ReadDepth"

    plot_depth(
        hist_depth, 
        "DEPTH", 
        args.ref_len,
        args.hh,
        args.hdv,
        depth_title
    )
    

    # plot ends normalized to depths
    norm_title = f"{sam.at[1, 'SEQNAME']}_DepthNormalizedEnds"

    plot_3_5_norm(
        hist_depth,
        "NORM5PRIME",
        "NORM3PRIME",
        args.ref_len,
        args.hh,
        args.hdv,
        norm_title
    )


if __name__ == "__main__":
    main()
