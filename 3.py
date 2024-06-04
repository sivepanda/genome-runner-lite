import requests

genome = "hg38"
chrom = "chr1"
start = 60825584
end = 222387434
features_of_interest = ["wgEncodeRegDnaseClustered" ,  "tRNAs" , "knownAlt" , "knownGene" , "cpgIslandExt"]

def fetch_bed_from_ucsc(genome, track, chrom, start, end, table="knownGene"):
    base_url = "http://genome.ucsc.edu/cgi-bin/hgTables"
    params = {
        "db": genome,
        "hgta_group": "allTracks",
        "hgta_track": track,
        "hgta_table": track,
        "hgta_regionType": "range",
        "position": f"{chrom}:{start}-{end}",
        "hgta_outFileName": "stdout",
        "hgta_outputType": "primaryTable",
        "hgta_doTopSubmit": "get output"
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        res = ""
        for line in response.text.split("\n"):
            split_line = line.strip().split("\t")
            if(len(split_line) > 1):
                res += "\t".join([split_line[1], *split_line[3:7]]) + "\n"
        return res
    else:
        response.raise_for_status()


# print(fetch_bed_from_ucsc(genome, features_of_interest[0], chrom, start, end))

print(fetch_bed_from_ucsc(genome, features_of_interest[0], chrom, start, end).split("\n")[10])
print(fetch_bed_from_ucsc(genome, features_of_interest[1], chrom, start, end).split("\n")[10])
