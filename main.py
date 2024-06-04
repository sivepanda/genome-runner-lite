import requests
import pybedtools
import os

# API https://genome.ucsc.edu/goldenPath/help/api.html
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
        for line in response.text.split("\n")[1:]:
            split_line = line.strip().split("\t")
            if(len(split_line) > 1):
                res += "\t".join(*[split_line[1:7]]) + "\n"
        return res 
    else:
        response.raise_for_status()

def fetch_tracks(genome, features_of_interest, chrom, start, end, table="knownGene"):
    for track in features_of_interest:
        bed_data = fetch_bed_from_ucsc(genome, track, chrom, start, end)
        fi_nm = track + "_" + genome + ".bed"
        file = open("./track/" + fi_nm, "w")
        file.write(bed_data)
        file.close()


# Testing with human genome 38 and a few tracks of interest on chromosome 1 
genome = "hg38"
chrom = "chr1"
start = 60825584
end = 222387434
features_of_interest = ["wgEncodeRegDnaseClustered" ,  "tRNAs" , "knownAlt" , "cpgIslandExt"]
# bed_data = fetch_tracks(genome, features_of_interest, chrom, start, end)

features_of_interest_bedtools = [] 
for feature in features_of_interest: 
    features_of_interest_bedtools.append(pybedtools.example_bedtool(os.path.join(os.getcwd(),"track", feature + "_" + genome + ".bed"))) 


reference = pybedtools.example_bedtool(os.path.join(os.getcwd(), "track", "lincrna_tucp_1.bed"))

overlaps = []

for feature in features_of_interest_bedtools:
    overlap = feature.intersect(reference)
    for a in overlap:
        print(a)
        print(a.length)
    overlaps.append(sum(f.length for f in overlap) / len(reference))


print(overlaps)
