import requests
import pybedtools
import os

# API https://genome.ucsc.edu/goldenPath/help/api.html

# TODO: Move functions into a local library and create runner file


# Fetch a single track as a BED file from UCSC
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

# Fetch ALL desired tracks from UCSC Table Browser
def fetch_tracks(genome, features_of_interest, chrom, start, end, table="knownGene"):
    print("Pulling data from UCSC Genome Browser...")
    for track in features_of_interest:
        bed_data = fetch_bed_from_ucsc(genome, track, chrom, start, end)
        fi_nm = track + "_" + genome + ".bed"
        file = open("./track/" + fi_nm, "w")
        file.write(bed_data)
        file.close()
    print("Completed pulling data.")


# Returns an array of bedtools objects for each of the features of interest. Pull data from UCSC if desired.
def create_bedtools(genomic_features, base, genomic_features_to_pull=[], pull_new_data=False):
    genomic_features_bedtools = []

    # TODO: Dynamically pull data if the bed file does not exist in ./track
    if pull_new_data == True:
        if len(genomic_features_to_pull) == 0:
            genomic_features_to_pull = genomic_features

        bed_data = fetch_tracks(genome, genomic_features_to_pull, chrom, start, end)

    print("Creating BEDTools objects...")
    for feature in genomic_features: 
        genomic_features_bedtools.append(pybedtools.example_bedtool(os.path.join(os.getcwd(), base , feature + "_" + genome + ".bed"))) 
    print("BEDTools created.")

    return genomic_features_bedtools

# Calculate the ratio of overlaps between a reference and an array of features of interest
def calculate_overlaps(reference, genomic_features_bedtools, minimum_overlap):
    overlaps = []
    print("Calculating overlaps...")

    for feature in genomic_features_bedtools:
        # TODO: Create shuffles to complete either chi-square or Monte Carlo significance tests
        
        overlap = feature.intersect(reference, f = minimum_overlap, u=True)
        overlaps.append(sum(f.length for f in overlap) / sum(f.length for f in feature))

    print("Completed calulcation.\n\n\n")
    return overlaps


# Testing with human genome 38 and a few tracks of interest on chromosome 1 
genome = "hg38"
chrom = "chr1"
start = 60825584
end = 222387434
genomic_features = ["wgEncodeRegDnaseClustered" ,  "tRNAs" , "knownAlt" , "cpgIslandExt" , "centromeres" , "lincrna_tucp"]
genomic_features_to_pull = genomic_features[:-1]
print(genomic_features_to_pull)
base = "track"
pull_new_data = True
min_overlap = 1

features_of_interest_bedtools = create_bedtools(genomic_features, base, genomic_features_to_pull, pull_new_data)

reference = pybedtools.example_bedtool(os.path.join(os.getcwd(), base, "cpgIslandExt_hg38.bed"))

overlaps = calculate_overlaps(reference, features_of_interest_bedtools, min_overlap)

# overlaps

print(overlaps)
