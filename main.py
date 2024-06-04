import requests
import pybedtools

# API https://genome.ucsc.edu/goldenPath/help/api.html
def fetch_bed_from_ucsc(genome, track, chrom, start, end, table="knownGene"):
    base_url = "http://genome.ucsc.edu/cgi-bin/hgTables"
    params = {
        "db": genome,
        "hgta_group": "allTracks",
        "hgta_track": track,
        "hgta_table": table,
        "hgta_regionType": "range",
        "position": f"{chrom}:{start}-{end}",
        "hgta_outFileName": "stdout",
        "hgta_outputType": "primaryTable",
        "hgta_doTopSubmit": "get output"
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        return response.text
    else:
        response.raise_for_status()

def fetch_tracks(genome, tracks_of_interest, chrom, start, end, table="knownGene"):
    for track in tracks_of_interest:
        bed_data = fetch_bed_from_ucsc(genome, track, chrom, start, end)
        fi_nm = track + "_" + genome + ".bed"
        file = open("./track/" + fi_nm, "w")
        file.write(bed_data)
        file.close()


# Testing with human genome 38 and a few tracks of interest on chromosome 1 
genome = "hg38"
chrom = "chr1"
start = 100000
end = 200000
tracks_of_interest = ["wgEncodeRegDnaseClustered" , "wgEncodeRegTfbsClustered" ,  "gptIslandExt" , "oregannoswitchDbTss" , "tRNAs" , "knownAlt" , "knownGeneExons" , "knownGene" , "cpgIslandExt"]
bed_data = fetch_tracks(genome, tracks_of_interest, chrom, start, end)
