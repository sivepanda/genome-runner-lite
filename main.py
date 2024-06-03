import requests


# API https://genome.ucsc.edu/goldenPath/help/api.html
def fetch_bed_from_ucsc(genome, chrom, start, end, table="knownGene"):
    base_url = "http://genome.ucsc.edu/cgi-bin/hgTables"
    params = {
        "db": genome,
        "hgta_group": "allTracks",
        "hgta_track": table,
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

# Example usage
genome = "hg38"
chrom = "chr1"
start = 100000
end = 200000
bed_data = fetch_bed_from_ucsc(genome, chrom, start, end)
print(bed_data)
