import requests
import random
import concurrent.futures
import numpy as np
from mpmath import mpf
from pull_data import create_bedtools
import pybedtools
import os
import threading

def progress_bar(amount, total, length=40):
    fraction = amount/total

    arrow = int(fraction * length - 1) * '-' + '>'
    padding = int(length - len(arrow)) * ' '

    ending = '\n' if amount == total else '\r'

    print(f'Progress: [{arrow}{padding}] {int(fraction*100)}%', end=ending)


def create_permutation(genome, reference, feature, minimum_overlap):
    # return sum(f.length for f in (feature.intersect(reference, f=minimum_overlap, u=True))) / (sum(f.length for f in feature))
    feature = feature.shuffle(genome=genome, chrom=True, seed=random.randint(1, 10000000))
    print(f"Thread {threading.get_ident()} is processing a permutation")
    return sum( f.length for f in  (feature.intersect(reference, f=minimum_overlap, u=True))  ) / ( sum( f.length for f in feature) )

# Complete permutations of a genome sequence to determine its p-val (posterior probability of alternate hypothesis)
def permutation_p_vals(genome, reference, feature, minimum_overlap, num_permutations, prior_prob_null, prior_prob_alt):
    observed_overlap = sum(f.length for f in (feature.intersect(reference, f=minimum_overlap, u=True))) / (sum(f.length for f in feature))
    permuted_overlap_ratios = np.zeros(num_permutations)

    max_workers = os.cpu_count()
    print("This program uses multithreaded operation. This program is using", max_workers, "threads.")
    
    print("Creating permutations...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(create_permutation, genome, reference, feature, minimum_overlap) for i in range(num_permutations)]
        
        permuted_overlap_ratios = [future.result() for future in concurrent.futures.as_completed(futures)]

    permuted_overlap_ratios = [n*1.0000000000 for n in permuted_overlap_ratios]

    print(permuted_overlap_ratios)
    print(observed_overlap)
    
    permuted_overlap_ratios = np.array(permuted_overlap_ratios)
    
    
    print("Calculating permutation statistics...")
    print(f'{ np.sum( permuted_overlap_ratios >= observed_overlap  ) } of the simulated situations are found to be located at the same location.')
    permutation_dist_given_null = mpf(np.sum( permuted_overlap_ratios >= observed_overlap ) / num_permutations ) 
    permutation_dist_given_alt = 1 - mpf( permutation_dist_given_null ) 

    # heuristic to avoid divide by zero errors
    if permutation_dist_given_null == 0:
        epsilon = 1e-10 
        permutation_dist_given_null  += epsilon
        permutation_dist_given_alt += epsilon
    
    print('null', permutation_dist_given_null)
    print('alt', permutation_dist_given_alt)

    baysean_factor = mpf(permutation_dist_given_alt) / mpf(permutation_dist_given_null)

    prior_odds = mpf(prior_prob_alt) / mpf( prior_prob_null )

    post_odds = mpf( baysean_factor ) * mpf( prior_odds )

    # Return posterior probability of H1
    return mpf( post_odds ) / mpf( post_odds + 1 )


# Calculate the ratio of overlaps between a reference and an array of features of interest
def calculate_overlaps(genome, reference, genomic_features_bedtools, minimum_overlap):
    overlaps = []
    
    print("Calculating overlaps...")

    for feature in genomic_features_bedtools:
        overlaps.append(permutation_p_vals(genome, reference, feature, minimum_overlap, 50, 0.5, 0.5))

    print("Completed calulcation.\n\n\n")
    return overlaps




# Testing with human genome 38 and a few tracks of interest on chromosome 1 
genome = "hg38"
chrom = "chr1"
start = 1 
end = 248956421
genomic_features = ["centromeres" , "cloneEndABC10" , "cloneEndABC11" , "cloneEndABC12" , "cloneEndABC13" , "cloneEndABC14" , "cloneEndABC16" , "cloneEndABC18" , "cloneEndABC20" , "cloneEndABC21" , "cloneEndABC22" , "cloneEndABC23" , "cloneEndABC24" , "cloneEndABC27" , "cloneEndABC7" , "cloneEndABC9" , "cloneEndbadEnds" , "cloneEndCH17" , "cloneEndCOR02" , "cloneEndCOR2A" , "cloneEndCTD" , "cloneEndmultipleMaps" , "cloneEndRP11" , "cloneEndWI2" , "coriellDelDup" , "cpgIslandExt" , "cpgIslandExtUnmasked" , "cytoBand" , "cytoBandIdeo" , "dgvMerged" , "fishClones" , "gap" , "geneReviews" , "genomicSuperDups" , "gold" , "gtexGene" , "gtexGeneV8" , "gwasCatalog" , "hg38ContigDiff" , "hgIkmc" , "iscaBenignGainCum" , "iscaBenignLossCum" , "iscaPathGainCum" , "iscaPathLossCum" , "knownAlt" , "lincRNAsCTAdipose" , "lincRNAsCTAdrenal" , "lincRNAsCTBrain" , "lincRNAsCTBrain_R" , "lincRNAsCTBreast" , "lincRNAsCTColon" , "lincRNAsCTForeskin_R" , "lincRNAsCTHeart" , "lincRNAsCThLF_r1" , "lincRNAsCThLF_r2" , "lincRNAsCTKidney" , "lincRNAsCTLiver" , "lincRNAsCTLung" , "lincRNAsCTLymphNode" , "lincRNAsCTOvary" , "lincRNAsCTPlacenta_R" , "lincRNAsCTProstate" , "lincRNAsCTSkeletalMuscle" , "lincRNAsCTTestes" , "lincRNAsCTTestes_R" , "lincRNAsCTThyroid" , "lincRNAsCTWhiteBloodCell" , "microsat" , "nestedRepeats" , "scaffolds" , "snp141Flagged" , "snp141Mult" , "snp142Flagged" , "snp142Mult" , "snp144Flagged" , "snp144Mult" , "snp146Flagged" , "snp146Mult" , "snp147Flagged" , "snp147Mult" , "snp150Flagged" , "snp150Mult" , "snp151Flagged" , "snp151Mult" , "snpediaText" , "stsMap" , "tRNAs" , "ucscGenePfam" , "ucscToINSDC" , "ucscToRefSeq" , "wgRna"] 
# genomic_features_to_pull = genomic_features[:-1]
base = "track"
pull_new_data = False
min_overlap = 1

features_of_interest_bedtools = create_bedtools(genomic_features, base, genome, chrom, start, end )

reference = pybedtools.example_bedtool(os.path.join(os.getcwd(), base, "hg38_lincRNAsCTColon.bed"))

overlaps = calculate_overlaps(genome, reference, features_of_interest_bedtools, min_overlap)

print(overlaps)
