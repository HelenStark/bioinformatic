def hamming_distance(pattern, k_mer):
    #anzahl der falschen Base in einem k-mer im Vergleich zu pattern
    hamming_distance = 0
    i=0
    while i < len(pattern):
        if pattern[i] != k_mer[i]:
            hamming_distance = hamming_distance + 1
        i = i+1
    return hamming_distance

def distance_pattern_strings (pattern, dna):
    k=len(pattern)  
    distance=0
    for text in dna: #für jeden dna strang 
        ha_distance=k+1 #max sind alle Basen verschieden,
                    #daher k+1 immer höher als das mögliche
        for i in range (0, len(text)-k+1):
            k_mer = text [i:i+k] #alle kmer in text 
            hamming_dist = hamming_distance(pattern, k_mer) #Variabelname nicht gleich Fkt.
            if ha_distance > hamming_dist: 
                ha_distance=hamming_dist
        distance = distance + ha_distance
    return distance
