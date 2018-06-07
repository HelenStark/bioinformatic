def number_to_symbol (base):
    switcher = {
        0:"A",
        1:"C",
        2:"G",
        3:"T",
        }
    return switcher.get(base, "fault")

def number_to_pattern (index, k):
    if k==1:
        return number_to_symbol (index)
    prefix_index = index //4
    r= index % 4
    symbol= number_to_symbol(r)
    prefix_pattern = number_to_pattern(prefix_index, k-1)
    return prefix_pattern+symbol

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

def median_string (dna,k):
    distance=k #da distance max.so groß ist wie die Länge der gesuchten Seq. entspricht es unendlich
    for i in range (0, 4**k-1):
        pattern=number_to_pattern(i,k)
        if distance >  distance_pattern_strings(pattern,dna):
            distance=distance_pattern_strings(pattern,dna)
            median=pattern
    return median


