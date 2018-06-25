def number_to_symbol (base):
    switcher = {
        0:"A",
        1:"C",
        2:"G",
        3:"T",
        }
    return switcher.get(base, "fault")

def hamming_distance(motif1, motif2):
    hamming_distance = 0
    i=0
    while i < len(motif1):
        if motif1[i] != motif2[i]:
            hamming_distance = hamming_distance + 1
        i = i+1
    return hamming_distance

def get_consensus(motifs):#aus Profil wird Konsensus seq. generiert
    profile=get_profile(motifs)
    consensus = ''
    for i in range(0, len(motifs[0])): #für jede Position
        liste2 = []
        for lis in profile.values(): #der Eintrag aus dem Profil(dictonary aus Listen) mit dem passenden Index wird in liste 2 gespeichert
            liste2.append(lis[i])
        
        ind = liste2.index(max(liste2)) #der index des maximums repräsentiert die Base an dieser Stelle
        base = number_to_symbol(ind) #der index wird in die base übersetzt
        consensus = consensus + base #der consensus string wird erweitert um Base
    return consensus

def get_profile(motifs):
    k=len(motifs[0]) #Vorraussetzung alle Motive gleich lang
    profile = {
        'A': [1/(len(motifs)+4)]*k, #pseudocounts -> starten bei 1 nicht bei 0
        'C': [1/(len(motifs)+4)]*k, #teilen um direkt auf profil zu kommen und nicht erst auf Count
        'G': [1/(len(motifs)+4)]*k, 
        'T': [1/(len(motifs)+4)]*k,
        }
    for motif in motifs: #für jedes Motiv in der liste motifs
        for index, base in enumerate(motif): #sowohl index als auch base werden für jede Position verwendet
            profile[base][index] += 1/(len(motifs)+4) #dem profil wird in der Liste mit der entsprechenden Base am Index entsprechender Ausdruck hinzugefügt
    return profile 

def get_best_motifs(profile,dna): #mithilfe des Profils wird jeweils das beste motif aus der jeweiligen seq. gesucht
    k=len(profile["A"])
    best_probs=[] #Liste der besten Wahrscheinlichkeit pro Sequenz
    motifs=[] #Liste der zugehörigen besten k_mere zu gegebenem Profil
    for sequence in dna:
        bestprob=0
        for i in range(0,len(sequence)-k+1):
            k_mer=sequence [i:i+k] #Definition was ist kmer
            prob=1
            index=0 #man könnte auch wie in Fkt. profile über enumerate index und base bekommen.
            for base in k_mer:
                p=profile[base][index] #Wahrscheinlichkeit für Base an Stelle index
                prob=prob*p #Multiplikation der einzelnen Wahrscheinlichkeiten für k-mer
                index+=1
            if prob>bestprob: #nur wenn Wahrscheinlichkeit k-mer höher/besser als vorher
                bestprob=prob
                best_k_mer=k_mer
        best_probs.append(bestprob) #Hinzufügen zu Listen und Ausgabe der Fkt.
        motifs.append(best_k_mer)
    return motifs #aus interesse könnte man sich auch die Wahrscheinlichkeiten best_probs ausgeben, ist für die Aufgabe aber nicht nötig

def score(motifs):
    consensus = get_consensus(motifs)
    score = 0
    for motif in motifs: #für jedes motif wird die hamming_distance zur consensussequenz berechnet
        hdist = hamming_distance(consensus, motif)
        score = score + hdist
    return score

def randomized_motif_search(dna,k,t):
    motifs=[]
    import random
    for sequence in dna:
        i=random.randint(0,len(sequence)-k) #zufällig wird die Startposition i gewählt
        k_mer=sequence[i:i+k] #und ein k-mer ausgesucht
        motifs.append(k_mer)
    best_motifs=motifs
    while True: #while schleife läuft solange, bis sich der score bzw. die Ausgabe nicht mehr verbessern.
        profile=get_profile(motifs)
        motifs=get_best_motifs(profile,dna)
        if score(motifs) < score(best_motifs): #wenn der score der motifs (werden aus profil generiert) kleiner ist als der score der best_motifs (in erstem durchlauf die random kmere)
            best_motifs=motifs #die neuen motifs werden als best_motifs genommen und ein neues profil generiert
        else:
            return best_motifs #zum schluss ausgabe der motive

dna=['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG','TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
dna2=["ACTTATATCTAGAGTAAAGCCCTGATTCCATTGACGCGATCCCTACCTCCATCATACTCCACAGGTTCTTCAATAGAACATGGGGAAAACTGAGGTACACCAGGTCTAACGGAGATTTCTGGCACTAACTACCCAAAATCGAGTGATTGAACTGACTTATATCTAGAGT",
    "AAAGCCCTGATTCCATTGACGCGATCCCTACCTCCATCATACTCCACAGGTTCTTCAATAGAACATGGGGAAAACTGAGGTACACCAGGTCTAACGGAGATTTCTGGCACTAACTACCCAAAATCCTCTCGATCACCGACGAGTGATTGAACTGACTTATATCTAGAGT",
    "CACTCCCGTCCGTCTGACGCCAGGTGCTCTACCCCGCTGATTGTCTGGTACATAGCAGCCTATAGATCACCGATGCAGAAACACTTCGAGGCAGCCGATTTCGCTTATCACAACGTGACGGAATTTGATAAACCACGTACTCTAATACCGTCACGGGCCCATCAACGAA",
    "ACAAGAACTGGTGGGGAGACTATGACACTCTAGCGGTCGCATAAGGGCCGGAAACCAGGACAAATCGATAAGATGAAGCGGGGATATAAGCCTTATACTGCGACTGGTTCCTTATATTATTTAGCCCCGATTGATCACCGATTAAAATATTCTGCGGTTTTCGAGACGG",
    "TAACCACACCTAAAATTTTTCTTGGTGAGATGGACCCCCGCCGTAAATATCAGGATTAAATGTACGGATACCCATGACCCTCCAGTCATCTACCTTCCCGTGGTGGTCGCTCAGCCTTGTGCAGACCGAACTAGCACCTGTCACATACAATGTTGCCCGCATAGATCGT",
    "ATCCGACAGAGGCAGTGAATAAGGTTTCGTTTCCTCAGAGAGTAGAACTGCGTGTGACCTTGCCTTCACCGACATCCGTTTCCAATTGAGCTTTTCAGGACGTTTAGGTAACTGATTGTCATTGCAATTGTCCGGGGGATTTAGATGGCCGGGTACCTCTCGGACTATA",
    "CCTTGTTGCCACCGATTCGCGAGCAACATCGGAGTGCTCTGATTCACGGCGATGCTCCACGAAGAGGACCGCGGCACGACACGCCCTGTACCTACGTTTCTGGATATCCTCCGGCGAGTTAATAGAGCAATACGACCTGGTCGTCGAGATCGTGTATCTAGCCCTACCT",
    "ATAGGTTAACGAATCAGGAGAGTTAATTTTACCTAGCTAGAGCGGACGGTGCCTGGCTGTATTCGCGTTTGACTTTCGGGCTCGCTGATAACTTGTGATCACCTTTTACGCTTACTGGATCCAACGATGGATCAAAGTTGAGAATTTCTGTGCCTTGGGTGTGAGCTGT",
    "CTGACGAAAGGACGGGCGGTGTACTTAGTTTGGGGTAAAATAGTTGGTATAATTCTGTGCGACAGACATTTGGTCAGGCCATACTGCCATATCGTGATGTAACTATCCACACTACGTCATAGGCCCTTGTGATCAATTAAACGTTCCTCATGCCAGGCTATCTGTTTAA",
    "GGCTTCGCGTTTAAGGCTGGATTAAGTACTCCGCCTTGTGATCTGTGATCCTCCGACCTGTGATCAGCAAGATTGGAACCTAGGTAGGCGGCGGGTCTACGCTGGCCCACAATCGTGAGTCCCCCACTCCGTAGGTTGTGGAATTTATAGACCCGCAAGGGGCACCACT",
    "AGGATGACACCCAGGATGAATCTGGATTAGGAACACCAACCCGACATATTTGTTACCGCTGCAGCATTTCGCTCTTGGACGCGTAACCCGAGATCCGTCTCGCGATCGTCACGGATCGGGATTATGCAGGCAATACCTTGTGATCACTCCGCGCTTGGTTTTGCTAGCG",
    "ACATCTCTAGTCACTTTTATTGAGCAGGTGGGCGGATTCATGATCCGGCTCTGTCGTACGTCCAACCACGGTGACATGTTCGGAGCTGTCGCCGTGGAGCAGAGATACATCGGATCTATCAATTTTACTAAGAGCAACTAGCCACGACAAACTGTGATCACCGATTGGA",
    "AATTTGCGTATCTCTAGGACTCCCTCATACAAATCAAAGCTTGGATGGGTAAGATGCCGCAGCAGCAGGTATCTCATATTGGCTATTAAGAGCCAGGCCCTATGGCCTTAGTATCACCGATCAGACGTCGCATGAGCGGGCCCGTTGTCCTATCTCTTTAGCTGCCGCA",
    "GAAGTAAAGGGGTTCCACTGCGTAGAGCGTGCCCCTCTGGTGTGCCGTACTGTTATGGTGATACAGCTTCCTTATACCCCTCGTAAAGCGGCTAATGGTCCTAATGAATGCCCTTGTGAAATCCGAATCGCTTTACAATTGCGTTCGGCGGAATGCAGTCACCAGTGTT",
    "TACACTACGCGTTATTTACTTTTACTGAGTCCTTGTCGCCACCGAACGAGGATTGTTCATTGTATCCGGAGATTAGGAGTTCGCATCGCTGACACAGCCAGTTCGTAGCAAATACCGCTGGCCCTGGGCACTCCAGATCAGAACTACTAGCCCTAAACTCTATGACACA",
    "TTGGGTCTCGATCCCTCTATGTTAAGCTGTTCCGTGGAGAATCTCCTGGGTTTTATGATTTGAATGACGAGAATTGGGAAGTCGGGATGTTGTGATCACCGCCGTTCGCTTTCATAAATGAACCCCTTTTTTTCAGCAGACGGTGGCCTTTCCCTTTCATCATTATACA",
    "TTTCAAGTTACTACCGCCCTCTAGCGATAGAACTGAGGCAAATCATACACCGTGATCACCGACCCATGGAGTTTGACTCAGATTTACACTTTTAGGGGAACATGTTTGTCGGTCAGAGGTGTCAATTATTAGCAGATATCCCCCAACGCAGCGAGAGAGCACGGAGTGA",
    "GATCCATTACCCTACGATATGTATATAGCGCCCTAGTACGGCTTCTCCCTTGCAGACACGCAGGCGCTGTGCGCTATCGGCTTCCTCGGACATTCCTGGATATAAGTAACGGCGAACTGGCTATCACTACCGCCGCTCCTTAAGCCTTGGTTTCACCGACGATTGTCGT", 
    "TAGTAGATTATTACCTGTGGACCGTTAGCTTCAAGACCGAAACGTTGGTGATGCTACTTAAATGTCAAGAGTTGCGAAGTTGGGCGAAGCACATCCGTACTCCCAAGTGGACGATCGATAGATCCATGGAGTTTCCATCCATCTTAATCCGCCCTTTGCATCACCGACG",
    "TACAAGGCACAAACGAGACCTGATCGAACGGTGCACGGTCGAGGCAGCGAGATAAATGTACATTGAGAGCACCTTGTGATTTACGACCTGCATCGAAGGTTTCTTGGCACCCACCTGTCGTCCGCCAGGGCAGAGCCGACATTATATGACGCTGATGTACGAAGCCCCT"]
scores=[]
motifs=[]
for x in range(0,1000): #randomized_motif_search wird in diesem Falle 1000 mal angesprochen
    result=randomized_motif_search(dna,8,5)
    sc=score(result) #der score des ergebnisses aus randomized_motif_search wird berechnet und gespeichert
    scores.append(sc)
    motifs.append(result) #die motive des ergebnisses werden in Liste motifs gespeichert

index=scores.index(min(scores)) #der index aus dem minimalen score (entspricht bester score) wird gespeichert
print(motifs[index],scores[index]) #die motive und der score mit dem besten Ergebnis werden ausgegeben.
