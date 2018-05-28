def symbol_to_number (base):
    switcher = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3,
    }
    #dictonary, welches die Basen ersetzt
    return switcher.get(base, "fault")

def pattern_to_number (pattern):
    if pattern == "":
        return 0
    symbol = pattern [-1]
    prefix = pattern [0:len(pattern)-1]
    return 4*pattern_to_number (prefix) + symbol_to_number (symbol)

def computing_frequencies(text,k):
    frequency_array = [0]*(4**k)
    for i in range (0,len(text)-k+1):
        pattern = text [i:i+k]
        j = pattern_to_number (pattern)
        frequency_array [j] = frequency_array [j] +1
    return frequency_array

