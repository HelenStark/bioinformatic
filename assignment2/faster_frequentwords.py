def symbol_to_number (base):
    switcher = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3,
    }
    #dictonary, welches die Basen ersetzt
    return switcher.get(base, "fault")

def number_to_symbol (base):
    switcher = {
        0:"A",
        1:"C",
        2:"G",
        3:"T",
        }
    return switcher.get(base, "fault")

def pattern_to_number (pattern):
    if pattern == "":
        return 0
    symbol = pattern [-1]
    prefix = pattern [0:len(pattern)-1]
    return 4*pattern_to_number (prefix) + symbol_to_number (symbol)
    
def number_to_pattern (index, k):
    if k==1:
        return number_to_symbol (index)
    prefix_index = index //4
    r= index % 4
    symbol= number_to_symbol(r)
    prefix_pattern = number_to_pattern(prefix_index, k-1)
    return prefix_pattern+symbol

def computing_frequencies(text,k):
    frequency_array = [0]*(4**k)
    for i in range (0,len(text)-k+1):
        pattern = text [i:i+k]
        j = pattern_to_number (pattern)
        frequency_array [j] = frequency_array [j] +1
    return frequency_array

def faster_frequentwords (text, k):
    frequent_patterns = []
    frequency_array = computing_frequencies (text,k)
    maxcount= max(frequency_array)
    for i in range (0, 4**k-1):
        if frequency_array [i] == maxcount:
            pattern= number_to_pattern (i,k)
            frequent_patterns.append(pattern)
    return list(set(frequent_patterns))

