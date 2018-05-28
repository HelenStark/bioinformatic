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
