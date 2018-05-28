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
