def count_p(s):
    ret = 0
    for c in s:
        if c == '(':
            ret += 1
    return ret

def find_overlap(s1, s2):
    assert len(s1) == len(s2)
    
    shared = 0
    for i in range(len(s1)):
        if s1[i] == '(' and s2[i] == '(':
            shared += 1
        elif s1[i] == ')' and s2[i] == ')':
            shared += 1
    return shared/2
