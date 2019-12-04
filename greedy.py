import queue

class GreedySolver:
    def __init__(self, s):
        """
        s: sequence of RNA
        """
        
        self.s = s
        # pairing dictionary
        self.T = {
            'A': 'U',
            'G': 'C',
            'C': 'G',
            'U': 'A',
        }

    def solve(self):
        """    
        Very simple algorithm. It reads the string left to right, and for each read it sees if it has
        seen a match for that letter earlier (which is maintained in a queue for that letter that matches). If the 
        queue is empty, no match is possible, so add the current letter to it's queue
        
        Conceptually, the max base pairs for any string is min(count(A), count(U)) + min(count(C), count(G))
        This algorithm finds the pairings greedily, without care for pseudoknots
        """
        letter_to_queue = {k: queue.Queue() for k in self.T}
        
        score = 0
        self.pairings = []
        
        for i, l in enumerate(self.s):
            match = self.T[l] # grab complement
            if not letter_to_queue[match].empty():
                # match current letter by popping it from the queue for the match letter
                idx = letter_to_queue[match].get()
                self.pairings.append((i, idx))
                score += 1
            else:
                # store index in the queue
                letter_to_queue[l].put(i)
        
        return score

    def dot_parentheses(self, prettify=True):
        result = ['-' for _ in range(len(self.s))] # gapped sequence
        for pair_idx1, pair_idx2 in self.pairings:
            if pair_idx1 > pair_idx2:
                pair_idx1, pair_idx2 = pair_idx2, pair_idx1
            result[pair_idx1] = '('
            result[pair_idx2] = ')'
            
        # the raw string
        result = ' '.join(result)
        if not prettify:
            return result
        
        # printing
        ref = ' '.join(self.s)
        middle = ''
        for _ in range(len(ref)):
            middle += '-'
        return ref + '\n' + middle + '\n' + result

    