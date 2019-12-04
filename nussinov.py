import queue
import itertools
import tqdm

from utils import AlteredList

class Node:
    # stored in the DP table while filling it out for easy traceback
    score = 0
    coord = None
    # two in case of bifurcation
    back_pointer1 = None
    back_pointer2 = None
    # useful booleans for backtrace and other functions
    bifurcate = False
    match = False

    def __repr__(self):
        return str(self.score)


class Nussinov:
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

        
    def solve(self, min_padding=0, tie_break_permute=None):
        """
        DP solution to nussinov. Runtime: |s|^3
        
        min_padding: padding between i and j at minimum
        tie_break_permute: pass a list that tells what order to prioritize the following if there's a tie
            - left (L)
            - down (D)
            - match (M)
            - bifurcate (B)
        """        
        
        # initializes variables for the class to use
        self.table = AlteredList(len(self.s), fill_value=Node)
        self.trace = None
        
        # we fill the table diagonally, I'll use a while loop cause it helps me
        cur_i, cur_j = 0, 1
        next_i_start, next_j_start = 0, 2
        
        # if over 50 characters, I'll show a progress bar
        need_bar = False if len(self.s) < 50 else True
        # progress bar, total found via some simple geometry on the table
        if need_bar:
            pbar = tqdm.tqdm(total= (len(self.s)-1)*len(self.s)/2)
        
        prev, count = 0, 0
        while cur_i != 0 or cur_j != len(self.s):

            bifurcation = None
            if cur_j - cur_i >= 2: # if bifurcation is possible
                options = [(self.table[cur_i,k].score + self.table[k+1,cur_j].score, (cur_i,k), (k+1,cur_j)) for k in range(cur_i+1, cur_j)]
                bifurcation = max(options, key = lambda x: x[0]) # choose max based on score (index 0 of tuple)
            else:
                bifurcation = (0, None, None)
            
            match = None
            # if we have a match in the allowable distance
            if self.T[self.s[cur_i]] == self.s[cur_j] and cur_j - cur_i - 1 >= min_padding:
                match = (1 + self.table[cur_i+1,cur_j-1].score, (cur_i+1, cur_j-1))
            else:
                match = (0, None)

            down = (self.table[cur_i+1,cur_j].score, (cur_i+1,cur_j))
            left = (self.table[cur_i,cur_j-1].score, (cur_i,cur_j-1))

            tie_break_order = None
            # tie breaking happens here
            if tie_break_permute is not None:
                tie_map = {'L': left, 'D': down, 'M': match, 'B': bifurcation}
                tie_break_order = [tie_map[l] for l in tie_break_permute]
                
            else:
                # default
                tie_break_order = [left, down, match, bifurcation]
                
            choice = max(tie_break_order, key=lambda x: x[0])
            new_node = Node()
            new_node.score = choice[0]
            new_node.coord = (cur_i, cur_j)
            if choice[1] == (cur_i+1,cur_j-1): # sees if we chose to match
                new_node.match = True
            
            if len(choice) == 3: # bifurcation
                new_node.back_pointer1 = choice[1]
                new_node.back_pointer2 = choice[2]
                new_node.bifurcate = True
            else: # regular
                new_node.back_pointer1 = choice[1]
            
            self.table[cur_i, cur_j] = new_node

            cur_i += 1
            cur_j += 1
            
            if need_bar:
                count += 1
                if count % 100 == 0:
                    pbar.update(count - prev)
                    prev = count

            if cur_j >= len(self.s):
                cur_i = next_i_start
                cur_j = next_j_start
                next_j_start += 1
                
        if need_bar:
            pbar.update(count - prev)    
            pbar.close()

        return self.table[0, len(self.s)-1].score # top right corner has the answer
    
    
    def backtrace(self):
        # backtraces and returns sequences of nodes that led to the final path
        if self.trace is not None: # already computed
            return self.trace
        
        start_node = self.table[0,len(self.s)-1]
        self.trace = []
        
        node_queue = queue.Queue()
        node_queue.put(start_node)
        
        while not node_queue.empty():
            cur_node = node_queue.get()
            self.trace.append(cur_node)
            
            if cur_node.back_pointer1 is not None:
                pointer = cur_node.back_pointer1
                node_queue.put(self.table[pointer])
                
            if cur_node.back_pointer2 is not None:
                pointer = cur_node.back_pointer2
                node_queue.put(self.table[pointer])
                
        return self.trace
    
    
    def dot_parentheses(self, prettify=True):
        """
        Creates a string showing the dot parenthesis version of the structure
        """
        if self.trace is None:
            self.backtrace()
        
        result = ['-' for _ in range(len(self.s))] # gapped sequence
        for node in self.trace:
            if not node.match:
                continue
            i, j = node.coord
            result[i] = '('
            result[j] = ')'
            
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
        
        
    def evaluate_tie_breaks(self, min_padding=0, prettify=True):
        """
        Passes all tie permutations to save and stores the dot_parenthesis computation
        """
        to_permute = ['L', 'D', 'M', 'B']
        perms = list(itertools.permutations(to_permute))
        
        results = {}
        for perm in perms:
            self.solve(min_padding, tie_break_permute=perm)
            dot_parenth = self.dot_parentheses(prettify)
            if dot_parenth in results:
                results[dot_parenth].append(perm)
            else:
                results[dot_parenth] = [perm]
        return results
        
