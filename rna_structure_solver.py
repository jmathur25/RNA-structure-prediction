from nussinov import Nussinov
from greedy import GreedySolver
from zuker_freier import ZukerFreier

class RNA_Structure_Solver:
    def __init__(self, mode):
        assert mode in ['greedy', 'nussinov', 'zuker']
        self.mode = mode
        
    def set_mode(self, mode):
        assert mode in ['greedy', 'nussinov', 'zuker']
        self.mode = mode
        
    def solve(self, sequence, min_padding=0, tie_break_permute=None, prettify=False):
        """
            sequence: string of RNA bases represented as A,U,G,C
            min_padding: minimum distance between two paired bases
            tie_break_permute: None or order in which to break ties
            prettify: True if pretty dot paranthesis is to be returned, otherwise False
            
            return: the number of pairings and the dot parantheses notation
        """
        
        if self.mode == 'greedy':
            gr = GreedySolver(sequence)
            score = gr.solve()
            dot_p = gr.dot_parentheses(prettify)
            return score, dot_p
        
        elif self.mode == 'nussinov':
            nv = Nussinov(sequence)
            score = nv.solve(min_padding, tie_break_permute)
            dot_p = nv.dot_parentheses(prettify)
            return score, dot_p
            
        elif self.mode == 'zuker':
            zk = ZukerFreier(sequence)
            score = zk.solve()
            dot_p = zk.dot_parentheses()
            return score, dot_p
        
        
        raise ValueError('mode is not set')
    