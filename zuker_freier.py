import tqdm

from zuker_freier_utils import FreierPair, Q_const, FreierDangleLeft, FreierDangleRight, M_const, FreierPair, Qi_const, FreierDangleLeft, FreierDangleRight, FreierIS1, FreierIS2, m_const

class ZukerFreier:
    def __init__(self, s):
        """
        s: sequence of RNA
        """
        
        self.s = s
        self._init_matrices()
        
    def solve(self, min_padding=3):
        global m_const
        m_const = min_padding

        # we fill the table diagonally, I'll use a while loop cause it makes life easier
        cur_i, cur_j = 0, 1
        next_i_start, next_j_start = 0, 2
        
        # if over 50 characters, I'll show a progress bar
        need_bar = False if len(self.s) < 50 else True
        # progress bar, total found via some simple geometry on the table
        if need_bar:
            pbar = tqdm.tqdm(total= (len(self.s)-1)*len(self.s)/2)

        prev, count = 0, 0
        while cur_i != 0 or cur_j != len(self.s):
            self._wx(cur_i, cur_j)

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

        return self._wx(0, len(self.s) - 1)
        
        
    def _wx(self, i, j):
        if (j - i < m_const) :
            return 0
    
        if self.W[i][j] != "-":
            return self.W[i][j]
            
        min_k = self._get_min_k(i, j)
        possibilities = [FreierPair(self.s, i,j) + self._vx(i,j), FreierDangleLeft(self.s, i, i+1, j-1) + FreierDangleRight(self.s, j, i+1, j-1) + FreierPair(self.s, i + 1, j - 1) + self._vx(i+1, j-1),
                         FreierDangleLeft(self.s, i, i+1, j) + FreierPair(self.s, i + 1, j) + self._vx(i+1, j), FreierDangleRight(self.s, j, i, j-1) + FreierPair(self.s, i, j - 1) + self._vx(i, j-1),
                        Q_const + self._wx(i + 1, j), Q_const + self._wx(i, j - 1), min_k[0]]
        
        self.W[i][j] = min(possibilities)
        self.TB[i][j] = possibilities.index(self.W[i][j])
        if self.TB[i][j]  == 6:
            self.TB[i][j] = min_k[1] * -1
        self.W[i][j] = round(self.W[i][j], 1)
        return self.W[i][j]
    
    def _vx(self, i,j):
        if (FreierPair(self.s, i,j) == float('inf')):
            return float('inf')
        if(j - i < m_const) :
            return float('inf')
        if self.V[i][j] != "-":
            return self.V[i][j]
        possibilities = [FreierIS1(i,j), self._optimal_IS2(i,j), self._optimal_IS3(i,j)]

        self.V[i][j] = round(min(possibilities), 2)
        return self.V[i][j]
    
    def _optimal_IS3(self, i, j) :
        minimum = float('inf')
        for k in range (i + 1, j - 1) :
            if k - i - 1 >= m_const and j - k - 1 >= m_const:
                multi = M_const + self._wxi(i + 1,k) + self._wxi(k + 1,j - 1)
                if minimum > multi :
                    minimum = multi
        return minimum
    
    def _wxi(self, i, j):
        if(j - i < m_const) :
            return 0
        if self.Wi[i][j] != "-":
            return self.Wi[i][j]
        min_k = self._get_min_k(i, j)
        possibilities = [FreierPair(self.s, i,j) + self._vx(i,j), FreierDangleLeft(self.s, i, i+1, j-1) + FreierDangleRight(self.s, j, i+1, j-1) + FreierPair(self.s, i + 1,j - 1) + self._vx(i+1,j-1),
                         FreierDangleLeft(self.s, i, i+1, j) + FreierPair(self.s, i + 1,j) + self._vx(i+1,j), FreierDangleRight(self.s, j, i, j-1) + FreierPair(self.s, i,j - 1) + self._vx(i,j-1),
                        Qi_const + self._wxi(i + 1,j), Qi_const + self._wxi(i,j - 1), min_k[0]]
        self.Wi[i][j] = min(possibilities)
        self.TB[i][j] = possibilities.index(self.Wi[i][j])
        if(self.TB[i][j]) == 6:
            self.TB[i][j] = min_k[1] * -1
        self.Wi[i][j] = round(self.Wi[i][j], 1)
        return self.Wi[i][j]

    def _optimal_IS2(self, i, j):
        if(FreierPair(self.s, i,j) == float('inf')):
            return float('inf')
        minimum = float('inf')
        for k in range(i + 1, j - 1):
            for l in range (k + 1, j):
                is2 = FreierIS2(self.s, i, j, k, l) + FreierPair(self.s, k,l) + self._vx(k, l)
                if(is2 < minimum):
                    minimum = is2
        return minimum
    
    def _get_min_k(self, i, j):
        minimum = float('inf')
        min_k = -1
        for k in range (i + 1, j - 1) :
            result = self._wx(i,k) + self._wx(k + 1,j)
            if minimum > result :
                minimum = result
                min_k = k
        return minimum, min_k
    
    def _init_matrices(self):
        self.V = []
        self.W = []
        self.Wi = []
        self.TB = []
        
        for x in range(len(self.s)):
            self.V.append([])
            self.W.append([])
            self.Wi.append([])
            self.TB.append([])
            for _ in range(len(self.s)):
                self.V[x].append('-')
                self.W[x].append('-')
                self.Wi[x].append('-')
                self.TB[x].append('-')
        for i in reversed(range(len(self.V))):
            for j in range(len(self.V[i])):
                if j - i <= m_const:
                    self.V[i][j] = float('inf')
                    self.W[i][j] = 0
                    self.Wi[i][j] = 0
                    self.TB[i][j] = 6
        
    def dot_parentheses(self):
        structure = []
        for i in range(len(self.s)):
            structure.append('.')
        i = 0
        j = len(self.s) - 1
        while i < j :
            if self.TB[i][j] == 0:
                structure[i] = '('
                structure[j] = ')'
                i+= 1
                j-= 1
            elif self.TB[i][j] == 1:
                i+= 1
                j-= 1   
            elif self.TB[i][j] == 2 or 4:
                i+= 1
            elif self.TB[i][j] == 3 or 5:
                j-= 1

        string = ''
        for i in range(len(structure)):
            string +=(structure[i])
        return string
        
    