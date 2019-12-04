# thin wrapper over typical Python list functionality to allow NumPy style indexing
class AlteredList:
    def __init__(self, length, fill_value=0):
        if callable(fill_value):
            self.list = [[fill_value() for _ in range(length)] for _ in range(length)]
        else:
            self.list = [[0 for _ in range(length)] for _ in range(length)]
    
    def __getitem__(self, i):
        if type(i) == int:
            return self.list[i]
        elif type(i) == list or type(i) == tuple:
            assert len(i) == 2
            return self.list[i[0]][i[1]]
            
        raise ValueError("cannot access list with an index of type {}".format(type(i)))
        
    
    def __setitem__(self, i, obj):
        if type(i) == int:
            self.list[i] = obj
            return
        
        elif type(i) == list or type(i) == tuple:
            assert len(i) == 2
            self.list[i[0]][i[1]] = obj
            return
            
        raise ValueError("cannot access list with an index of type {}".format(type(i)))
        
        
    def __len__(self):
        return len(self.list)
