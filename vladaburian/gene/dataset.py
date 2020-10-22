

class Dataset2D:
    def __init__(self, rows, cols):
        self._rows = rows
        self._cols = cols
        self._data = [[{} for _ in range(self._cols)] for _ in range(self._rows)]

        self._rowsDict = {k:i for i,k in enumerate(self._rows)}
        self._colsDict = {k:i for i,k in enumerate(self._cols)}

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            raise RuntimeError()

        if len(key) == 2:


    def __setitem__(self, key, value):
        pass


"""
a = A()

a[1]
a[1:5]
a[:]
a['a']
a[1,5]
a[:5,4]
a[:5,:]
a[52,:]
a[52,slice(None,None,None)]
"""
