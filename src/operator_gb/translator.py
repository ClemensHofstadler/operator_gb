from __future__ import absolute_import

from .auxiliary import flatten
import re

class Translator:

    def __init__(self,vars):
        if isinstance(vars[0],list):
            X = flatten(vars)
        else:
            X = vars
            
        if len(X) > 128:
            raise ValueError("So far, only <= 128 variables are supported")
        
        self.__d = {str(v) : chr((i + 97) % 128) for i,v in enumerate(X)}
        self.__d_inv = {v:k for k,v in self.__d.items()}
############################################################################                    
    def d(self): return self.__d
    def d_inv(self): return self.__d_inv
    def gens(self): return list(self.__d.keys())
    def internal_gens(self): return list(self.__d_inv.keys())
############################################################################                    
    def __repr__(self):
        s = "Translator "
        for k,v in self.d().items():
            s += k + "<->" + v + ", "
        return s[:-2]
    
############################################################################                        
    def __call__(self, string, to_internal=True):
        rep_dict = self.__d if to_internal else self.__d_inv
        pattern = re.compile("|".join([re.escape(k) for k in sorted(rep_dict,key=len,reverse=True)]), flags=re.DOTALL)
        return pattern.sub(lambda x: rep_dict[x.group(0)], string)
            


    
