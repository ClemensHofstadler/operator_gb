from __future__ import absolute_import

from .auxiliary import flatten
import re

def _internal_alphabet():
    r"""
    The single characters used to represent variables internally.

    Monomials are (de)serialised as strings using ``*`` as separator and ``^``
    for powers, so neither of these two characters, nor any digit, may be used
    as an internal variable. Starts with ``a``, ``b``, ``c``, ... to keep the
    historical layout for small numbers of variables.
    """
    forbidden = set('*^')
    letters = [chr(i) for i in range(97,123)] + [chr(i) for i in range(65,91)]
    punctuation = [chr(i) for i in range(33,127)
                   if not chr(i).isalnum() and chr(i) not in forbidden]
    latin1 = [chr(i) for i in range(161,256) if not chr(i).isdigit()]
    return letters + punctuation + latin1

INTERNAL_ALPHABET = _internal_alphabet()

class Translator:

    def __init__(self,vars):
        if isinstance(vars[0],list):
            X = flatten(vars)
        else:
            X = vars
            
        if len(X) > len(INTERNAL_ALPHABET):
            raise ValueError("So far, only <= %d variables are supported"
                             % len(INTERNAL_ALPHABET))
        
        self.__d = {str(v) : INTERNAL_ALPHABET[i] for i,v in enumerate(X)}
        self.__d_inv = {v:k for k,v in self.__d.items()}
        self.__pattern = self.__compile(self.__d)
        self.__pattern_inv = self.__compile(self.__d_inv)
############################################################################                    
    @staticmethod
    def __compile(rep_dict):
        keys = sorted(rep_dict,key=len,reverse=True)
        return re.compile("|".join([re.escape(k) for k in keys]),flags=re.DOTALL)
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
        pattern = self.__pattern if to_internal else self.__pattern_inv
        return pattern.sub(lambda x: rep_dict[x.group(0)], string)
            


    
