from ctypes import *
from os import *

class blasfeo_dmat(Structure):
    _fields_ = [    
        ("m", c_int),                  # rows
	("n", c_int),                  # cols
	("pm", c_int),                 # packed number or rows
	("cn", c_int),                 # packed number or cols
	("pA", POINTER(c_double)),     # pointer to a pm*pn array of doubles, the first is aligned to cache line size
	("dA", POINTER(c_double)),     # pointer to a min(m,n) (or max???) array of doubles
	("use_dA", c_int),             # flag to tell if dA can be used
	("memsize", c_int)             # size of needed memory
        ]




class blasfeo_dvec(Structure):
    _fields_ = [    
        ("m", c_int),                   # size
	("pm", c_int),                  # packed size
        ("pa", POINTER(c_double)),      # pointer to a pm array of doubles, the first is aligned to cache line size
	("memsize", c_int)              # size of needed memory
        ]
