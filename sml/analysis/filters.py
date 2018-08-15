from abc import ABCMeta, abstractmethod

import numpy as np

__all__ = [
    'DifferenceImaging',
    'WaveletFilter'
]

class Base(object):
    @abstractmethod
    def __call__(self, curr_frame):
        pass

class WaveletFilter(Base):
    pass
