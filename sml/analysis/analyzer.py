import logging

import imageio

from . import filters, detectors, estimators

logger = logging.getLogger(__name__)

class Analyzer(object):
    """
    Analyze localization events in provided frames.
    """
    def __init__(self, data):
        pass
