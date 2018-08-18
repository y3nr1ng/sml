import logging

import imageio

from . import filters, detectors, estimators, events

logger = logging.getLogger(__name__)

filter_dispatch = {
}

class Analyzer(object):
    """
    Analyze localization events in provided frames.
    """
    def __init__(self, filters='difference'):
        #TODO retrieve processor primitives
        self.events = events.Events()

        self.detector = detectors.LocalMaximum(4, 10)
        
        #TODO distributive computing

    def process_frame(self, frame):
        candidates = self.detector(frame)
        return candidates
