import logging

import imageio

from . import filters, detectors, estimators, events

logger = logging.getLogger(__name__)

class Analyzer(object):
    """
    Analyze localization events in provided frames.
    """
    def __init__(self, filters='difference'):
        #TODO retrieve processor primitives
        #self.events = events.Events(2)

        self.detector = detectors.DifferenceImaging(4, 10)
        #self.estimator = estimators.GaussianMLE()

        #TODO distributive computing

    def process_frame(self, frame):
        candidates = self.detector(frame)
        coords = self.estimator(frame, candidates)
        return coords
