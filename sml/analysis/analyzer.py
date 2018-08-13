import logging

import imageio

from . import filters, detectors, estimators, events

logger = logging.getLogger(__name__)

filter_dispatch = {
    'difference': filters.DifferenceImaging
}

class Analyzer(object):
    """
    Analyze localization events in provided frames.
    """
    def __init__(self, data, filters='difference'):
        self.data = data
        self.filters = filter_dispatch.get(filters)
        self.detectors = None
        self.estimators = None
        self.events = events.Event()
