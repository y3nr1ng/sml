import logging
import imageio

from . import filters, detectors, estimators, events

logger = logging.getLogger(__name__)

filter_dispatch = {
    'difference': detectors.DifferenceImaging,
}
estimator_dispatch = {
    'GaussianMLE': estimators.GaussianMLE,
    'Phasor': estimators.Phasor
}

class Analyzer(object):
    """
    Analyze localization events in provided frames.
    """
    def __init__(self, settings, filters='difference'):
        #TODO retrieve processor primitives
        #self.events = events.Events(2)

        self.detector = filter_dispatch.get(filters)(
                            frameRange = (settings['FRAME_ID_START'],
                                          settings['FRAME_ID_END']),
                            ROI = ((settings['FRAME_X1'],settings['FRAME_Y1']),
                                   (settings['FRAME_X2'],settings['FRAME_Y2'])),
                            radius = (settings['X_FIND_PIXELS'],
                                      settings['Y_FIND_PIXELS']),
                            thresMax = settings['I_THRES_MAX'],
                            thresMin = settings['I_THRES_MIN'])
#        self.estimator = estimator_dispatch.get('GaussianMLE')

        #TODO distributive computing

    def process(self, frame):
        candidates = self.detector(frame)
        print("# of candidates: %d" % len(candidates))
#        coords = self.estimator(frame, candidates)
#        return coords
