import logging

from ._src import get_centerline


__version__ = "0.2"

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
