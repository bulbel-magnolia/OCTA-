"""Auditable rectangular-tail quantification for OCTA speckle variance."""

from .geometry import VesselGeometry
from .quantification import QuantificationResult, quantify_frame

__all__ = ["QuantificationResult", "VesselGeometry", "quantify_frame"]
__version__ = "0.1.0"
