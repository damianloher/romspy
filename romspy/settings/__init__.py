from .settings import Settings
from .forcing import forcing_adjustments
from .clim import clim_adjustments

forcing_settings = Settings(forcing_adjustments, None)
clim_settings = Settings(clim_adjustments, None)
