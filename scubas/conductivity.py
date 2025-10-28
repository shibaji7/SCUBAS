#!/usr/bin/env python

"""
conductivity.py: Earth conductivity profile utilities backed by the LITHO1.0 model.
"""

from __future__ import annotations

__author__ = "Murphy, B.; Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import shutil
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Union

import numpy as np
import pandas as pd
from loguru import logger
from pyproj import Geod
from scipy.interpolate import RegularGridInterpolator
from scipy.io import netcdf_file

from scubas.datasets import Site
from scubas.utils import RecursiveNamespace


class ConductivityProfile:
    """
    Build conductivity (resistivity) profiles from the LITHO1.0 reference model.

    The class exposes convenience helpers to retrieve conductivity profiles for
    points, bins, and ensembles of geographical coordinates, optionally
    returning SCUBAS ``Site`` instances ready for downstream modelling.
    """

    def __init__(
        self,
        conductivity_params: Optional[Mapping[str, Any]] = None,
    ) -> None:
        """
        Parameters
        ----------
        conductivity_params :
            Optional mapping that overrides the default model configuration
            (e.g. file name, interpolation method, resistivity values). See
            :mod:`scubas.conductivity` source for supported entries.
        """
        default_conductivity_params: Dict[str, Any] = {
            "earth_model": "LITHO1.0.nc",
            "grid_interpolation_method": "nearest",
            "seawater_resistivity": 0.3,
            "sediment_resistivity": 3.0,
            "crust_resistivity": 3000.0,
            "lithosphere_resistivity": 1000.0,
            "asthenosphere_resistivity": 100.0,
            "transition_zone_resistivity": 10.0,
            "lower_mantle_resistivity": 1.0,
            "transition_zone_top": 410.0,
            "transition_zone_bot": 660.0,
            "profile_max_depth": 1000.0,
            "uri": "http://ds.iris.edu/files/products/emc/emc-files/LITHO1.0.nc",
        }
        if conductivity_params:
            merged_params = {**default_conductivity_params, **dict(conductivity_params)}
        else:
            merged_params = default_conductivity_params.copy()

        self.cprop = RecursiveNamespace(**merged_params)
        self.earth_model = self.cprop.earth_model

        # NOTE these are specified in terms of resistivity, in ohm-m
        self.seawater_resistivity = self.cprop.seawater_resistivity
        self.sediment_resistivity = self.cprop.sediment_resistivity
        self.crust_resistivity = self.cprop.crust_resistivity
        self.lithosphere_resistivity = self.cprop.lithosphere_resistivity
        self.asthenosphere_resistivity = self.cprop.asthenosphere_resistivity
        self.transition_zone_resistivity = self.cprop.transition_zone_resistivity
        self.lower_mantle_resistivity = self.cprop.lower_mantle_resistivity

        # fixed layer depths (below the lithosphere), in km
        self.transition_zone_top = self.cprop.transition_zone_top
        self.transition_zone_bot = self.cprop.transition_zone_bot
        self.profile_max_depth = self.cprop.profile_max_depth

        # things that control how this code works in detail...
        # Options: "nearest" or "linear"
        self.grid_interpolation_method = self.cprop.grid_interpolation_method

        # Load netcdf file
        self.load_earth_model()
        return

    def load_earth_model(self) -> None:
        """
        Load the LITHO1.0 model into cached interpolator callables.

        Raises
        ------
        RuntimeError
            If the Earth model cannot be downloaded or parsed.
        """
        config_dir = Path(".scubas_config")
        config_dir.mkdir(parents=True, exist_ok=True)
        filename = config_dir / self.earth_model
        if not filename.exists():
            uri = self.cprop.uri
            try:
                self._download_earth_model(uri, filename)
            except urllib.error.URLError as exc:
                logger.error(
                    "Unable to download Earth model '%s' from %s", filename, uri
                )
                raise RuntimeError(
                    f"Failed to download Earth model from {uri}"
                ) from exc
            except OSError as exc:
                logger.error("Unable to write Earth model file '%s'", filename)
                raise RuntimeError(
                    f"Failed to persist downloaded Earth model at {filename}"
                ) from exc
        try:
            with netcdf_file(str(filename)) as f:
                latitude = np.copy(f.variables["latitude"][:])
                longitude = np.copy(f.variables["longitude"][:])
                # base of lithosphere (top of asthenosphere)
                asthenospheric_mantle_top_depth = np.copy(
                    f.variables["asthenospheric_mantle_top_depth"][:]
                )
                # top/bottom of mantle lithosphere
                lithosphere_bottom_depth = np.copy(f.variables["lid_bottom_depth"][:])
                lithosphere_top_depth = np.copy(f.variables["lid_top_depth"][:])
                # crustal layers
                lower_crust_bottom_depth = np.copy(
                    f.variables["lower_crust_bottom_depth"][:]
                )
                lower_crust_top_depth = np.copy(f.variables["lower_crust_top_depth"][:])
                middle_crust_bottom_depth = np.copy(
                    f.variables["middle_crust_bottom_depth"][:]
                )
                middle_crust_top_depth = np.copy(
                    f.variables["middle_crust_top_depth"][:]
                )
                upper_crust_bottom_depth = np.copy(
                    f.variables["upper_crust_bottom_depth"][:]
                )
                upper_crust_top_depth = np.copy(f.variables["upper_crust_top_depth"][:])
                # sediment layers
                lower_sediments_bottom_depth = np.copy(
                    f.variables["lower_sediments_bottom_depth"][:]
                )
                lower_sediments_top_depth = np.copy(
                    f.variables["lower_sediments_top_depth"][:]
                )
                middle_sediments_bottom_depth = np.copy(
                    f.variables["middle_sediments_bottom_depth"][:]
                )
                middle_sediments_top_depth = np.copy(
                    f.variables["middle_sediments_top_depth"][:]
                )
                upper_sediments_bottom_depth = np.copy(
                    f.variables["upper_sediments_bottom_depth"][:]
                )
                upper_sediments_top_depth = np.copy(
                    f.variables["upper_sediments_top_depth"][:]
                )
                # water levels
                water_bottom_depth = np.copy(f.variables["water_bottom_depth"][:])
                water_top_depth = np.copy(f.variables["water_top_depth"][:])
        except (IOError, OSError) as exc:
            logger.error("Unable to read Earth model netCDF file '%s'", filename)
            raise RuntimeError(f"Failed to load Earth model from {filename}") from exc
        self.lithosphere_model = {
            "latitude": latitude,
            "longitude": longitude,
            "asthenospheric_mantle_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                asthenospheric_mantle_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lithosphere_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                lithosphere_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lithosphere_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                lithosphere_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lower_crust_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                lower_crust_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lower_crust_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                lower_crust_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "middle_crust_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                middle_crust_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "middle_crust_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                middle_crust_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "upper_crust_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                upper_crust_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "upper_crust_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                upper_crust_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lower_sediments_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                lower_sediments_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "lower_sediments_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                lower_sediments_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "middle_sediments_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                middle_sediments_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "middle_sediments_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                middle_sediments_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "upper_sediments_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                upper_sediments_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "upper_sediments_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                upper_sediments_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "water_bottom_depth": RegularGridInterpolator(
                (latitude, longitude),
                water_bottom_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
            "water_top_depth": RegularGridInterpolator(
                (latitude, longitude),
                water_top_depth,
                method=self.grid_interpolation_method,
                bounds_error=True,
            ),
        }
        return

    @staticmethod
    def _download_earth_model(uri: str, destination: Path) -> None:
        """
        Download and persist the Earth model file.

        Parameters
        ----------
        uri :
            Remote URI to the LITHO1.0 netCDF file.
        destination :
            Filesystem path where the downloaded file is written.
        """
        request = urllib.request.Request(uri)
        with urllib.request.urlopen(request) as response, destination.open("wb") as fh:
            shutil.copyfileobj(response, fh)

    def get_interpolation_points(
        self, pt0: Sequence[float], pt1: Sequence[float]
    ) -> np.ndarray:
        """
        Construct a sequence of latitude/longitude interpolation points.

        Parameters
        ----------
        pt0 :
            Starting coordinate expressed as ``[latitude, longitude]``.
        pt1 :
            Ending coordinate expressed as ``[latitude, longitude]``.

        Returns
        -------
        numpy.ndarray
            Array of points ordered as ``[[lat, lon], ...]`` suitable for
            interpolation against the LITHO1.0 grid.
        """
        globe = Geod(ellps="WGS84")

        # somewhat janky way of figuring out how many interp points to have between
        # the bin edges... basic idea is that we want one point spaced every ~1 deg,
        # since that"s the resolution of LITHO1.0
        npts = int(round(np.sqrt((pt1[1] - pt0[1]) ** 2.0 + (pt1[0] - pt0[0]) ** 2.0)))
        if npts < 1:
            npts = 1

        # obtain equally spaced points along a geodesic (doesn"t include end points)
        latlons = globe.npts(pt0[1], pt0[0], pt1[1], pt1[0], npts)

        # add the bin edges to the lat/lon list... note this array is now [lon, lat]
        interpolation_points = np.vstack(
            (np.array([pt0[1], pt0[0]]), np.array(latlons), np.array([pt1[1], pt1[0]]))
        )

        # now swap back to [lat, lon]
        interpolation_points = np.vstack(
            (interpolation_points[:, 1], interpolation_points[:, 0])
        ).T

        return interpolation_points

    def get_water_layer(
        self, lithosphere_model: Mapping[str, RegularGridInterpolator], pts: np.ndarray
    ) -> float:
        """
        Interpolate and average seawater thickness for the given coordinates.

        Parameters
        ----------
        lithosphere_model :
            Mapping of pre-configured interpolators keyed by layer name.
        pts :
            ``(n, 2)`` array of latitude/longitude points.

        Returns
        -------
        float
            Average thickness of the seawater layer in kilometres.
        """
        water_top_values = lithosphere_model["water_top_depth"](pts)
        water_bot_values = lithosphere_model["water_bottom_depth"](pts)

        water_top = np.nanmean(water_top_values)
        water_bot = np.nanmean(water_bot_values)
        water_thk = water_bot - water_top

        # sanity check... water_top can be slightly different than zero,
        # but shouldn"t be that much different
        if np.absolute(water_top) > 0.01:
            logger.warning(f"PROBLEM: water_top doesn't make sense: {water_top}")

        return water_thk

    def get_sediment_layer(
        self, lithosphere_model: Mapping[str, RegularGridInterpolator], pts: np.ndarray
    ) -> float:
        """
        Interpolate and average sediment thickness for the given coordinates.

        Parameters
        ----------
        lithosphere_model :
            Mapping of pre-configured interpolators keyed by layer name.
        pts :
            ``(n, 2)`` array of latitude/longitude points.

        Returns
        -------
        float
            Average thickness of the sediment layer in kilometres.
        """
        upper_sed_top_values = lithosphere_model["upper_sediments_top_depth"](pts)
        upper_sed_bot_values = lithosphere_model["upper_sediments_bottom_depth"](pts)
        upper_sed_thk_values = upper_sed_bot_values - upper_sed_top_values

        middle_sed_top_values = lithosphere_model["middle_sediments_top_depth"](pts)
        middle_sed_bot_values = lithosphere_model["middle_sediments_bottom_depth"](pts)
        middle_sed_thk_values = middle_sed_bot_values - middle_sed_top_values

        lower_sed_top_values = lithosphere_model["lower_sediments_top_depth"](pts)
        lower_sed_bot_values = lithosphere_model["lower_sediments_bottom_depth"](pts)
        lower_sed_thk_values = lower_sed_bot_values - lower_sed_top_values

        sed_thk_values = np.nansum(
            np.vstack(
                (upper_sed_thk_values, middle_sed_thk_values, lower_sed_thk_values)
            ),
            axis=0,
        )

        sed_thk = np.nanmean(sed_thk_values)

        return sed_thk

    def get_crust_layer(
        self, lithosphere_model: Mapping[str, RegularGridInterpolator], pts: np.ndarray
    ) -> float:
        """
        Interpolate and average crust thickness for the given coordinates.

        Parameters
        ----------
        lithosphere_model :
            Mapping of pre-configured interpolators keyed by layer name.
        pts :
            ``(n, 2)`` array of latitude/longitude points.

        Returns
        -------
        float
            Average thickness of the crustal layer in kilometres.
        """
        upper_crust_top_values = lithosphere_model["upper_crust_top_depth"](pts)
        upper_crust_bot_values = lithosphere_model["upper_crust_bottom_depth"](pts)
        upper_crust_thk_values = upper_crust_bot_values - upper_crust_top_values

        middle_crust_top_values = lithosphere_model["middle_crust_top_depth"](pts)
        middle_crust_bot_values = lithosphere_model["middle_crust_bottom_depth"](pts)
        middle_crust_thk_values = middle_crust_bot_values - middle_crust_top_values

        lower_crust_top_values = lithosphere_model["lower_crust_top_depth"](pts)
        lower_crust_bot_values = lithosphere_model["lower_crust_bottom_depth"](pts)
        lower_crust_thk_values = lower_crust_bot_values - lower_crust_top_values

        crust_thk_values = np.nansum(
            np.vstack(
                (
                    upper_crust_thk_values,
                    middle_crust_thk_values,
                    lower_crust_thk_values,
                )
            ),
            axis=0,
        )

        crust_thk = np.nanmean(crust_thk_values)

        return crust_thk

    def get_lithosphere_layer(
        self, lithosphere_model: Mapping[str, RegularGridInterpolator], pts: np.ndarray
    ) -> float:
        """
        Interpolate and average mantle lithosphere thickness.

        Parameters
        ----------
        lithosphere_model :
            Mapping of pre-configured interpolators keyed by layer name.
        pts :
            ``(n, 2)`` array of latitude/longitude points.

        Returns
        -------
        float
            Average thickness of the lithosphere in kilometres.
        """
        lithosphere_top_values = lithosphere_model["lithosphere_top_depth"](pts)
        lithosphere_bot_values = lithosphere_model["lithosphere_bottom_depth"](pts)
        asthenosphere_top_values = lithosphere_model["asthenospheric_mantle_top_depth"](
            pts
        )

        litho_top = np.nanmean(lithosphere_top_values)
        litho_bot = np.nanmean(lithosphere_bot_values)
        litho_thk = litho_bot - litho_top
        astheno_top = np.nanmean(asthenosphere_top_values)

        # sanity check... make sure top of asthenosphere is the same as bottom of the lithosphere
        if litho_bot != astheno_top:
            logger.warning(
                f"PROBLEM: lithosphere-asthenosphere dont line up: {litho_bot}, {astheno_top}"
            )

        return litho_thk

    def get_upper_mantle_layer(
        self, lithosphere_model: Mapping[str, RegularGridInterpolator], pts: np.ndarray
    ) -> float:
        """
        Interpolate and average upper mantle thickness.

        Parameters
        ----------
        lithosphere_model :
            Mapping of pre-configured interpolators keyed by layer name.
        pts :
            ``(n, 2)`` array of latitude/longitude points.

        Returns
        -------
        float
            Average thickness of the upper mantle in kilometres.
        """
        asthenosphere_top_values = lithosphere_model["asthenospheric_mantle_top_depth"](
            pts
        )
        astheno_top = np.nanmean(asthenosphere_top_values)
        astheno_thk = self.transition_zone_top - astheno_top
        return astheno_thk

    def get_transition_zone_layer(self) -> float:
        """
        Return the fixed mantle transition zone thickness.

        Returns
        -------
        float
            Thickness of the transition zone in kilometres.
        """
        return self.transition_zone_bot - self.transition_zone_top

    def get_lower_mantle_layer(self) -> float:
        """
        Return the fixed lower mantle thickness.

        Returns
        -------
        float
            Thickness of the lower mantle in kilometres.
        """
        return self.profile_max_depth - self.transition_zone_bot

    def _compile_profile_(self, pts: np.ndarray) -> pd.DataFrame:
        """
        Compile a conductivity profile for a set of interpolation points.

        Parameters
        ----------
        pts :
            ``(n, 2)`` array of latitude/longitude points.

        Returns
        -------
        pandas.DataFrame
            Dataframe with ``thickness`` (km), ``resistivity`` (ohm-m), and
            human-readable ``name`` for each layer.
        """
        # now progress through the model layers from near-surface to deep Earth...
        # first the water layer
        water_thk = self.get_water_layer(self.lithosphere_model, pts)
        # sediment layer
        sed_thk = self.get_sediment_layer(self.lithosphere_model, pts)
        # crust layer
        crust_thk = self.get_crust_layer(self.lithosphere_model, pts)
        # mantle lithosphere layer
        litho_thk = self.get_lithosphere_layer(self.lithosphere_model, pts)
        # asthenosphere layer
        astheno_thk = self.get_upper_mantle_layer(self.lithosphere_model, pts)
        # transition zone layer
        tz_thk = self.get_transition_zone_layer()
        # lower mantle layer
        lm_thk = self.get_lower_mantle_layer()
        resistivity_profile = np.array(
            [
                [water_thk, self.seawater_resistivity, "Seawater"],
                [sed_thk, self.sediment_resistivity, "Sediment"],
                [crust_thk, self.crust_resistivity, "Crust"],
                [litho_thk, self.lithosphere_resistivity, "Lithosphere"],
                [astheno_thk, self.asthenosphere_resistivity, "Upper Mantle"],
                [tz_thk, self.transition_zone_resistivity, "Transition Zone"],
                [lm_thk, self.lower_mantle_resistivity, "Lower Mantle"],
            ]
        )
        rf = pd.DataFrame()
        rf["thickness"], rf["resistivity"], rf["name"] = (
            np.array(resistivity_profile[:, 0]).astype(float),
            resistivity_profile[:, 1],
            resistivity_profile[:, 2],
        )
        return rf

    @staticmethod
    def compile_profile(
        latlon: Sequence[float],
        kind: str = "rounded",
        to_site: bool = True,
        site_name: str = "",
        site_description: str = "",
        **conductivity_params: Any,
    ) -> Union[pd.DataFrame, Site]:
        """
        Compile a conductivity profile for a single latitude/longitude pair.

        Parameters
        ----------
        latlon :
            Coordinate pair expressed as ``[latitude, longitude]``.
        kind :
            ``"rounded"`` will round the coordinates to the nearest integer,
            ``"exact"`` keeps supplied precision.
        to_site :
            When ``True`` the return value is a :class:`scubas.datasets.Site`
            populated with conductivity and thickness information.
        site_name :
            Optional site identifier.
        site_description :
            Optional human readable description.
        **conductivity_params :
            Override values passed to :class:`ConductivityProfile`.

        Returns
        -------
        Union[pandas.DataFrame, scubas.datasets.Site]
            SCUBAS site instance or raw conductivity dataframe depending on
            ``to_site``.
        """
        cp = ConductivityProfile(**conductivity_params)
        latlon_array = np.asarray(latlon, dtype=float)
        if kind == "rounded":
            latlon_array = np.rint(latlon_array)
        logger.info(f"Lat-lon: {latlon_array}")
        profile = cp._compile_profile_(latlon_array)
        logger.info(f"Compiled Profile \n {profile}")
        if to_site:
            profile = Site.init(
                1.0 / profile["resistivity"].to_numpy(dtype=float),
                profile["thickness"].to_numpy(dtype=float),
                profile["name"],
                site_description,
                site_name,
            )
        return profile

    @staticmethod
    def compile_profiles(
        latlons: Sequence[Sequence[float]],
        kind: str = "rounded",
        to_site: bool = True,
        site_names: Optional[Sequence[str]] = None,
        site_descriptions: Optional[Sequence[str]] = None,
        **conductivity_params: Any,
    ) -> List[Union[pd.DataFrame, Site]]:
        """
        Compile conductivity profiles for multiple coordinates.

        Parameters
        ----------
        latlons :
            Iterable of ``[latitude, longitude]`` coordinate pairs.
        kind :
            ``"rounded"`` will round the coordinates to the nearest integer,
            ``"exact"`` keeps supplied precision.
        to_site :
            When ``True`` output contains :class:`scubas.datasets.Site` objects.
        site_names :
            Optional sequence of site names corresponding to ``latlons``.
        site_descriptions :
            Optional sequence of site descriptions corresponding to ``latlons``.
        **conductivity_params :
            Override values passed to :class:`ConductivityProfile`.

        Returns
        -------
        list
            List of conductivity profiles as dataframes or ``Site`` instances.
        """
        cp = ConductivityProfile(**conductivity_params)
        profiles = []
        resolved_site_names = list(site_names or [])
        resolved_site_descriptions = list(site_descriptions or [])
        for i, latlon in enumerate(latlons):
            latlon_array = np.asarray(latlon, dtype=float)
            if kind == "rounded":
                latlon_array = np.rint(latlon_array)
            logger.info(f"Lat-lon: {latlon_array}")
            profile = cp._compile_profile_(latlon_array)
            logger.info(f"Compiled Profile \n {profile}")
            if to_site:
                site_name = (
                    resolved_site_names[i] if i < len(resolved_site_names) else ""
                )
                site_description = (
                    resolved_site_descriptions[i]
                    if i < len(resolved_site_descriptions)
                    else ""
                )
                profile = Site.init(
                    1.0 / profile["resistivity"].to_numpy(dtype=float),
                    profile["thickness"].to_numpy(dtype=float),
                    profile["name"],
                    site_description,
                    site_name,
                )
            profiles.append(profile)
        return profiles

    @staticmethod
    def compile_bined_profile(
        bined_latlon: Sequence[Sequence[float]],
        to_site: bool = True,
        site_name: str = "",
        site_description: str = "",
        **conductivity_params: Any,
    ) -> Union[pd.DataFrame, Site]:
        """
        Compile a conductivity profile for a pair of binned coordinates.

        Parameters
        ----------
        bined_latlon :
            Sequence with two coordinate pairs ``[[lat0, lon0], [lat1, lon1]]``.
        to_site :
            When ``True`` return value is a :class:`scubas.datasets.Site`.
        site_name :
            Optional site identifier.
        site_description :
            Optional human readable description.
        **conductivity_params :
            Override values passed to :class:`ConductivityProfile`.

        Returns
        -------
        Union[pandas.DataFrame, scubas.datasets.Site]
            Conductivity profile as dataframe or ``Site`` instance.
        """
        cp = ConductivityProfile(**conductivity_params)
        start = np.asarray(bined_latlon[0], dtype=float)
        end = np.asarray(bined_latlon[1], dtype=float)
        ipts = cp.get_interpolation_points(start, end)
        profile = cp._compile_profile_(ipts)
        logger.info(f"Compiled Profile \n {profile}")
        if to_site:
            profile = Site.init(
                1.0 / profile["resistivity"].to_numpy(dtype=float),
                profile["thickness"].to_numpy(dtype=float),
                profile["name"],
                site_description,
                site_name,
            )
        return profile

    @staticmethod
    def compile_bined_profiles(
        bined_latlons: np.ndarray,
        to_site: bool = True,
        site_names: Optional[Sequence[str]] = None,
        site_descriptions: Optional[Sequence[str]] = None,
        **conductivity_params: Any,
    ) -> List[Union[pd.DataFrame, Site]]:
        """
        Compile conductivity profiles for an ordered set of bin coordinates.

        Parameters
        ----------
        bined_latlons :
            ``(n, 2)`` array of coordinates defining bin edges.
        to_site :
            When ``True`` output contains :class:`scubas.datasets.Site` objects.
        site_names :
            Optional sequence of site names corresponding to each bin.
        site_descriptions :
            Optional sequence of site descriptions corresponding to each bin.
        **conductivity_params :
            Override values passed to :class:`ConductivityProfile`.

        Returns
        -------
        list
            List of conductivity profiles as dataframes or ``Site`` instances.
        """
        cp = ConductivityProfile(**conductivity_params)
        profiles = []
        bined_latlons_array = np.asarray(bined_latlons, dtype=float)
        nbins = len(bined_latlons_array) - 1
        resolved_site_names = list(site_names or [])
        resolved_site_descriptions = list(site_descriptions or [])
        for i in range(nbins):
            ipts = cp.get_interpolation_points(
                bined_latlons_array[i, :], bined_latlons_array[i + 1, :]
            )
            profile = cp._compile_profile_(ipts)
            profile.thickness = profile.thickness * 1e3
            logger.info(f"Compiled Profile \n {profile}")
            if to_site:
                site_name = (
                    resolved_site_names[i] if i < len(resolved_site_names) else ""
                )
                site_description = (
                    resolved_site_descriptions[i]
                    if i < len(resolved_site_descriptions)
                    else ""
                )
                profile = Site.init(
                    1.0 / profile["resistivity"].to_numpy(dtype=float),
                    profile["thickness"].to_numpy(dtype=float),
                    profile["name"],
                    site_description,
                    site_name,
                )
            profiles.append(profile)
        return profiles

    @staticmethod
    def compile_mcmc_bined_profiles(
        bined_latlons: Sequence[Sequence[float]],
        n: int = 100,
        continental_shelves_thickness_range: Optional[Sequence[float]] = None,
        to_site: bool = True,
        site_names: Optional[Sequence[str]] = None,
        site_descriptions: Optional[Sequence[str]] = None,
        random_seed: int = 0,
        **conductivity_params: Any,
    ) -> List[List[Union[pd.DataFrame, Site]]]:
        """
        Compile Monte Carlo conductivity profiles for a set of binned coordinates.

        Parameters
        ----------
        bined_latlons :
            Iterable of ``[latitude, longitude]`` coordinates defining bin edges.
        n :
            Number of Monte Carlo realisations to produce.
        continental_shelves_thickness_range :
            Optional ``[min, max]`` range for perturbing the seawater layer at
            the end bins (kilometres).
        to_site :
            When ``True`` output contains :class:`scubas.datasets.Site` objects.
        site_names :
            Optional sequence of site names corresponding to each bin.
        site_descriptions :
            Optional sequence of site descriptions corresponding to each bin.
        random_seed :
            Random seed used for reproducibility.
        **conductivity_params :
            Override values passed to :class:`ConductivityProfile`.

        Returns
        -------
        list
            Nested list of Monte Carlo profiles, each entry containing the
            profiles for one draw across all bins.
        """
        np.random.seed(random_seed)
        cp = ConductivityProfile(**conductivity_params)
        bined_latlons_array = np.asarray(bined_latlons, dtype=float)
        thickness_range = (
            list(continental_shelves_thickness_range)
            if continental_shelves_thickness_range is not None
            else [0.01, 1.0]
        )
        resolved_site_names = list(site_names or [])
        resolved_site_descriptions = list(site_descriptions or [])
        mcmc_profiles, profiles = [], []
        for i, latlon in enumerate(bined_latlons_array):
            profile = cp._compile_profile_(latlon)
            profile.fillna(0, inplace=True)
            profiles.append(profile)
        for _ in range(n):
            mc_profiles = []
            for i in range(len(profiles) - 1):
                profile = profiles[i].copy()
                profile.thickness = (
                    np.random.uniform(
                        np.array(profiles[i].thickness),
                        np.array(profiles[i + 1].thickness),
                    )
                    * 1e3
                )
                if (i == 0) or (i == len(profiles) - 2):
                    profile.thickness[0] = np.random.uniform(
                        thickness_range[0],
                        thickness_range[1],
                    )
                logger.info(f"Compiled Profile \n {profile}")
                if to_site:
                    site_name = (
                        resolved_site_names[i] if i < len(resolved_site_names) else ""
                    )
                    site_description = (
                        resolved_site_descriptions[i]
                        if i < len(resolved_site_descriptions)
                        else ""
                    )
                    profile = Site.init(
                        1.0 / profile["resistivity"].to_numpy(dtype=float),
                        profile["thickness"].to_numpy(dtype=float),
                        profile["name"],
                        site_description,
                        site_name,
                    )
                mc_profiles.append(profile)
            mcmc_profiles.append(mc_profiles)
        return mcmc_profiles
