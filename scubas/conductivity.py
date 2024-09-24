#!/usr/bin/env python

"""
    conductiviy.py: Module is used to implement Earth conductivity methods
"""

__author__ = "Murphy, B.; Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os

import numpy as np
import pandas as pd
from loguru import logger
from pyproj import Geod
from scipy.interpolate import RegularGridInterpolator
from scipy.io import netcdf_file

from scubas.datasets import Site
from scubas.utils import RecursiveNamespace


class ConductivityProfile(object):
    """
    Class is dedicated to create conductivity profiles
    from LITHO1.0 model.

    1.  LITHO1.0 as the basis for constructing the resistivity profiles
        http://ds.iris.edu/ds/products/emc-litho10/
    """

    def __init__(
        self,
        conductivity_params=dict(
            earth_model="LITHO1.0.nc",
            grid_interpolation_method="nearest",
            seawater_resistivity=0.3,
            sediment_resistivity=3.0,
            crust_resistivity=3000.0,
            lithosphere_resistivity=1000.0,
            asthenosphere_resistivity=100.0,
            transition_zone_resistivity=10.0,
            lower_mantle_resistivity=1.0,
            transition_zone_top=410.0,
            transition_zone_bot=660.0,
            profile_max_depth=1000.0,
            uri="http://ds.iris.edu/files/products/emc/emc-files/LITHO1.0.nc",
        ),
    ):
        self.cprop = RecursiveNamespace(**conductivity_params)
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

    def load_earth_model(self):
        """
        Returns a dict containing the pertinent parts of the Earth model
        Dict entries are callable interpolation objects
        """
        os.makedirs(".scubas_config/", exist_ok=True)
        filename = ".scubas_config/" + self.earth_model
        if not os.path.exists(filename):
            uri = self.cprop.uri
            os.system(f"wget -O {filename} {uri}")
        with netcdf_file(filename) as f:
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
            middle_crust_top_depth = np.copy(f.variables["middle_crust_top_depth"][:])
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

    def get_interpolation_points(self, pt0, pt1):
        """ """
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

    def get_water_layer(self, lithosphere_model, pts):
        """
        Interpolates and averages seawater layer thickness within the given bin
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

    def get_sediment_layer(self, lithosphere_model, pts):
        """
        Interpolates and averages sediment layer thickness within the given bin
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

    def get_crust_layer(self, lithosphere_model, pts):
        """
        Interpolates and averages crust layer thickness within the given bin
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

    def get_lithosphere_layer(self, lithosphere_model, pts):
        """
        Interpolates and averages mantle lithosphere thickness within the given bin
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

    def get_upper_mantle_layer(self, lithosphere_model, pts):
        """
        Interpolates and averages upper mantle thickness within the given bin
        """
        asthenosphere_top_values = lithosphere_model["asthenospheric_mantle_top_depth"](
            pts
        )
        astheno_top = np.nanmean(asthenosphere_top_values)
        astheno_thk = self.transition_zone_top - astheno_top
        return astheno_thk

    def get_transition_zone_layer(self):
        """
        Return fixed mantle transition zone thickness
        """
        return self.transition_zone_bot - self.transition_zone_top

    def get_lower_mantle_layer(self):
        """
        Return fixed lower mantle thickness
        """
        return self.profile_max_depth - self.transition_zone_bot

    def _compile_profile_(self, pts):
        """
        Compile profile for a list of points
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
        latlon,
        kind="rounded",
        to_site=True,
        site_name="",
        site_description="",
        **conductivity_params,
    ):
        """
        Compile profiles for a set of latlons
        """
        cp = ConductivityProfile(**conductivity_params)
        if kind == "rounded":
            latlon = np.rint(latlon)
        logger.info(f"Lat-lon: {latlon}")
        profile = self._compile_profile_(latlon)
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
        latlons,
        kind="rounded",
        to_site=True,
        site_names=[],
        site_descriptions=[],
        **conductivity_params,
    ):
        """
        Compile profiles for a set of latlons
        """

        cp = ConductivityProfile(**conductivity_params)
        profiles = []
        for i, latlon in enumerate(latlons):
            if kind == "rounded":
                latlon = np.rint(latlon)
            logger.info(f"Lat-lon: {latlon}")
            profile = cp._compile_profile_(latlon)
            logger.info(f"Compiled Profile \n {profile}")
            if to_site:
                site_name = site_names[i] if i < len(site_names) else ""
                site_description = (
                    site_descriptions[i] if i < len(site_descriptions) else ""
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
        bined_latlon,
        to_site=True,
        site_name="",
        site_description="",
        **conductivity_params,
    ):
        """
        Compile profiles for a binned latlons
        """
        cp = ConductivityProfile(**conductivity_params)
        ipts = cp.get_interpolation_points(bined_latlon[0], bined_latlon[1])
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
        bined_latlons,
        to_site=True,
        site_names=[],
        site_descriptions=[],
        **conductivity_params,
    ):
        """
        Compile profiles for a binned latlons
        """
        cp = ConductivityProfile(**conductivity_params)
        profiles = []
        nbins = len(bined_latlons) - 1
        for i in range(nbins):
            ipts = cp.get_interpolation_points(
                bined_latlons[i, :], bined_latlons[i + 1, :]
            )
            profile = cp._compile_profile_(ipts)
            logger.info(f"Compiled Profile \n {profile}")
            if to_site:
                site_name = site_names[i] if i < len(site_names) else ""
                site_description = (
                    site_descriptions[i] if i < len(site_descriptions) else ""
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
    def compile_bined_profiles(
        bined_latlons,
        n=100,
        continental_shelves_thickness_range=[0.01, 1.0],
        to_site=True,
        site_names=[],
        site_descriptions=[],
        random_seed=0,
        **conductivity_params,
    ):
        """
        Compile MCMC profiles for a binned latlons
        """
        np.random.seed(random_seed)
        cp = ConductivityProfile(**conductivity_params)
        mcmc_profiles, profiles = [], []
        for i, latlon in enumerate(bined_latlons):
            profile = cp._compile_profile_(latlon)
            profile.fillna(0, inplace=True)
            profiles.append(profile)
        for _ in range(n):
            mc_profiles = []
            for i in range(len(profiles) - 1):
                profile = profiles[i].copy()
                profile.thickness = np.random.uniform(
                    np.array(profiles[i].thickness), np.array(profiles[i + 1].thickness)
                )
                if (i == 0) or (i == len(profiles) - 2):
                    profile.thickness[0] = np.random.uniform(
                        continental_shelves_thickness_range[0],
                        continental_shelves_thickness_range[1]
                    )
                logger.info(f"Compiled Profile \n {profile}")
                if to_site:
                    site_name = site_names[i] if i < len(site_names) else ""
                    site_description = (
                        site_descriptions[i] if i < len(site_descriptions) else ""
                    )
                    profile = Site.init(
                        1.0 / profile["resistivity"].to_numpy(dtype=float),
                        profile["thickness"].to_numpy(dtype=float) * 1e3,
                        profile["name"],
                        site_description,
                        site_name,
                    )
                mc_profiles.append(profile)
            mcmc_profiles.append(mc_profiles)
        return mcmc_profiles
