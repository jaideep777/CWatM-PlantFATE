from datetime import datetime

import numpy as np
import pandas as pd

from pypfate import Patch as patch
from pypfate import Clim

# Module will be moved to plantFATE package eventually and called from there

class PFatePatch:
    def __init__(self, param_file, acclim_forcing_file):
        self.patch = patch(str(param_file))
        self.time_unit_base = self.process_time_units()
        self.tcurrent = 0
        self.time = []
        self.swp = []
        self.swc = []
        self.trans = []
        self.acclimation_forcing = self.read_acclimation_file(acclim_forcing_file)

    def read_acclimation_file(self, file):
        df = pd.read_csv(file)
        alldates = df['date'].map(lambda x: datetime.strptime(x, "%Y/%m/%d") - self.time_unit_base)
        alldates = alldates.map(lambda x: x.days - 1)
        df['date_jul'] = alldates
        return df

    def runstep(
        self,
        soil_water_potential,
        vapour_pressure_deficit,
        photosynthetic_photon_flux_density,
        temperature,
    ):
        self.patch.update_climate(
            368.9, # co2 - need to make it better
            temperature,
            vapour_pressure_deficit * 1000,
            photosynthetic_photon_flux_density,
            soil_water_potential,
        )
        index_acclim = self.acclimation_forcing.index[self.acclimation_forcing['date_jul'] == self.tcurrent].tolist()
        self.patch.update_climate_acclim(self.tcurrent,
                                         368.9,
                                         self.acclimation_forcing.loc[index_acclim, 'temp.C.'],
                                         self.acclimation_forcing.loc[index_acclim, 'vpd'],
                                         self.calculate_photosynthetic_photon_flux_density(self.acclimation_forcing.loc[index_acclim, 'shortwave.W.m2.']),
                                         soil_water_potential)
        self.patch.simulate_to(self.tcurrent)
        trans = self.patch.props.fluxes.trans / 365
        # return evapotranspiration, soil_specific_depletion_1, soil_specific_depletion_2, soil_specific_depletion_3
        return trans, 0, 0, 0


    def first_step(
        self,
        tstart,
        soil_moisture_layer_1,  # ratio [0-1]
        soil_moisture_layer_2,  # ratio [0-1]
        soil_moisture_layer_3,  # ratio [0-1]
        soil_tickness_layer_1,  # m
        soil_tickness_layer_2,  # m
        soil_tickness_layer_3,  # m
        soil_moisture_wilting_point_1,  # ratio [0-1]print(type(time_unit_base))
        soil_moisture_wilting_point_2,  # ratio [0-1]
        soil_moisture_wilting_point_3,  # ratio [0-1]
        soil_moisture_field_capacity_1,  # ratio [0-1]
        soil_moisture_field_capacity_2,  # ratio [0-1]
        soil_moisture_field_capacity_3,  # ratio [0-1]
        temperature,  # degrees Celcius, mean temperature
        relative_humidity,  # percentage [0-100]
        shortwave_radiation,  # W/m2, daily mean
    ):
        (
            soil_water_potential,
            vapour_pressure_deficit,
            photosynthetic_photon_flux_density,
            temperature,
        ) = self.get_plantFATE_input(
            soil_moisture_layer_1,  # ratio [0-1]
            soil_moisture_layer_2,  # ratio [0-1]
            soil_moisture_layer_3,  # ratio [0-1]
            soil_tickness_layer_1,  # m
            soil_tickness_layer_2,  # m
            soil_tickness_layer_3,  # m
            soil_moisture_wilting_point_1,  # ratio [0-1]
            soil_moisture_wilting_point_2,  # ratio [0-1]
            soil_moisture_wilting_point_3,  # ratio [0-1]
            soil_moisture_field_capacity_1,  # ratio [0-1]
            soil_moisture_field_capacity_2,  # ratio [0-1]
            soil_moisture_field_capacity_3,  # ratio [0-1]
            temperature,  # degrees Celcius, mean temperature
            relative_humidity,  # percentage [0-100]
            shortwave_radiation)  # W/m2, daily mean

        newclim = Clim()
        newclim.tc = temperature #- 273.15  # C
        newclim.ppfd_max = photosynthetic_photon_flux_density * 4
        newclim.ppfd = photosynthetic_photon_flux_density
        newclim.vpd = vapour_pressure_deficit * 1000  # kPa -> Pa
        newclim.swp = soil_water_potential  # MPa


        # Convert time to proper units according to the time unit base
        datestart = datetime(tstart.year, tstart.month, tstart.day)
        datediff = datestart - self.time_unit_base
        datediff = datediff.days - 1



        self.swp.append(soil_water_potential)
        self.swc.append(soil_moisture_layer_2)
        self.trans.append(0)
        self.time.append(datestart)
        self.patch.init(datediff, datediff + 1000)
        self.tcurrent = datediff
        self.patch.update_climate(368.9, temperature, vapour_pressure_deficit * 1000, photosynthetic_photon_flux_density, soil_water_potential)
        # TODO: change to be using panda input
        index_acclim = self.acclimation_forcing.index[self.acclimation_forcing['date_jul'] == self.tcurrent].tolist()
        self.patch.update_climate_acclim(self.tcurrent,
                                         368.9,
                                         self.acclimation_forcing.loc[index_acclim, 'temp.C.'],
                                         self.acclimation_forcing.loc[index_acclim, 'vpd'],
                                         self.calculate_photosynthetic_photon_flux_density(self.acclimation_forcing.loc[index_acclim, 'shortwave.W.m2.']),
                                         soil_water_potential)

    def process_time_units(self):
        time_unit = self.patch.config.time_unit
        time_unit = time_unit.split()
        if time_unit[0] != 'days' or time_unit[1] != 'since':
            print("wrong plantfate unit; cwatm supports only daily timescale")
            return
        else:
            time_unit = time_unit[2].split("-")
            return datetime(int(time_unit[0]),
                            int(time_unit[1]),
                            int(time_unit[2]))

    # def read_acclimation_file(self):
    #     acclimation_file = self.patch.config.

    def simulate(self):
        self.patch.simulate()

    def calculate_soil_water_potential_MPa(
        self,
        soil_moisture,  # [m]
        soil_moisture_wilting_point,  # [m]
        soil_moisture_field_capacity,  # [m]
        soil_tickness,  # [m]
        wilting_point=-1500,  # kPa
        field_capacity=-33,  # kPa
    ):
        # https://doi.org/10.1016/B978-0-12-374460-9.00007-X (eq. 7.16)
        soil_moisture_fraction = soil_moisture / soil_tickness
        assert soil_moisture_fraction >= 0 and soil_moisture_fraction <= 1
        del soil_moisture
        soil_moisture_wilting_point_fraction = (
            soil_moisture_wilting_point / soil_tickness
        )
        assert (
            soil_moisture_wilting_point_fraction >= 0
            and soil_moisture_wilting_point_fraction <= 1
        )
        del soil_moisture_wilting_point
        soil_moisture_field_capacity_fraction = (
            soil_moisture_field_capacity / soil_tickness
        )
        assert (
            soil_moisture_field_capacity_fraction >= 0
            and soil_moisture_field_capacity_fraction <= 1
        )
        del soil_moisture_field_capacity

        n_potential = -(
            np.log(wilting_point / field_capacity)
            / np.log(
                soil_moisture_wilting_point_fraction
                / soil_moisture_field_capacity_fraction
            )
        )
        assert n_potential >= 0
        a_potential = (
            1.5 * 10**6 * soil_moisture_wilting_point_fraction**n_potential
        )
        assert a_potential >= 0
        soil_water_potential = -a_potential * soil_moisture_fraction ** (-n_potential)
        return soil_water_potential / 1_000_000  # Pa to MPa

    def calculate_vapour_pressure_deficit_kPa(self, temperature, relative_humidity):
        assert (
            temperature < 100
        )  # temperature is in Celsius. So on earth should be well below 100.
        assert (
            temperature > -100
        )  # temperature is in Celsius. So on earth should be well above -100.
        assert (
            relative_humidity >= 1 and relative_humidity <= 100
        )  # below 1 is so rare that it shouldn't be there at the resolutions of current climate models, and this catches errors with relative_humidity as a ratio [0-1].
        # https://soilwater.github.io/pynotes-agriscience/notebooks/vapor_pressure_deficit.html
        saturated_vapour_pressure = 0.611 * np.exp(
            (17.502 * temperature) / (temperature + 240.97)
        )  # kPa
        actual_vapour_pressure = (
            saturated_vapour_pressure * relative_humidity / 100
        )  # kPa
        vapour_pressure_deficit = saturated_vapour_pressure - actual_vapour_pressure
        return vapour_pressure_deficit

    def calculate_photosynthetic_photon_flux_density(self, shortwave_radiation, xi=0.5):
        # https://search.r-project.org/CRAN/refmans/bigleaf/html/Rg.to.PPFD.html
        photosynthetically_active_radiation = shortwave_radiation * xi
        photosynthetic_photon_flux_density = (
            photosynthetically_active_radiation * 4.6
        )  #  W/m2 -> umol/m2/s
        return photosynthetic_photon_flux_density

    def get_plantFATE_input(
        self,
        soil_moisture_layer_1,  # m
        soil_moisture_layer_2,  # m
        soil_moisture_layer_3,  # m
        soil_tickness_layer_1,  # m
        soil_tickness_layer_2,  # m
        soil_tickness_layer_3,  # m
        soil_moisture_wilting_point_1,  # m
        soil_moisture_wilting_point_2,  # m
        soil_moisture_wilting_point_3,  # m
        soil_moisture_field_capacity_1,  # m
        soil_moisture_field_capacity_2,  # m
        soil_moisture_field_capacity_3,  # m
        temperature,  # degrees Celcius, mean temperature
        relative_humidity,  # percentage [0-100]
        shortwave_radiation,  # W/m2, daily mean
    ):
        assert (
            temperature < 100
        )  # temperature is in Celsius. So on earth should be well below 100.
        assert relative_humidity >= 0 and relative_humidity <= 100

        soil_water_potential = self.calculate_soil_water_potential_MPa(
            soil_moisture_layer_1 + soil_moisture_layer_2 + soil_moisture_layer_3,
            soil_moisture_wilting_point_1
            + soil_moisture_wilting_point_2
            + soil_moisture_wilting_point_3,
            soil_moisture_field_capacity_1
            + soil_moisture_field_capacity_2
            + soil_moisture_field_capacity_3,
            soil_tickness_layer_1 + soil_tickness_layer_2 + soil_tickness_layer_3,
        )

        vapour_pressure_deficit = self.calculate_vapour_pressure_deficit_kPa(
            temperature, relative_humidity
        )

        photosynthetic_photon_flux_density = (
            self.calculate_photosynthetic_photon_flux_density(shortwave_radiation)
        )

        return (
            soil_water_potential,
            vapour_pressure_deficit,
            photosynthetic_photon_flux_density,
            temperature,
        )

    def step(
        self,
        curr_time,
        soil_moisture_layer_1,  # ratio [0-1]
        soil_moisture_layer_2,  # ratio [0-1]
        soil_moisture_layer_3,  # ratio [0-1]
        soil_tickness_layer_1,  # m
        soil_tickness_layer_2,  # m
        soil_tickness_layer_3,  # m
        soil_moisture_wilting_point_1,  # ratio [0-1]
        soil_moisture_wilting_point_2,  # ratio [0-1]
        soil_moisture_wilting_point_3,  # ratio [0-1]
        soil_moisture_field_capacity_1,  # ratio [0-1]
        soil_moisture_field_capacity_2,  # ratio [0-1]
        soil_moisture_field_capacity_3,  # ratio [0-1]
        temperature,  # degrees Celcius, mean temperature
        relative_humidity,  # percentage [0-100]
        shortwave_radiation,  # W/m2, daily mean
    ):
        (
            soil_water_potential,
            vapour_pressure_deficit,
            photosynthetic_photon_flux_density,
            temperature,
        ) = self.get_plantFATE_input(
            soil_moisture_layer_1,  # ratio [0-1]
            soil_moisture_layer_2,  # ratio [0-1]
            soil_moisture_layer_3,  # ratio [0-1]
            soil_tickness_layer_1,  # m
            soil_tickness_layer_2,  # m
            soil_tickness_layer_3,  # m
            soil_moisture_wilting_point_1,  # ratio [0-1]
            soil_moisture_wilting_point_2,  # ratio [0-1]
            soil_moisture_wilting_point_3,  # ratio [0-1]
            soil_moisture_field_capacity_1,  # ratio [0-1]
            soil_moisture_field_capacity_2,  # ratio [0-1]
            soil_moisture_field_capacity_3,  # ratio [0-1]
            temperature,  # degrees Celcius, mean temperature
            relative_humidity,  # percentage [0-100]
            shortwave_radiation,  # W/m2, daily mean
        )

        curr_time_dt = datetime(curr_time.year, curr_time.month, curr_time.day)
        timediff = curr_time_dt - self.time_unit_base
        self.tcurrent = timediff.days - 1


        (
            evapotranspiration,
            soil_specific_depletion_1,
            soil_specific_depletion_2,
            soil_specific_depletion_3,
        ) = self.runstep(
            soil_water_potential,
            vapour_pressure_deficit,
            photosynthetic_photon_flux_density,
            temperature
        )

        self.swp.append(soil_water_potential)
        self.swc.append(soil_moisture_layer_2)
        self.time.append(curr_time_dt)
        self.trans.append(evapotranspiration)

        soil_specific_depletion_1 = (
            np.nan
        )  # this is currently not calculated in plantFATE, so just setting to np.nan to avoid confusion
        soil_specific_depletion_2 = (
            np.nan
        )  # this is currently not calculated in plantFATE, so just setting to np.nan to avoid confusion
        soil_specific_depletion_3 = (
            np.nan
        )  # this is currently not calculated in plantFATE, so just setting to np.nan to avoid confusion

        if curr_time.month == 12 and curr_time.day == 31:
            df = pd.DataFrame(data={'Time': self.time,
                                    'SWP': self.swp,
                                    'SWC': self.swc,
                                    'Trans': self.trans})
            # write the DataFrame to a CSV file
            df.to_csv('Soil_Water.csv')

        evapotranspiration = evapotranspiration / 1000  # kg H2O/m2/day to m/day

        return (
            evapotranspiration,
            soil_specific_depletion_1,
            soil_specific_depletion_2,
            soil_specific_depletion_3,
        )

    def finalize(self):

        # create a Pandas DataFrame from the dictionary
        df = pd.DataFrame(data={'Time': range(1, len(self.swp)),
                'SWP': self.swp})
        print(df)
        # write the DataFrame to a CSV file
        df.to_csv('SWP.csv')
        print(self.swp)
        print("closing patch")
        self.patch.close()

    @property
    def n_individuals(self):
        return sum(self.patch.props.structure.n_ind_vec)

    @property
    def biomass(self):
        return sum(self.patch.cwm.biomass_vec)  # kgC / m2
