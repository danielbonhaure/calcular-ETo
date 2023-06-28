#!/usr/bin/env python

"""
Documentation:
¯¯¯¯¯¯¯¯¯¯¯¯¯
 • FAO Irrigation and Drainage Paper No. 56 (Fao56.pdf)


Download links:
¯¯¯¯¯¯¯¯¯¯¯¯¯¯
 • https://www.researchgate.net/publication/284300773_FAO_Irrigation_and_drainage_paper_No_56        --> (Fao56.pdf)


Glossary (* source data, • derived data):
¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
 * tmax    -> maximum temperature [°C]            externally provided, MANDATORY
 * tmin    -> minimum temperature [°C]            externally provided, MANDATORY
 • tmean   -> mean daily air temperature [°C]                                   --> pag: 33 equation: 19 (Fao56.pdf)
 • Delta   -> slope vapour pressure curve [kPa °C-1]                            --> pag: 37 equation: 13 (Fao56.pdf)
 • P       -> atmospheric pressure [kPa]                                        --> pag: 31 equation: 7  (Fao56.pdf)
 • gamma   -> psychrometric constant [kPa °C-1]                                 --> pag: 32 equation: 8  (Fao56.pdf)
 • eotmin  -> saturation vapour pressure at the air temperature T [kPa]         --> pag: 36 equation: 11 (Fao56.pdf)
 • eotmax  -> saturation vapour pressure at the air temperature T [kPa]         --> pag: 36 equation: 11 (Fao56.pdf)
 • es      -> saturation vapour pressure [kPa]                                  --> pag: 36 equation: 12 (Fao56.pdf)
 * hr      -> relative humidity [%]               externally provided, SUGGESTED
 * tdew    -> dewpoint temperature [°C]           externally provided, SUGGESTED
       obs -> it can be filled (or even estimated) from tmin                    --> pag: 261 equation: 6.6 (Fao56.pdf)
 • ea      -> actual vapour pressure [kPa]
       op1 -> with dew point                                                    --> pag: 37 equation: 14 (Fao56.pdf)
       op2 -> with relative humidity (RH)                                       --> pag: 39 equation: 19 (Fao56.pdf)
 • delta   -> solar declination [rad]                                           --> pag: 46 equation: 24 (Fao56.pdf)
 • fi      -> latitude [rad]                                                    --> pag: 46 equation: 22 (Fao56.pdf)
 • omega   -> sunset hour angle [rad]                                           --> pag: 46 equation: 25 (Fao56.pdf)
 • dr      -> inverse relative distance Earth-Sun                               --> pag: 46 equation: 23 (Fao56.pdf)
 • Ra      -> extraterrestrial radiation for daily periods [MJ m-2day-1]        --> pag: 46 equation: 21 (Fao56.pdf)
 • N       -> maximum possible duration of sunshine or daylight hours [hour]    --> pag: 48 equation: 34 (Fao56.pdf)
 * helio   -> heliophany [hours/day]              externally provided, SUGGESTED
 * nub     -> cloudiness or nubosity [%]          externally provided, SUGGESTED
 • Rs      -> solar or shortwave radiation [MJ m-2 day-1]
       op1 -> with average sunshine hours (heliophany), n/N = helio/N           --> pag: 50 equation: 35 (Fao56.pdf)
       op2 -> with cloudiness (nubosity), n/N = 1 - (nub/8)                     --> pag: 50 equation: 35 (Fao56.pdf)
       op3 -> with tmin and tmax                                                --> pag: 60 equation: 50 (Fao56.pdf)
 • Rso     -> calculated clear-sky radiation [MJ m-2 day-1]                     --> pag: 51 equation: 36 (Fao56.pdf)
 • Rns     -> net solar or shortwave radiation [MJ m-2 day-1]                   --> pag: 51 equation: 38 (Fao56.pdf)
 • Rnl     -> net outgoing long-wave radiation [MJ m-2 day-1]                   --> pag: 52 equation: 39 (Fao56.pdf)
 • Rn      -> net radiation at the crop surface [MJ m-2 day-1]                  --> pag: 53 equation: 40 (Fao56.pdf)
 * vmed    -> average wind speed [m/s]            externally provided, SUGGESTED
 • u2      -> wind speed at 2 m height [m/s]                                    --> pag: 56 equation: 47 (Fao56.pdf)
 • G       -> soil heat flux [MJ m-2 day-1]                                     --> pag: 54 equation: 42 (Fao56.pdf)
 • Rad     -> left term of the numerator of the equation6, at page 65
           -> 0.408 * Delta * (Rn - G)                                          --> pag: 54 equation: 42 (Fao56.pdf)
 • Aero    -> right term of the numerator of the equation 6, at page 65
           -> gamma * (900 / (T + 273)) * u2 * (es - ea), T=tmean               --> pag: 54 equation: 42 (Fao56.pdf)
 • Denom   -> the denominator of the equation 6, at page 65
           -> Delta + gamma * (1 + 0.34 * u2)                                   --> pag: 54 equation: 42 (Fao56.pdf)
 • eto     -> reference evapotranspiration [mm day-1]                           --> pag: 65 equation: 6  (Fao56.pdf)
       obs -> (Rad+Aero)/Denom, reformulation of the equation 6, at page 65


MANDATORY DATA, externally provided:
¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
 • tmax           -> It can be filled externally (see: https://github.com/danielbonhaure/imputar-faltantes)
 • tmin           -> It can be filled externally (see: https://github.com/danielbonhaure/imputar-faltantes)


SUGGESTED DATA, externally provided:
¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
 • hr or tdew     -> Tdew can be estimated with the equation 6.6, page 261 (Fao56.pdf)
 • helio or nub   -> Only required for compute Rs, but Rs can be estimated also using tmax and tmin
 • vmed           -> It can be filled with the monthly average or 0.5

"""

import logging
import argparse
import numpy as np
import pandas as pd

from pathlib import Path
from argparse import Namespace
from dataclasses import dataclass


# Conf logging
logging.basicConfig(format='%(asctime)s -- %(levelname)7s -- %(message)s',
                    datefmt='%Y/%m/%d %I:%M:%S %p', level=logging.INFO)


@dataclass
class WeatherStation:
    lat_dd: float
    lon_dd: float
    elev_m: int


def calcular_eto(input_csv: Path, output_csv: Path, station_data: WeatherStation) -> None:
    """ Python script to calculate ETo. """

    # Read input data
    df = pd.read_csv(input_csv.absolute(), parse_dates=["date"]).set_index("date", drop=False)

    # Check MANDATORY data
    assert 'tmax' in df.columns
    assert not df.tmax.hasnans
    assert 'tmin' in df.columns
    assert not df.tmin.hasnans

    # Check SUGGESTED data
    if 'hr' not in df.columns:
        df['hr'] = np.nan
    if 'tdew' not in df.columns:
        df['tdew'] = np.nan
    if 'nub' not in df.columns:
        df['nub'] = np.nan
    if 'helio' not in df.columns:
        df['helio'] = np.nan

    # Calculate ETo
    df['tmean'] = (df.tmin + df.tmax) / 2
    df['Delta'] = (4098 * (0.6108 * np.exp((17.27 * df.tmean) / (df.tmean + 237.3)))) / np.power((df.tmean + 237.3), 2)
    df['elev'] = station_data.elev_m
    df['P'] = 101.3 * np.power(((293 - (0.0065 * df.elev)) / 293), 5.26)
    df['gamma'] = 0.000665 * df.P
    df['eotmin'] = 0.6108 * np.exp((17.27 * df.tmin) / (df.tmin + 237.3))
    df['eotmax'] = 0.6108 * np.exp((17.27 * df.tmax) / (df.tmax + 237.3))
    df['es'] = (df.eotmin + df.eotmax) / 2
    df['hr'] = df.hr                          # externally provided, there is no method to estimate missing values
    df['tdew'] = df.tdew.fillna(df.tmin - 2)  # externally provided, filled with equation 6.6, page 261 (Fao56.pdf)
    df['ea'] = np.where(df.hr.notna(), (df.hr / 100) * df.es, 0.6108 * np.exp((17.27 * df.tdew) / (df.tdew + 237.3)))
    df['es-ea'] = df.es - df.ea
    df['jday'] = df.index.dayofyear
    df['delta'] = 0.409 * np.sin((((2 * np.pi) / 365) * df.jday) - 1.39)
    df['latitude'] = station_data.lat_dd
    df['longitude'] = station_data.lon_dd
    df['fi'] = (np.pi / 180) * df.latitude
    df['omega'] = np.arccos(-np.tan(df.fi) * np.tan(df.delta))
    df['dr'] = 1 + 0.033 * np.cos(((2 * np.pi) / 365) * df.jday)
    df['Ra'] = ((24 * 60) / np.pi) * 0.0820 * df.dr * (
            (df.omega * np.sin(df.fi) * np.sin(df.delta)) + (np.cos(df.fi) * np.cos(df.delta) * np.sin(df.omega)))
    df['N'] = (24 / np.pi) * df.omega
    df['helio'] = df.helio  # externally provided, only required for compute Rs, but Rs can be externally provided
    df['nub'] = df.nub      # externally provided, only required for compute Rs, but Rs can be externally provided
    if 'Rs' not in df.columns:  # Rs is only calculated if not yet available in the input dataframe
        df['Rs'] = np.where(df.helio.notna(), (0.25 + (0.5 * (df.helio / df.N))) * df.Ra,
                            np.where(df.nub.notna(), (0.25 + (0.5 * (1 - (df.nub / 8)))) * df.Ra,
                                     0.16 * np.sqrt(np.absolute(df.tmax - df.tmin)) * df.Ra))
    df['Rso'] = (0.25 + 0.5) * df.Ra
    df['Rns'] = (1 - 0.23) * df.Rs
    df['Rnl'] = 4.903 * pow(10, -9) * ((np.power(df.tmax + 273.16, 4) + np.power(df.tmin + 273.16, 4)) / 2) * (
            0.34 - (0.14 * np.sqrt(df.ea))) * ((1.35 * (df.Rs / df.Rso)) - 0.35)
    df['Rn'] = df.Rns - df.Rnl
    df['vmed'] = df.vmed.fillna(0.5)  # filling method can be improved using monthly average instead of 0.5
    df['u2'] = df.vmed * (4.87 / np.log((67.8 * 10) - 5.42))
    df['G'] = 0
    df['Rad'] = 0.408 * df.Delta * (df.Rn - df.G)
    df['Aero'] = df.gamma * (900 / (df.tmean + 273)) * df.u2 * (df.es - df.ea)
    df['Denom'] = df.Delta + (df.gamma * (1 + (0.34 * df.u2)))
    df['eto'] = (df.Rad + df.Aero) / df.Denom

    # Save all process data
    all_output_csv = str(input_csv).replace(input_csv.suffix, f'_all{input_csv.suffix}')
    df.to_csv(all_output_csv, index=False)

    # Save only ETo (without missing values)
    faoeto = df.eto.asfreq('D')
    faoeto.to_csv(output_csv.absolute(), index=True)


def parse_args() -> Namespace:
    """ Function to parse script arguments. """

    parser = argparse.ArgumentParser(description='Python script to calculate ETo.')
    parser.add_argument('-l', '--station-lon-dd', dest='station_lon_dd',
                        help='Station longitude in decimal degrees (dd)',
                        required=True, type=float)
    parser.add_argument('-a', '--station-lat-dd', dest='station_lat_dd',
                        help='Station latitude in decimal degrees (dd)',
                        required=True, type=float)
    parser.add_argument('-e', '--station-elev-m', dest='station_elev_m',
                        help='Station elevation above sea level in meters (m)',
                        required=True, type=float)
    parser.add_argument('-i', '--input-csv', dest='input_csv', default='station_data.csv',
                        help='Path to a CSV file that contains the weather station\'s meteorological data. '
                             'Mandatory meteorological daily data: tmax, tmin. '
                             'Suggested meteorological daily data: hr or tdew, helio or nub or Rs, vmed.',
                        required=False, type=Path)
    parser.add_argument('-o', '--output-csv', dest='output_csv', default='faoeto.csv',
                        help='Path to the output CSV file that will be created.',
                        required=False, type=Path)
    return parser.parse_args()


if __name__ == '__main__':
    logging.info('Start ETo calculation.')

    # Get script args
    args = parse_args()

    # Set Weather Station Data
    station = WeatherStation(
        lon_dd=args.station_lon_dd,
        lat_dd=args.station_lat_dd,
        elev_m=args.station_elev_m
    )

    try:

        # Compute ETo
        calcular_eto(
            input_csv=args.input_csv,
            output_csv=args.output_csv,
            station_data=station
        )

    except Exception as e:
        logging.error("Exception occurred", exc_info=True)
        raise  # see: http://www.markbetz.net/2014/04/30/re-raising-exceptions-in-python/

    else:
        logging.info('ETo calculation completed successfully.')

