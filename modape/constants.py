"""Constants for modape"""

import re
from typing import Any

from osgeo import gdal

REGEX_PATTERNS = {
    "product": re.compile(r"^(?:VNP|M[O|Y]D)\d{2}\w{1}\d{1}"),
    "version": re.compile(r".+\.(\d{3})\..+"),
    "tile": re.compile(r"h\d+v\d+"),
    "date": re.compile(r".+A(\d{7}).+"),
    "processing_timestamp": re.compile(r".+(\d{13}).+"),
    "AquaTerra": re.compile(r"^M[Y|O]D\d{2}.+"),
    "Aqua": re.compile(r"^MYD\d{2}.+"),
    "Terra": re.compile(r"^MYD\d{2}.+"),
    "VIM": re.compile(r"^(?:VNP|M[O|Y|X]D)13.+"),
    "LST": re.compile(r"^M[Y|O]D11.+"),
    "VIMLST": re.compile(r"^(?:VNP|M[O|Y]D)[11|13].+"),
}

VAM_PRODUCT_CODES = dict(zip(["VIM", "VEM", "LTD", "LTN"], ["NDVI", "EVI", "LST_Day", "LST_Night"]))

PRODUCT_SRS_DICT: dict[str, dict[str, Any]] = {
    "VNP13A2": {
        "ProjectionWKT": """\
PROJCRS["unnamed",
    BASEGEOGCRS["Unknown datum based upon the custom spheroid",
        DATUM["Not specified (based on custom spheroid)",
            ELLIPSOID["Custom spheroid",6371007.181,0,
                LENGTHUNIT["metre",1,
                    ID["EPSG",9001]]]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433,
                ID["EPSG",9122]]]],
    CONVERSION["unnamed",
        METHOD["Sinusoidal"],
        PARAMETER["Longitude of natural origin",0,
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8802]],
        PARAMETER["False easting",0,
            LENGTHUNIT["Meter",1],
            ID["EPSG",8806]],
        PARAMETER["False northing",0,
            LENGTHUNIT["Meter",1],
            ID["EPSG",8807]]],
    CS[Cartesian,2],
        AXIS["easting",east,
            ORDER[1],
            LENGTHUNIT["Meter",1]],
        AXIS["northing",north,
            ORDER[2],
            LENGTHUNIT["Meter",1]]]""",
        "AxisMapping": "1,2",
        "PixelSize": tuple([xy * 926.6254330558334 for xy in [1, -1]]),
        "Origin": {"h": dict(index=18, offset=0), "v": dict(index=9, offset=0)},
    }
}

# Mapping for <NASA Product>_<VAM Subdataset> to physical encoding of the data: value range, nodata
PRODUCT_SDS_DICT: dict[str, dict[str, int | list[int]]] = {
    "VNP13A2_NDVI": {
        "Name": "//HDFEOS/GRIDS/VIIRS_Grid_16Day_VI_1km/Data_Fields/1_km_16_days_NDVI",
        "DataType": gdal.GDT_Int16,
        "ValueRange": (-10000, 10000),
        "NoDataValue": (-15000, -13000),
        "Size": (1200, 1200),
        "BlockSize": (1200, 1),
    }
}

TEMPORAL_DICT = {
    "VNP13": {
        "temporalresolution": 8,
        "tshift": 8,
    },
    "MXD13": {
        "temporalresolution": 8,
        "tshift": 8,
    },
    "MOD13": {
        "temporalresolution": 16,
        "tshift": 8,
    },
    "MYD13": {
        "temporalresolution": 16,
        "tshift": 8,
    },
    "MOD11": {
        "temporalresolution": 8,
        "tshift": 4,
    },
    "MYD11": {
        "temporalresolution": 8,
        "tshift": 4,
    },
}

LST_NAME_LUD = {
    "LTD": {"MOD": "TDT", "MYD": "TDA"},
    "LTN": {"MOD": "TNT", "MYD": "TNA"},
}

TEMPINT_LABELS = {
    5: "p",
    10: "d",
}

DATE_LABELS = {
    5: dict(zip([3, 8, 13, 18, 23, 28], ["p1", "p2", "p3", "p4", "p5", "p6"])),
    10: dict(zip([5, 15, 25], ["d1", "d2", "d3"])),
}
