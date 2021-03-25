"""Constants for modape"""
import re

REGEX_PATTERNS = {
    "product": re.compile(r"^M[O|Y]D\d{2}\w{1}\d{1}"),
    "version": re.compile(r".+\.(\d{3})\..+"),
    "tile": re.compile(r"h\d+v\d+"),
    "date": re.compile(r".+A(\d{7}).+"),
    "processing_timestamp": re.compile(r".+(\d{13}).+"),
    "AquaTerra": re.compile(r"^M[Y|O]D\d{2}.+"),
    "Aqua": re.compile(r"^MYD\d{2}.+"),
    "Terra": re.compile(r"^MYD\d{2}.+"),
    "VIM": re.compile(r"^M[Y|O|X]D13.+"),
    "LST": re.compile(r"^M[Y|O]D11.+"),
    "VIMLST": re.compile(r"^M[Y|O]D[11|13].+"),
}

VAM_PRODUCT_CODES = dict(
        zip(["VIM", "VEM", "LTD", "LTN"],
            ["NDVI", "EVI", "LST_Day", "LST_Night"])
)

TEMPORAL_DICT = {
    "MXD13":{
        "temporalresolution": 8,
        "tshift": 8,
    },
    "MOD13":{
        "temporalresolution": 16,
        "tshift": 8,
    },
    "MYD13":{
        "temporalresolution": 16,
        "tshift": 8,
    },
    "MOD11":{
        "temporalresolution": 8,
        "tshift": 4,
    },
    "MYD11":{
        "temporalresolution": 8,
        "tshift": 4,
    },
}

LST_NAME_LUD = {
    "LTD": {
        "MOD": "TDT",
        "MYD": "TDA"
    },
    "LTN": {
        "MOD": "TNT",
        "MYD": "TNA"
    },
}

TEMPINT_LABELS = {
    5: "p",
    10: "d",
}

DATE_LABELS = {
    5: dict(zip([3, 8, 13, 18, 23, 28], ["p1", "p2", "p3", "p4", "p5", "p6"])),
    10: dict(zip([5, 15, 25], ["d1", "d2", "d3"]))
}
