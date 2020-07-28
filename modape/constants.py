"""Constants for modape"""
import re

REGEX_PATTERNS = {
    "product": re.compile(r'M\w{6}'),
    "version": re.compile(r'.+\.(\d{3})\..+'),
    "tile": re.compile(r'h\d+v\d+'),
    "date": re.compile(r'.+A(\d{7}).+'),
    "processing_timestamp": re.compile(r'.+(\d{13}).+')
}
