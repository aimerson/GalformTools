#! /usr/bin/env python                                                                                                                                                                                                      
class GalformError(Exception):
    """Base class for exceptions in this module."""
    pass


class FileError(GalformError):
    def __init__(self, message, errors=None):
        # Call the base class constructor with the parameters it needs
        super(FileError, self).__init__(message)
