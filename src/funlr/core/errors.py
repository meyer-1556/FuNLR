"""Shared exceptions for FuNLR."""


class FunlrError(Exception):
    """Base class for all FuNLR errors."""


class ConfigError(FunlrError):
    """Invalid or incomplete configuration."""


class ToolNotFound(FunlrError):
    """A required external executable could not be resolved."""


class StageError(FunlrError):
    """Fatal error during a pipeline stage."""
