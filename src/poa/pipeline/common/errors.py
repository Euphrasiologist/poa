"""Errors a pipeline stage raises when it can't do its job.

A stage's library function raises StageError for an expected failure (a
missing input, a tool that failed); its command-line wrapper prints it as
"[err] ..." and exits 1, and `poa run` reports it as that stage failing.
Anything else that escapes a stage is a bug.
"""


class StageError(Exception):
    pass
