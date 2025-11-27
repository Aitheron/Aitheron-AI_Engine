import os
import logging

import sentry_sdk
from sentry_sdk.integrations.logging import LoggingIntegration
from dotenv import load_dotenv

_SENTRY_INITIALIZED = False

load_dotenv()

def _init_sentry():
    global _SENTRY_INITIALIZED

    if _SENTRY_INITIALIZED:
        return

    dsn = os.getenv("SENTRY_DSN")
    if not dsn:
        return

    sentry_logging = LoggingIntegration(
        level=logging.INFO,
        event_level=logging.INFO,
    )

    sentry_sdk.init(
        dsn=dsn,
        integrations=[sentry_logging],
        traces_sample_rate=1.0,
        environment=os.getenv("SENTRY_ENVIRONMENT", "production"),
    )

    _SENTRY_INITIALIZED = True

def get_logger(name: str = "app"):
    _init_sentry()

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    return logger
