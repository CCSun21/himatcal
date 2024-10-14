from __future__ import annotations

import logging
from datetime import datetime, timezone

from maggma.stores import MongoURIStore
from monty.json import jsanitize

from himatcal import SETTINGS


def save_to_db(label, info=None, database="himat", collection_name="job"):
    """
    Save a document to the database.
    """
    if info is None:
        info = {}
    store = MongoURIStore(
        uri=SETTINGS.MONGODB_URI,
        database=database,
        collection_name=collection_name,
    )

    base_info = {
        "label": label,
        "time": jsanitize(datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")),
    }

    document = {**base_info, **info}
    with store:
        store.update([document])
    logging.debug(f"{label} has been saved to the database!")


def load_from_db(label, database="himat", collection_name="job"):
    """
    Load a document from the database.
    """
    store = MongoURIStore(
        uri=SETTINGS.MONGODB_URI,
        database=database,
        collection_name=collection_name,
    )

    with store:
        documents = store.query(criteria={"label": label})
    if len(documents) == 0:
        logging.debug(f"{label} not found in the database!")
        return None
    else:
        logging.debug(f"{label} has been loaded from the database!")
        return documents
