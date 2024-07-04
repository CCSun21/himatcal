from maggma.stores import MongoURIStore
from monty.serialization import loadfn
from himatcal import SETTINGS
from datetime import datetime
from monty.json import jsanitize

def save_to_db(label, info={}, database='himat',collection_name="job"):
    store = MongoURIStore(
        uri=SETTINGS.MONGODB_URI,
        database=database,
        collection_name=collection_name,
    )

    base_info = {
        "label": label,
        "time": jsanitize(datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
    }

    document = {**base_info, **info}
    with store:
        store.update([document])
    print(f"{label} has been saved to the database!")


def load_from_db(label, database='himat',collection_name="job"):
    store = MongoURIStore(
        uri=SETTINGS.MONGODB_URI,
        database=database,
        collection_name=collection_name,
    )

    with store:
        documents = store.query(criteria={"label": label})
    if len(documents) == 0:
        print(f"{label} not found in the database!")
        return None
    else:
        print(f"{label} has been loaded from the database!")
        return documents