import json
import os

from jsonschema import validate
import pytest

from sourmash.schema import SBT_SCHEMAS, MINHASH_SCHEMAS


def test_validate_data():
    dirpath = os.path.join(os.path.dirname(__file__), 'test-data')

    for root, dirs, files in os.walk(dirpath):
        for test_data_file in files:
            try:
                with open(os.path.join(root, test_data_file), 'r') as f:
                    instance = json.load(f)
            except json.decoder.JSONDecodeError:
                print(f"skipped {test_data_file}: Not JSON")
                continue

            if test_data_file.endswith(".sig"):
                validate(instance=instance, schema=MINHASH_SCHEMAS["0.4"])
                print(f"validated {test_data_file}")
            elif test_data_file.endswith(".sbt.json"):
                try:
                    validate(instance=instance, schema=SBT_SCHEMAS[instance["version"]])
                    print(f"validated {test_data_file}")
                except KeyError:
                    print(f"need to write SBT {instance['version']} schema for {test_data_file}")
            else:
                #raise Exception("Couldn't decide schema to use")
                print(f"skipped {test_data_file}")
