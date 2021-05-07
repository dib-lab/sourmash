import json
import pkg_resources


SBT_SCHEMAS = {
  5: json.loads(pkg_resources.resource_string('sourmash', 'schema/sbt.schema')),
  4: json.loads(pkg_resources.resource_string('sourmash', 'schema/sbt.schema'))
}

MINHASH_SCHEMAS = {
  '0.4': json.loads(pkg_resources.resource_string('sourmash', 'schema/minhash.schema'))
}
