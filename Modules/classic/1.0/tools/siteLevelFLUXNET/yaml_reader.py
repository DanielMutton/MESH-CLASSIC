import sys
import yaml
with open(sys.argv[1], "r") as f:
  d = yaml.safe_load(f)
print(d[sys.argv[2]])
