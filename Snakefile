rule try:
  input:
    "src/data/history.data"
  output:
    "src/data/try.npy"
  conda:
    "environment.yml"
  cache:
    True
  script:
    "src/scripts/convert_MESA_output.py"
