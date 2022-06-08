rule try:
  input:
    "src/data/MESA_output.tar.gz"
  output:
    "src/data/try.npy"
  conda:
    "environment.yml"
  # cache:
  #   True
  script:
    "src/scripts/convert_MESA_output.py"
