rule try:
     output:
       "src/data/try.npy"
     cache:
        True
     input:
       "src/scripts/convert_MESA_output.py"
