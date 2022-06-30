rule MESA_output:
    input:
        "src/data/MESA_output.tar.gz"
    output: # needs to be one thing only -- directory allows to bypass this limit, but the cache will be a tarball
        directory("src/data/MESA_output/")
    conda:
        "environment.yml"
    cache:
        True
    script:
        "src/scripts/convert_MESA_output.py"
