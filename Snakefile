rule MESA_output:
    input:
        "src/data/MESA_output.tar.gz"
    output:
        directory("src/data/MESA_output/")
    conda:
        "environment.yml"
    cache:
        True
    script:
        "src/scripts/convert_MESA_output.py"
